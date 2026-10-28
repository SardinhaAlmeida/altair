"""Regression checks for provenance, generated inputs, and failure reporting."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


HELPER = Path(__file__).with_name("run_summary.py")
PIPELINE = Path(__file__).with_name("pipeline.sh")


class RunSummaryTests(unittest.TestCase):
    def invoke(self, *arguments):
        subprocess.run([sys.executable, str(HELPER), *map(str, arguments)],
                       check=True, capture_output=True, text=True)

    def test_analysis_fingerprints_inputs_and_counts_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            host, pathogen = root / "host.fa", root / "pathogen.fa"
            host.write_text(">host\nACGT\n")
            pathogen.write_text(">pathogen\nAACG\n")
            summary = root / "results" / "run_summary.json"
            self.invoke("start", summary, "--workflow", "analyse", "--name", "example",
                        "--host", host, "--pathogen", pathogen,
                        "--kmin", 2, "--kmax", 3, "--started-at", "2026-10-18T12:00:00Z")
            data = json.loads(summary.read_text())
            self.assertEqual(data["status"], "running")
            self.assertEqual(data["inputs"]["host_fasta"]["sha256"],
                             hashlib.sha256(host.read_bytes()).hexdigest())
            (summary.parent / "pathogen.fa-k2.eg").write_text("0\t1\tAA\n")
            (summary.parent / "pathogen.fa-k3.eg").touch()
            (summary.parent / "merged_gc.csv").write_text("0,1,AA,0,2\n")
            (summary.parent / "merged_mraw_gc.csv").touch()
            self.invoke("finish", summary, "--exit-code", 0, "--elapsed-seconds", 7)
            data = json.loads(summary.read_text())
            self.assertEqual(data["status"], "completed")
            self.assertEqual(data["elapsed_seconds"], 7)
            self.assertEqual(data["output_counts"]["raw_occurrences_by_k"], {"2": 1, "3": 0})
            self.assertEqual(data["output_counts"]["missing_raw_k"], [])
            self.assertEqual(data["output_counts"]["raw_occurrences"], 1)
            self.assertEqual(data["output_counts"]["mraw_occurrences"], 0)

    def test_generated_fasta_is_recorded_separately_from_original_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pathogen = root / "example_segments.fasta"
            pathogen.write_text(">old\nAAAA\n")
            summary = root / "paired_output" / "example_run_summary.json"
            self.invoke("start", summary, "--workflow", "prepare", "--name", "example",
                        "--pathogen", pathogen, "--started-at", "2026-10-18T12:00:00Z")
            original_hash = hashlib.sha256(pathogen.read_bytes()).hexdigest()
            pathogen.write_text(">new\nACGT\n")
            (summary.parent / "pairs.csv").write_text("Accession_A,Accession_B\nA,B\n")
            self.invoke("finish", summary, "--exit-code", 1, "--elapsed-seconds", 2)
            data = json.loads(summary.read_text())
            self.assertEqual(data["status"], "failed")
            self.assertEqual(data["exit_code"], 1)
            self.assertEqual(data["inputs"]["pathogen_fasta"]["sha256"], original_hash)
            self.assertEqual(data["prepared_artifacts"]["pathogen_fasta"]["sha256"],
                             hashlib.sha256(pathogen.read_bytes()).hexdigest())
            self.assertEqual(data["output_counts"]["pairs.csv"], 1)

    def test_pipeline_records_failure_and_preserves_exit_code(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = root / "dataset"
            dataset.mkdir()
            host = root / "host.fa"
            host.write_text(">host\nACGT\n")
            environment = dict(os.environ)
            environment["PATH"] = f"{Path(sys.executable).parent}:{environment['PATH']}"
            result = subprocess.run(
                ["bash", str(PIPELINE), "complete", str(host), str(dataset),
                 "example", "private@example.com", "2", "3", "secret-api-key"],
                env=environment, capture_output=True, text=True)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn("metadata CSV not found", result.stderr)
            summary = root / "dataset" / "results" / "example" / "run_summary.json"
            data = json.loads(summary.read_text())
            self.assertEqual(data["inputs"]["host_fasta"]["path"], str(host))
            self.assertEqual(data["status"], "failed")
            self.assertEqual(data["exit_code"], result.returncode)
            self.assertNotIn("private@example.com", summary.read_text())
            self.assertNotIn("secret-api-key", summary.read_text())


if __name__ == "__main__":
    unittest.main()
