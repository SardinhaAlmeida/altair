"""Run: ALTAIR_TEST_BINARY=/path/to/AltaiR python3 -m unittest discover \
    -s pipelines/tests -p test_filter_mraws.py -v
"""

import csv
import importlib.util
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    "filter_mraws", Path(__file__).with_name("Filter_mRAWs.py"))
FILTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FILTER)
BINARY = os.environ.get("ALTAIR_TEST_BINARY")


class FilterUnitTests(unittest.TestCase):
    def test_minimality_requires_both_factors_and_absent_full_word(self):
        self.assertEqual(FILTER.classify("ACG", {"ACG"}), "retained")
        self.assertEqual(FILTER.classify("ACG", {"ACG", "AC"}),
                         "rejected_prefix_absent")
        self.assertEqual(FILTER.classify("ACG", {"ACG", "CG"}),
                         "rejected_suffix_absent")
        self.assertEqual(FILTER.classify("ACG", {"ACG", "AC", "CG"}),
                         "rejected_both_absent")
        self.assertEqual(FILTER.classify("ACG", set()), "rejected_present_in_host")

    def test_bad_input_and_in_place_output_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            data, host = p / "raw.csv", p / "host.fa"
            host.write_text(">host\nACGT\n")
            for row in ["1,1,ACN,50,3\n", "1,1,ACG,nan,3\n",
                        "1,1,ACG,50,4\n", "1,0,ACG,50,3\n",
                        "genome_idx,position,raw,gc_pct,kmer_len\n"]:
                data.write_text(row)
                with self.assertRaises(ValueError):
                    FILTER.read_candidates(data)
            with self.assertRaisesRegex(ValueError, "distinct"):
                FILTER.filter_candidates(data, data, host, "unused")


@unittest.skipUnless(BINARY, "Set ALTAIR_TEST_BINARY for integration tests")
class AltairIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.host = self.root / "host.fa"
        self.binary = str(Path(BINARY).resolve())

    def write_host(self, records):
        self.host.write_text("".join(
            f">host_{index}_{'x' * 150}\n{sequence}\n"
            for index, sequence in enumerate(records)))

    def run_filter(self, rows, ignore_ir=True):
        input_path, output = self.root / "raw.csv", self.root / "mraw.csv"
        with input_path.open("w", newline="") as handle:
            csv.writer(handle).writerows(rows)
        original = input_path.read_bytes()
        report = FILTER.filter_candidates(input_path, output, self.host,
                                          self.binary, ignore_ir)
        self.assertEqual(input_path.read_bytes(), original)
        self.assertEqual(json.loads(output.with_suffix(".summary.json").read_text()), report)
        with output.open(newline="") as handle:
            actual = list(csv.reader(handle))
        return actual, report

    def test_exhaustive_small_forward_words_against_independent_membership(self):
        sequence = "AAAAAACGTACGTAAAAAA"
        self.write_host([sequence])
        words = ["".join(chars) for k in (2, 3, 4)
                 for chars in itertools.product("ACGT", repeat=k)]
        rows = [["1", str(i + 1), word, "50.00", str(len(word))]
                for i, word in enumerate(words)]
        expected = [row for row in rows if row[2] not in sequence
                    and row[2][:-1] in sequence and row[2][1:] in sequence]
        self.assertTrue(expected)
        # Duplicate occurrences must survive, with coordinates and GC untouched.
        rows.append(["2", "78", expected[0][2], "50.00", expected[0][4]])
        expected.append(rows[-1])
        actual, report = self.run_filter(rows)
        self.assertEqual(actual, expected)
        self.assertEqual(report["retained_occurrences"], len(expected))

    def test_membership_matches_direct_altair_with_ambiguities_and_boundaries(self):
        self.write_host(["AAAAACGNNNTACGTTAAAA", "CCCCCTGCAACCCCC"])
        target_sequence = "ACGTACGTTGCAACGATCGATTTTGGGGCCCCAAAAACGT"
        target = self.root / "target.fa"
        target.write_text(f">target_{'x' * 150}\n{target_sequence}\n")
        for ignore_ir in (True, False):
            for k in (2, 3, 4):
                command = [self.binary, "raw", "-f", "-min", str(k), "-max", str(k)]
                if ignore_ir:
                    command.append("-i")
                subprocess.run(command + [str(self.host), str(target)], cwd=self.root,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                               check=True)
                direct = {line.split("\t")[2].strip() for line in
                          Path(f"{target}-k{k}.eg").read_text().splitlines()}
                words = {target_sequence[i:i+k]
                         for i in range(len(target_sequence) - k)}
                queried = FILTER.absent_queries(words, self.host, self.binary,
                                                ignore_ir, self.root)
                self.assertEqual(queried, direct)

    def test_empty_input_and_no_retained_rows(self):
        self.write_host(["AAAAAACGTACGTAAAAAA"])
        for rows in ([], [["1", "1", "ACG", "66.67", "3"]]):
            actual, report = self.run_filter(rows)
            self.assertEqual(actual, [])
            self.assertEqual(report["retained_occurrences"], 0)

    def test_failed_executable_does_not_replace_existing_output(self):
        self.write_host(["AAAAAACGTACGTAAAAAA"])
        input_path, output = self.root / "raw.csv", self.root / "mraw.csv"
        input_path.write_text("1,1,ACG,66.67,3\n")
        output.write_text("previous result\n")
        with self.assertRaises(RuntimeError):
            FILTER.filter_candidates(input_path, output, self.host, "/bin/false")
        self.assertEqual(output.read_text(), "previous result\n")


if __name__ == "__main__":
    unittest.main()
