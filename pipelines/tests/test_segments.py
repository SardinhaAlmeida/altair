"""Run: python3 -m unittest discover -s pipelines/tests -p test_GetSegments.py -v"""

import csv
import importlib.util
import io
import random
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "segments", Path(__file__).with_name("GetSegments.py"))
SEGMENTS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SEGMENTS)

COLUMNS = ["Accession", "Species", "Isolate", "Segment", "GenBank_Title", "Length",
           "Nuc_Completeness", "Country"]


def mutate(sequence, rate, rng):
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b
                   for b in sequence)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class ClassificationTests(unittest.TestCase):
    def test_segment_is_primary_and_title_is_fallback(self):
        cases = [
            (("A", "virus isolate X, complete sequence"), "A"),
            (("DNA_A", ""), "A"),
            (("DNA B", ""), "B"),
            (("segment A", ""), "A"),
            (("DNA1", ""), "A"),
            (("DNA 2", ""), "B"),
            (("", "virus X DNA, segment: B, complete sequence"), "B"),
            (("", "virus X component A, complete sequence"), "A"),
            (("", "virus X segment DNA2, complete sequence"), "B"),
        ]
        for (segment, title), expected in cases:
            with self.subTest(segment=segment, title=title):
                self.assertEqual(SEGMENTS.classify_component(segment, title)[0], expected)

    def test_unclear_records_stay_unknown_with_a_reason(self):
        cases = [
            (("1", ""), "unrecognised Segment"),
            (("RNA1", ""), "unrecognised Segment"),
            (("A", "virus X segment B, complete sequence"), "Segment says A"),
            (("", "virus X DNA-A and DNA-B"), "both components"),
            (("", "Cassava mosaic alphasatellite DNA-1, complete sequence"), "satellite"),
            (("", "virus X, complete genome"), "no component"),
        ]
        for (segment, title), reason in cases:
            with self.subTest(segment=segment, title=title):
                component, _, note = SEGMENTS.classify_component(segment, title)
                self.assertEqual(component, "UNKNOWN")
                self.assertIn(reason, note)

    def test_title_core_pairs_components_of_one_isolate_only(self):
        core = SEGMENTS.title_pair_core
        self.assertEqual(core("SLCMV CaKT2 DNA, segment: A, complete sequence"),
                         core("SLCMV CaKT2 DNA, segment: B, complete sequence"))
        self.assertNotEqual(core("SLCMV CaKT2 DNA, segment: A, complete sequence"),
                            core("SLCMV CaKT6 DNA, segment: B, complete sequence"))
        self.assertEqual(core("SLCMV strain Erode:2011 segment DNA-A, complete sequence"),
                         core("SLCMV strain Erode:2011B segment DNA-B, complete sequence"))

    def test_kmers_ignore_strand_rotation_and_ambiguous_windows(self):
        rng = random.Random(1)
        sequence = "".join(rng.choice("ACGT") for _ in range(300))
        reference = SEGMENTS.kmer_set(sequence)
        self.assertEqual(SEGMENTS.kmer_set(reverse_complement(sequence)), reference)
        self.assertEqual(SEGMENTS.kmer_set(sequence[100:] + sequence[:100]), reference)
        _, size_with_n = SEGMENTS.kmer_set(sequence[:150] + "N" + sequence[151:])
        self.assertLess(size_with_n, reference[1])


class PipelineTests(unittest.TestCase):
    """Two passes on a synthetic dataset of two A-like and two B-like families."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        rng = random.Random(7)
        base_a = "".join(rng.choice("ACGT") for _ in range(2700))
        base_b = "".join(rng.choice("ACGT") for _ in range(2650))
        self.sequences = {
            "A1": mutate(base_a, 0.04, rng), "A2": mutate(base_a, 0.04, rng),
            "B1": mutate(base_b, 0.04, rng), "B2": mutate(base_b, 0.04, rng),
            "U1": mutate(base_a, 0.04, rng),                       # DNA-A without label
            "U2": reverse_complement(mutate(base_b, 0.04, rng)),   # DNA-B, other strand
            "U3": "".join(rng.choice("ACGT") for _ in range(2700)),  # unrelated
        }
        self.sequences["NC_000001"] = self.sequences["A1"]          # RefSeq copy of A1
        self.rows = [
            ("A1", "iso1", "A", "virus iso1 segment A, complete sequence"),
            ("B1", "iso1", "B", "virus iso1 segment B, complete sequence"),
            ("A2", "iso2", "DNA-A", "virus iso2 DNA-A, complete sequence"),
            ("B2", "iso3", "DNA-B", "virus iso3 DNA-B, complete sequence"),
            ("U1", "iso4", "", "virus iso4, complete genome"),
            ("U2", "iso2", "", "virus iso2, complete genome"),
            ("U3", "iso5", "", "virus iso5, complete genome"),
            ("NC_000001", "iso1", "A", "virus iso1 segment A, complete sequence"),
        ]
        self.write_csv(self.dir / "sequences.csv", COLUMNS, [
            [acc, "Virus x", iso, seg, title, len(self.sequences[acc]), "complete", "X"]
            for acc, iso, seg, title in self.rows])
        with open(self.dir / "all.fasta", "w") as handle:
            for acc, seq in self.sequences.items():
                handle.write(f">{acc}.1 test\n{seq}\n")

    def tearDown(self):
        self.tmp.cleanup()

    @staticmethod
    def write_csv(path, header, rows):
        with open(path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    def run_segments(self, *extra):
        argv = [str(self.dir / "sequences.csv"), *extra]
        with redirect_stdout(io.StringIO()):
            SEGMENTS.run(SEGMENTS.parse_args(argv))

    def read(self, name, out="paired_output"):
        with open(self.dir / out / name, newline="") as handle:
            return list(csv.DictReader(handle))

    def accessions(self, name, out="paired_output"):
        return [row["Accession"] for row in self.read(name, out)]

    def test_first_pass_lists_unknowns_for_download(self):
        self.run_segments()
        self.assertEqual(self.accessions("unknown_component.csv"), ["U1", "U2", "U3"])
        self.assertFalse((self.dir / "paired_output" / "similarity_report.csv").exists())
        # Before removal, A1 and its RefSeq copy share a key with B1: ambiguous.
        self.assertCountEqual(self.accessions("ambiguous.csv"), ["A1", "B1", "NC_000001"])

    def test_second_pass_labels_by_similarity_and_removes_refseq_copy(self):
        self.run_segments("--fasta", str(self.dir / "all.fasta"))
        pairs = self.read("pairs.csv")
        self.assertEqual([(r["Accession"], r["Component"]) for r in pairs],
                         [("A1", "A"), ("B1", "B"), ("A2", "A"), ("U2", "B")])
        self.assertEqual(self.accessions("unpaired_A.csv"), ["U1"])
        self.assertEqual(self.accessions("unpaired_B.csv"), ["B2"])
        self.assertEqual(self.accessions("unknown_component.csv"), ["U3"])
        excluded = self.read("excluded_records.csv")
        self.assertEqual([(r["Accession"], r["Reason"]) for r in excluded],
                         [("NC_000001", "RefSeq copy: identical to A1")])
        notes = {r["Accession"]: r["Note"] for r in self.read("similarity_report.csv")}
        self.assertTrue(notes["U3"].startswith("not labelled"))
        self.assertEqual(notes["A2"], "label agrees with similarity")

    def test_manual_files_override_and_pair(self):
        self.write_csv(self.dir / "manual_components.csv", ["Accession", "Component"],
                       [["U3", "B"]])
        self.write_csv(self.dir / "manual_pairs.csv", ["Accession_A", "Accession_B"],
                       [["U1", "B2"]])
        self.run_segments("--fasta", str(self.dir / "all.fasta"))
        pairs = self.accessions("pairs.csv")
        self.assertEqual(pairs[:2], ["U1", "B2"])           # manual pairs come first
        self.assertIn("U3", self.accessions("unpaired_B.csv"))
        self.assertEqual(self.accessions("unknown_component.csv"), [])

    def test_manual_reference_to_removed_refseq_copy_fails(self):
        self.write_csv(self.dir / "manual_pairs.csv", ["Accession_A", "Accession_B"],
                       [["NC_000001", "B2"]])
        with self.assertRaises(SystemExit) as caught, redirect_stdout(io.StringIO()):
            SEGMENTS.run(SEGMENTS.parse_args(
                [str(self.dir / "sequences.csv"), "--fasta", str(self.dir / "all.fasta")]))
        self.assertIn("RefSeq copy: identical to A1", str(caught.exception))

    def add_record(self, accession, sequence, isolate, segment, title, host="H1"):
        self.sequences[accession] = sequence
        with open(self.dir / "all.fasta", "a") as handle:
            handle.write(f">{accession}.1 test\n{sequence}\n")
        with open(self.dir / "sequences.csv", "a", newline="") as handle:
            csv.writer(handle).writerow([accession, "Virus x", isolate, segment, title,
                                         len(sequence), "complete", "X", host])

    def rewrite_with_hosts(self, hosts):
        """Add a Host column; hosts maps accession -> host (default H1)."""
        with open(self.dir / "sequences.csv", newline="") as handle:
            rows = list(csv.reader(handle))
        rows[0].append("Host")
        for row in rows[1:]:
            row.append(hosts.get(row[0], "H1"))
        self.write_csv(self.dir / "sequences.csv", rows[0], rows[1:])

    def test_quality_filters_run_before_similarity_and_are_reported(self):
        self.rewrite_with_hosts({"B2": "H2", "U3": ""})
        rng = random.Random(3)
        n_rich = self.sequences["A2"][:2000] + "N" * 700
        self.add_record("Q1", n_rich, "iso6", "A", "virus iso6 segment A, complete sequence")
        self.add_record("Q2", self.sequences["A2"][:1400], "iso7", "A",
                        "virus iso7 segment A, complete sequence")
        self.add_record("Q3", mutate(self.sequences["A1"], 0.04, rng), "iso8", "",
                        "UNVERIFIED: virus iso8, complete genome")
        (self.dir / "segments_options.txt").write_text(
            "# selection used for the test\n--exclude-unverified\n--host=h1\n"
            "--max-ambiguous=0.01\n--max-length-deviation=0.10\n")
        self.run_segments("--fasta", str(self.dir / "all.fasta"))
        reasons = {r["Accession"]: r["Reason"] for r in self.read("excluded_records.csv")}
        self.assertEqual(reasons["Q3"], "UNVERIFIED GenBank record")
        self.assertEqual(reasons["B2"], "host not selected: H2")
        self.assertEqual(reasons["U3"], "host not recorded")
        self.assertTrue(reasons["Q1"].startswith("ambiguous bases: 25.9%"))
        self.assertTrue(reasons["Q2"].startswith("length: 1400 bp"))
        self.assertEqual(reasons["NC_000001"], "RefSeq copy: identical to A1")
        reported = {r["Accession"] for r in self.read("similarity_report.csv")}
        self.assertFalse(reported & set(reasons))      # excluded: not labelled, not references
        self.assertEqual(self.accessions("unknown_component.csv"), [])

    def test_first_pass_does_not_download_metadata_exclusions(self):
        self.rewrite_with_hosts({"U3": "H2"})
        self.run_segments("--host", "H1")
        self.assertEqual(self.accessions("unknown_component.csv"), ["U1", "U2"])
        self.assertEqual(self.accessions("excluded_records.csv"), ["U3"])

    def test_excluded_record_in_manual_file_fails(self):
        self.rewrite_with_hosts({"B2": "H2"})
        self.write_csv(self.dir / "manual_pairs.csv", ["Accession_A", "Accession_B"],
                       [["U1", "B2"]])
        with self.assertRaises(SystemExit) as caught, redirect_stdout(io.StringIO()):
            self.run_segments("--host", "H1")
        self.assertIn("B2 (host not selected: H2)", str(caught.exception))

    def test_out_dir_is_respected(self):
        self.run_segments("--out-dir", str(self.dir / "elsewhere"))
        self.assertTrue((self.dir / "elsewhere" / "pairs.csv").is_file())
        self.assertFalse((self.dir / "paired_output").exists())

    def test_missing_required_column_fails(self):
        self.write_csv(self.dir / "sequences.csv", ["Accession", "GenBank_Title"],
                       [["A1", "virus segment A"]])
        with self.assertRaises(SystemExit) as caught, redirect_stdout(io.StringIO()):
            self.run_segments()
        self.assertIn("Segment", str(caught.exception))

    def test_dataset_without_component_labels_is_refused(self):
        self.write_csv(self.dir / "sequences.csv", COLUMNS, [
            ["X1", "Virus x", "i1", "", "virus i1, complete genome", 6000, "complete", "X"]])
        with self.assertRaises(SystemExit) as caught, redirect_stdout(io.StringIO()):
            self.run_segments()
        self.assertIn("bipartite", str(caught.exception))


if __name__ == "__main__":
    unittest.main()