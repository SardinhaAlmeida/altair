#!/usr/bin/env python3
"""Prepare unpadded DNA components and a layout table for display coordinates."""
import argparse
import csv
import re
from pathlib import Path

WRAP = 80
MIN_GAP = 500
STEP = 250
LAYOUT_FIELDS = ("record_idx", "genome_idx", "component", "accession", "length", "b_offset")


def read_fasta(path):
    records = {}
    header = None
    seq_chunks = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    accession = extract_accession(header)
                    records[accession] = {
                        "header": header,
                        "seq": "".join(seq_chunks).upper().replace(" ", "")
                    }

                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

        if header is not None:
            accession = extract_accession(header)
            records[accession] = {
                "header": header,
                "seq": "".join(seq_chunks).upper().replace(" ", "")
            }

    return records

def normalize_accession(accession):
    accession = (accession or "").strip()
    return re.sub(r"\.\d+$", "", accession)

def extract_accession(header):
    return normalize_accession(header.split()[0])

def classify_component(segment, title):
    text = f"{segment} {title}".lower()

    if re.search(r"\b(dna[\s\-]?a|component a|segment a|dna1)\b", text):
        return "A"

    if re.search(r"\b(dna[\s\-]?b|component b|segment b|dna2)\b", text):
        return "B"

    return "UNKNOWN"


def normalize_title(title):
    x = (title or "").strip()
    x = re.sub(r"\s+", " ", x)
    return x


def normalize_common_title(title_a="", title_b=""):
    """
    Build a common title by removing explicit A/B component wording.
    """
    x = normalize_title(title_a or title_b)

    x = re.sub(r"(?i)\bcomponent\s+[ab]\b", "", x)
    x = re.sub(r"(?i)\bsegment\s+dna[\s\-]?[ab]\b", "", x)
    x = re.sub(r"(?i)\bsegment\s+[ab]\b", "", x)
    x = re.sub(r"(?i)\bdna[\s\-]?[ab]\b", "", x)
    x = re.sub(r"(?i)\bdna[12]\b", "", x)
    x = re.sub(r"(?i),?\s*complete\s+(sequence|genome)\b", "", x)
    x = re.sub(r"\s+", " ", x).strip(" ,;-")

    return x


def compute_b_start(max_a_len, step=STEP):
    """Return the display offset (B position 1 is offset + 1), not a base count to write."""
    return ((max_a_len + MIN_GAP + step - 1) // step) * step


def wrap_seq(seq, width=WRAP):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def read_csv_rows(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def get_accession(row):
    return normalize_accession(row.get("Accession"))

def get_title(row):
    return row.get("GenBank_Title", "") or ""


def get_segment(row):
    return row.get("Segment", "") or ""


def get_component(row):
    """Use the Component column written by GetSegments.py; classify older CSVs from metadata."""
    component = (row.get("Component") or "").strip().upper()
    if component in ("A", "B"):
        return component
    return classify_component(get_segment(row), get_title(row))


def load_pairs_csv(path):
    """
    Load complete A/B pairs.

    This keeps the strict behavior of the original JoinSegments.py:
    pairs.csv must be ordered as:
      A row
      B row
      A row
      B row
      ...
    """
    rows = read_csv_rows(path)

    if len(rows) % 2 != 0:
        raise ValueError("pairs.csv has an odd number of data rows. It must be A/B ordered in pairs.")

    pairs = []

    for i in range(0, len(rows), 2):
        row_a = rows[i]
        row_b = rows[i + 1]

        comp_a = get_component(row_a)
        comp_b = get_component(row_b)

        if comp_a != "A" or comp_b != "B":
            raise ValueError(
                f"Rows {i + 2} and {i + 3} in pairs.csv are not ordered as A then B.\n"
                f"  Row {i + 2}: Accession={get_accession(row_a)} Component={comp_a}\n"
                f"  Row {i + 3}: Accession={get_accession(row_b)} Component={comp_b}"
            )

        pairs.append((row_a, row_b))

    return pairs


def load_unpaired_csv(path, expected_component):
    """
    Load unpaired DNA-A or DNA-B rows.

    The CSV should have at least:
      Accession
      Segment and/or GenBank_Title

    If a row cannot be classified from Segment/GenBank_Title, it is still accepted
    as the expected component because the file itself is unpaired_A.csv or unpaired_B.csv.
    """
    rows = read_csv_rows(path)
    accepted = []
    warnings = []

    for csv_line, row in enumerate(rows, start=2):
        accession = get_accession(row)
        component = get_component(row)

        if not accession:
            warnings.append(f"{path.name}: line {csv_line} ignored because Accession is empty.")
            continue

        if component not in (expected_component, "UNKNOWN"):
            warnings.append(
                f"{path.name}: line {csv_line} has Accession={accession}, "
                f"but looks like component {component}; expected {expected_component}. "
                f"It was ignored."
            )
            continue

        accepted.append(row)

    return accepted, warnings


def prepare_components(pairs_csv, unpaired_a_csv, unpaired_b_csv,
                       input_fasta, output_fasta, metadata_path=None):
    input_fasta, output_fasta = Path(input_fasta), Path(output_fasta)
    metadata_path = Path(metadata_path) if metadata_path else output_fasta.with_suffix(".layout.tsv")
    inputs = [Path(pairs_csv), Path(unpaired_a_csv), Path(unpaired_b_csv), input_fasta]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(f"Input not found: {path}")
    if output_fasta.resolve() == metadata_path.resolve() or any(
        output.resolve() == path.resolve()
        for output in (output_fasta, metadata_path) for path in inputs
    ):
        raise ValueError("FASTA, layout and input paths must be distinct")

    pairs = load_pairs_csv(Path(pairs_csv))
    unpaired_a, warnings_a = load_unpaired_csv(Path(unpaired_a_csv), "A")
    unpaired_b, warnings_b = load_unpaired_csv(Path(unpaired_b_csv), "B")
    records = read_fasta(input_fasta)
    genomes = []
    missing = []

    # Keep an available component even when its expected partner was not downloaded.
    for row_a, row_b in pairs + [(row, None) for row in unpaired_a] + [(None, row) for row in unpaired_b]:
        components = []
        for component, row in (("A", row_a), ("B", row_b)):
            if row is None:
                continue
            accession = get_accession(row)
            record = records.get(accession)
            if record is None or not record["seq"]:
                missing.append(accession)
                continue
            components.append((component, accession, record["seq"]))
        if components:
            title = normalize_common_title(get_title(row_a or {}), get_title(row_b or {}))
            genomes.append((title, components))

    if not genomes:
        raise ValueError("No usable DNA-A or DNA-B records found in FASTA.")

    max_a_len = max((len(seq) for _, parts in genomes
                     for component, _, seq in parts if component == "A"), default=0)
    b_offset = compute_b_start(max_a_len)
    layout = []
    with output_fasta.open("w", encoding="utf-8") as fasta:
        for genome_idx, (title, components) in enumerate(genomes, start=1):
            status = "PAIR" if len(components) == 2 else f"{components[0][0]}_ONLY"
            for component, accession, seq in components:
                record_idx = len(layout) + 1
                fasta.write(f">{accession} | genome_idx={genome_idx} | component={component} | "
                            f"status={status} | {title}\n{wrap_seq(seq)}\n")
                layout.append(dict(zip(LAYOUT_FIELDS, (
                    record_idx, genome_idx, component, accession, len(seq), b_offset))))
    with metadata_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LAYOUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(layout)

    print(f"Wrote {len(layout)} unpadded components for {len(genomes)} genomes: {output_fasta}")
    print(f"Layout metadata: {metadata_path}")
    print(f"Common DNA-B first display position: {b_offset + 1}")
    print(f"Visual gap after longest DNA-A: {b_offset - max_a_len} positions")
    if missing:
        print("Missing or empty accessions: " + ", ".join(sorted(set(missing))))
    for warning in warnings_a + warnings_b:
        print(f"Warning: {warning}")
    return layout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pairs_csv", type=Path)
    parser.add_argument("unpaired_a_csv", type=Path)
    parser.add_argument("unpaired_b_csv", type=Path)
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("--metadata", type=Path,
                        help="Output layout TSV (default: output FASTA with .layout.tsv suffix)")
    args = parser.parse_args()
    prepare_components(args.pairs_csv, args.unpaired_a_csv, args.unpaired_b_csv,
                       args.input_fasta, args.output_fasta, args.metadata)


if __name__ == "__main__":
    main()