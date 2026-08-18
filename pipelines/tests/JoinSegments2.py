#!/usr/bin/env python3
import csv
import math
import re
import sys
from pathlib import Path

WRAP = 80
PAD_BASE = "A"
STEP = 500


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
    """
    DNA-B starts at the next common multiple of STEP,
    always at least STEP bases after the longest DNA-A.
    """
    return math.ceil((max_a_len + step) / step) * step


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

        comp_a = classify_component(get_segment(row_a), get_title(row_a))
        comp_b = classify_component(get_segment(row_b), get_title(row_b))

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
        component = classify_component(get_segment(row), get_title(row))

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


def main():
    if len(sys.argv) != 6:
        print(
            "Usage:\n"
            "  python3 JoinSegments2.py pairs.csv unpaired_A.csv unpaired_B.csv input_sequences.fasta output_joined.fasta"
        )
        sys.exit(1)

    pairs_csv = Path(sys.argv[1])
    unpaired_a_csv = Path(sys.argv[2])
    unpaired_b_csv = Path(sys.argv[3])
    input_fasta = Path(sys.argv[4])
    output_fasta = Path(sys.argv[5])

    for path, label in [
        (pairs_csv, "pairs.csv"),
        (unpaired_a_csv, "unpaired_A.csv"),
        (unpaired_b_csv, "unpaired_B.csv"),
        (input_fasta, "input FASTA"),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    pairs = load_pairs_csv(pairs_csv)
    unpaired_a_rows, warnings_a = load_unpaired_csv(unpaired_a_csv, "A")
    unpaired_b_rows, warnings_b = load_unpaired_csv(unpaired_b_csv, "B")
    warnings = warnings_a + warnings_b

    fasta_records = read_fasta(input_fasta)

    missing_accessions = []
    usable_a_lengths = []
    usable_b_lengths = []

    for row_a, row_b in pairs:
        acc_a = get_accession(row_a)
        acc_b = get_accession(row_b)

        if acc_a in fasta_records:
            usable_a_lengths.append(len(fasta_records[acc_a]["seq"]))
        else:
            missing_accessions.append(acc_a)

        if acc_b in fasta_records:
            usable_b_lengths.append(len(fasta_records[acc_b]["seq"]))
        else:
            missing_accessions.append(acc_b)

    for row_a in unpaired_a_rows:
        acc_a = get_accession(row_a)

        if acc_a in fasta_records:
            usable_a_lengths.append(len(fasta_records[acc_a]["seq"]))
        else:
            missing_accessions.append(acc_a)

    for row_b in unpaired_b_rows:
        acc_b = get_accession(row_b)

        if acc_b in fasta_records:
            usable_b_lengths.append(len(fasta_records[acc_b]["seq"]))
        else:
            missing_accessions.append(acc_b)

    if not usable_a_lengths and not usable_b_lengths:
        raise ValueError("No usable DNA-A or DNA-B records found in FASTA.")

    max_a_len = max(usable_a_lengths) if usable_a_lengths else 0
    max_b_len = max(usable_b_lengths) if usable_b_lengths else 0

    b_start = compute_b_start(max_a_len, STEP)
    common_total_len = b_start + max_b_len

    joined_entries = []

    complete_pairs_written = 0
    unpaired_a_written = 0
    unpaired_b_written = 0
    skipped_pairs = 0
    skipped_unpaired_a = 0
    skipped_unpaired_b = 0

    # Complete pairs:
    # DNA-A + A-gap until B_start + DNA-B + final A-padding.
    for row_a, row_b in pairs:
        acc_a = get_accession(row_a)
        acc_b = get_accession(row_b)

        if acc_a not in fasta_records or acc_b not in fasta_records:
            skipped_pairs += 1
            continue

        seq_a = fasta_records[acc_a]["seq"]
        seq_b = fasta_records[acc_b]["seq"]

        gap_len = b_start - len(seq_a)
        right_pad_len = common_total_len - (len(seq_a) + gap_len + len(seq_b))

        if gap_len < 0 or right_pad_len < 0:
            raise ValueError(f"Invalid padding for pair {acc_a} / {acc_b}.")

        joined_seq = seq_a + (PAD_BASE * gap_len) + seq_b + (PAD_BASE * right_pad_len)

        common_title = normalize_common_title(get_title(row_a), get_title(row_b))
        if not common_title:
            common_title = f"{acc_a}_{acc_b}"

        header = (
            f"{common_title} | status=PAIR | "
            f"A={acc_a} | B={acc_b} | "
            f"A_len={len(seq_a)} | B_start={b_start} | gap={gap_len} | "
            f"B_len={len(seq_b)} | right_pad={right_pad_len} | joined_len={len(joined_seq)}"
        )

        joined_entries.append((header, joined_seq))
        complete_pairs_written += 1

    # Unpaired DNA-A:
    # DNA-A on the left + A-padding until common_total_len.
    for row_a in unpaired_a_rows:
        acc_a = get_accession(row_a)

        if acc_a not in fasta_records:
            skipped_unpaired_a += 1
            continue

        seq_a = fasta_records[acc_a]["seq"]
        right_pad_len = common_total_len - len(seq_a)

        if right_pad_len < 0:
            raise ValueError(f"Invalid right padding for unpaired DNA-A {acc_a}.")

        joined_seq = seq_a + (PAD_BASE * right_pad_len)

        common_title = normalize_common_title(get_title(row_a), "")
        if not common_title:
            common_title = acc_a

        header = (
            f"{common_title} | status=A_ONLY | "
            f"A={acc_a} | B=missing | "
            f"A_len={len(seq_a)} | B_start={b_start} | "
            f"right_pad={right_pad_len} | joined_len={len(joined_seq)}"
        )

        joined_entries.append((header, joined_seq))
        unpaired_a_written += 1

    # Unpaired DNA-B:
    # A-padding until B_start + DNA-B on the right + final A-padding.
    for row_b in unpaired_b_rows:
        acc_b = get_accession(row_b)

        if acc_b not in fasta_records:
            skipped_unpaired_b += 1
            continue

        seq_b = fasta_records[acc_b]["seq"]

        left_pad_len = b_start
        right_pad_len = common_total_len - (left_pad_len + len(seq_b))

        if right_pad_len < 0:
            raise ValueError(f"Invalid right padding for unpaired DNA-B {acc_b}.")

        joined_seq = (PAD_BASE * left_pad_len) + seq_b + (PAD_BASE * right_pad_len)

        common_title = normalize_common_title("", get_title(row_b))
        if not common_title:
            common_title = acc_b

        header = (
            f"{common_title} | status=B_ONLY | "
            f"A=missing | B={acc_b} | "
            f"B_start={b_start} | left_pad={left_pad_len} | "
            f"B_len={len(seq_b)} | right_pad={right_pad_len} | joined_len={len(joined_seq)}"
        )

        joined_entries.append((header, joined_seq))
        unpaired_b_written += 1

    with open(output_fasta, "w", encoding="utf-8") as out:
        for header, seq in joined_entries:
            out.write(f">{header}\n")
            out.write(wrap_seq(seq) + "\n")

    print("Done.")
    print(f"Complete pairs in pairs.csv: {len(pairs)}")
    print(f"Unpaired DNA-A rows in unpaired_A.csv: {len(unpaired_a_rows)}")
    print(f"Unpaired DNA-B rows in unpaired_B.csv: {len(unpaired_b_rows)}")
    print(f"Complete pairs written: {complete_pairs_written}")
    print(f"Unpaired DNA-A written: {unpaired_a_written}")
    print(f"Unpaired DNA-B written: {unpaired_b_written}")
    print(f"Skipped complete pairs because of missing FASTA accession: {skipped_pairs}")
    print(f"Skipped unpaired DNA-A because of missing FASTA accession: {skipped_unpaired_a}")
    print(f"Skipped unpaired DNA-B because of missing FASTA accession: {skipped_unpaired_b}")
    print(f"Common DNA-B start position: {b_start}")
    print(f"Common output length: {common_total_len}")
    print(f"Wrote: {output_fasta}")

    if missing_accessions:
        missing_unique = sorted(set(x for x in missing_accessions if x))
        print(f"Missing accessions in FASTA: {len(missing_unique)}")
        for acc in missing_unique[:20]:
            print(f"  - {acc}")
        if len(missing_unique) > 20:
            print(f"  ... and {len(missing_unique) - 20} more")

    if warnings:
        print("Warnings:")
        for warning in warnings[:30]:
            print(f"  - {warning}")
        if len(warnings) > 30:
            print(f"  ... and {len(warnings) - 30} more")


if __name__ == "__main__":
    main()
