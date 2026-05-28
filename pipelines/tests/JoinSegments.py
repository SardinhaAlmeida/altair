# #!/usr/bin/env python3
# import argparse
# import math
# import os
# import re
# from collections import defaultdict
# from typing import Dict, Tuple, Optional, List

# # Matches any of:
# # - DNA-A, DNA A, DNA_A
# # - segment DNA-A, segment DNA A, segment DNA_A
# # - segment A
# # Same for B.

# SEG_A_PAT = re.compile(
#     r"(?i)\b(?:"
#     r"dna\s*[-_ ]?\s*a"                                  # DNA-A / DNA A / DNA_A
#     r"|(?:segment|seg|component|comp)\s+(?:dna\s*)?[-_ ]?\s*a"  # segment DNA-A / segment DNA A / segment A
#     r")\b"
# )

# SEG_B_PAT = re.compile(
#     r"(?i)\b(?:"
#     r"dna\s*[-_ ]?\s*b"
#     r"|(?:segment|seg|component|comp)\s+(?:dna\s*)?[-_ ]?\s*b"
#     r")\b"
# )

# # Remove segment tokens to build a stable "pair key" across A/B headers
# STRIP_TOKENS_PAT = re.compile(
#     r"(?i)\b(?:"
#     r"dna\s*[-_ ]?\s*[ab]"
#     r"|(?:segment|seg|component|comp)\s+(?:dna\s*)?[-_ ]?\s*[ab]"
#     r")\b"
# )

# # Matches leading accession token like PQ451365.1, PV832049.1, NC_003861.1, etc.
# LEADING_ACCESSION_PAT = re.compile(r"^\s*\S+\.\d+\s+")

# FASTA_EXTS = (".fasta", ".fa", ".fna")  # case-insensitive
# PAD_BASE = "A"  # fixed as requested
# WRAP = 80       # fixed as requested


# # -------------------------
# # FASTA IO
# # -------------------------
# def read_fasta(path: str):
#     header = None
#     chunks: List[str] = []
#     with open(path, "r", encoding="utf-8") as f:
#         for line in f:
#             line = line.strip()
#             if not line:
#                 continue
#             if line.startswith(">"):
#                 if header is not None:
#                     yield header, "".join(chunks)
#                 header = line[1:].strip()
#                 chunks = []
#             else:
#                 chunks.append(line)
#         if header is not None:
#             yield header, "".join(chunks)

# def write_fasta(path: str, records: List[Tuple[str, str]]):
#     with open(path, "w", encoding="utf-8") as f:
#         for h, s in records:
#             f.write(f">{h}\n")
#             for i in range(0, len(s), WRAP):
#                 f.write(s[i:i+WRAP] + "\n")


# # -------------------------
# # Pairing + joining logic
# # -------------------------
# def sanitize_seq(seq: str) -> str:
#     return re.sub(r"\s+", "", seq).upper()

# def detect_segment(header: str) -> Optional[str]:
#     a = bool(SEG_A_PAT.search(header))
#     b = bool(SEG_B_PAT.search(header))
#     if a and not b:
#         return "A"
#     if b and not a:
#         return "B"
#     return None

# def normalize_key(header: str) -> str:
#     h = header.strip()

#     # 1) Remove leading accession (first token like PQ451365.1)
#     h = LEADING_ACCESSION_PAT.sub("", h)

#     # 2) Remove segment/DNA A/B markers (covers: segment DNA-A, segment DNA A, DNA-B, segment B, etc.)
#     h = STRIP_TOKENS_PAT.sub(" ", h)

#     # 3) Normalize punctuation / whitespace
#     h = re.sub(r"[|,/;()\[\]{}]", " ", h)   # keep commas from causing odd spacing
#     h = re.sub(r"\s+", " ", h).strip().lower()

#     return h


# def compute_gap_len(a_len: int, step: int) -> Tuple[int, int]:
#     """
#     Auto pad_to = next multiple of step that is > a_len (gap never 0).
#     Examples with step=1000:
#       a_len=2800 -> pad_to=3000 gap=200
#       a_len=3000 -> pad_to=4000 gap=1000 (never 0)
#     """
#     if step <= 0:
#         raise ValueError("step must be > 0")
#     next_mult = int(math.ceil(a_len / step) * step)
#     if next_mult == a_len:
#         next_mult += step  # force non-zero gap
#     gap = next_mult - a_len
#     return next_mult, gap

# def join_one_fasta(in_path: str, out_path: str, step: int) -> Dict[str, object]:
#     """
#     Discards incomplete pairs.
#     Returns a report dict including excluded IDs.
#     """
#     pairs = defaultdict(lambda: {"A": None, "B": None, "Ah": None, "Bh": None})
#     total_records = 0
#     unknown_segment_records = 0

#     for header, seq in read_fasta(in_path):
#         total_records += 1
#         seg = detect_segment(header)
#         seq = sanitize_seq(seq)

#         if seg is None:
#             unknown_segment_records += 1
#             continue

#         key = normalize_key(header)
#         if pairs[key][seg] is None:  # keep first occurrence if duplicates exist
#             pairs[key][seg] = seq
#             pairs[key]["Ah" if seg == "A" else "Bh"] = header

#     out_records: List[Tuple[str, str]] = []
#     excluded_ids: List[str] = []

#     for key, p in pairs.items():
#         a = p["A"]
#         b = p["B"]
#         if not (a and b):
#             excluded_ids.append(key)
#             continue

#         a_len = len(a)
#         b_len = len(b)
#         pad_to, gap_len = compute_gap_len(a_len, step)
#         joined = a + (PAD_BASE * gap_len) + b

#         out_header = (
#             f"{key} | joined=A+gap({PAD_BASE})+B | "
#             f"A={a_len} pad_to={pad_to} gap={gap_len} B={b_len} joined={len(joined)}"
#         )
#         out_records.append((out_header, joined))

#     write_fasta(out_path, out_records)

#     return {
#         "records_in": total_records,
#         "unknown_segment_records": unknown_segment_records,
#         "pair_keys_total": len(pairs),
#         "paired_out": len(out_records),
#         "incomplete_excluded": len(excluded_ids),
#         "excluded_ids": excluded_ids,
#     }


# # -------------------------
# # Batch traversal
# # -------------------------
# def find_variant_dirs(root: str) -> List[str]:
#     return [
#         os.path.join(root, d)
#         for d in sorted(os.listdir(root))
#         if os.path.isdir(os.path.join(root, d))
#     ]

# def find_fasta_files_in_dir(dirpath: str) -> List[str]:
#     out = []
#     for name in os.listdir(dirpath):
#         p = os.path.join(dirpath, name)
#         if not os.path.isfile(p):
#             continue
#         if name.lower().endswith(FASTA_EXTS):
#             out.append(p)
#     return sorted(out)

# def main():
#     ap = argparse.ArgumentParser(
#         description="For each variant folder under pathogens/, join Begomovirus DNA-A and DNA-B with an auto poly-A gap and write *_joined.fasta."
#     )
#     ap.add_argument("pathogens_root", help="Path to pathogens root (contains variant subfolders)")
#     ap.add_argument("--step", type=int, default=1000, help="Boundary step for gap placement (default: 1000)")
#     ap.add_argument("--print-excluded", action="store_true",
#                     help="Print all excluded IDs (can be long). If off, prints only first 30.")
#     args = ap.parse_args()

#     root = os.path.abspath(args.pathogens_root)
#     if not os.path.isdir(root):
#         raise SystemExit(f"[ERROR] Not a directory: {root}")
    
#     # --- Cleanup: remove previously generated outputs ---
#     removed = 0
#     for dirpath, _, files in os.walk(root):
#         for name in files:
#             if name.lower().endswith("_joined.fasta"):
#                 p = os.path.join(dirpath, name)
#                 try:
#                     os.remove(p)
#                     removed += 1
#                 except OSError as e:
#                     print(f"[warn] Could not delete {p}: {e}")

#     print(f"[i] Cleanup: removed {removed} *_joined.fasta file(s)")


#     variant_dirs = find_variant_dirs(root)
#     if not variant_dirs:
#         raise SystemExit(f"[ERROR] No variant subfolders found in: {root}")

#     processed = 0

#     for vdir in variant_dirs:
#         fasta_files = find_fasta_files_in_dir(vdir)
#         if not fasta_files:
#             continue

#         variant = os.path.basename(vdir)
#         print(f"\n[VARIANT] {variant}")

#         for in_path in fasta_files:
#             base, _ = os.path.splitext(in_path)
#             out_path = base + "_joined.fasta"

#             report = join_one_fasta(in_path, out_path, step=args.step)

#             excluded = report["excluded_ids"]
#             excluded_n = report["incomplete_excluded"]

#             print(f"  File: {os.path.basename(in_path)}")
#             print(f"    Input records (FASTA entries): {report['records_in']}")
#             print(f"    Pair keys detected:           {report['pair_keys_total']}")
#             print(f"    Output joined genomes:        {report['paired_out']}")
#             print(f"    Incomplete excluded:          {excluded_n}")
#             print(f"    Unknown-segment records:      {report['unknown_segment_records']}")
#             print(f"    Wrote: {os.path.basename(out_path)}")

#             if excluded_n > 0:
#                 if args.print_excluded:
#                     shown = excluded
#                 else:
#                     shown = excluded[:30]
#                 print(f"    Excluded IDs (showing {len(shown)}{' of ' + str(excluded_n) if len(shown) < excluded_n else ''}):")
#                 for k in shown:
#                     print(f"      - {k}")

#             processed += 1

#     print(f"\n[i] Done. FASTA files processed: {processed}")

# if __name__ == "__main__":
#     main()
#!/usr/bin/env python3
import csv
import math
import re
import sys
from pathlib import Path

WRAP = 80
PAD_BASE = "A"
STEP = 500  # DNA-B starts at a common multiple of 500, at least 500 after longest DNA-A


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


def extract_accession(header):
    # Takes first token from FASTA header
    return header.split()[0]


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


def normalize_common_title(title_a, title_b):
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
    # next multiple of 500, but always at least 500 after longest A
    return math.ceil((max_a_len + step) / step) * step


def wrap_seq(seq, width=WRAP):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def load_pairs_csv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if len(rows) % 2 != 0:
        raise ValueError("pairs.csv has an odd number of data rows. It must be A/B ordered in pairs.")

    pairs = []

    for i in range(0, len(rows), 2):
        row1 = rows[i]
        row2 = rows[i + 1]

        comp1 = classify_component(row1.get("Segment", ""), row1.get("GenBank_Title", ""))
        comp2 = classify_component(row2.get("Segment", ""), row2.get("GenBank_Title", ""))

        if comp1 != "A" or comp2 != "B":
            raise ValueError(
                f"Rows {i+2} and {i+3} are not ordered as A then B.\n"
                f"  Row {i+2}: Accession={row1.get('Accession','')} Component={comp1}\n"
                f"  Row {i+3}: Accession={row2.get('Accession','')} Component={comp2}"
            )

        pairs.append((row1, row2))

    return pairs


def main():
    if len(sys.argv) != 4:
        print("Usage: python3 JoinSegments.py pairs.csv input_sequences.fasta output_joined.fasta")
        sys.exit(1)

    pairs_csv = Path(sys.argv[1])
    input_fasta = Path(sys.argv[2])
    output_fasta = Path(sys.argv[3])

    if not pairs_csv.exists():
        raise FileNotFoundError(f"pairs.csv not found: {pairs_csv}")
    if not input_fasta.exists():
        raise FileNotFoundError(f"input FASTA not found: {input_fasta}")

    pairs = load_pairs_csv(pairs_csv)
    fasta_records = read_fasta(input_fasta)

    joined_entries = []
    missing_accessions = []
    a_lengths = []

    # first pass: collect A lengths
    for row_a, row_b in pairs:
        acc_a = (row_a.get("Accession") or "").strip()
        acc_b = (row_b.get("Accession") or "").strip()

        if acc_a not in fasta_records or acc_b not in fasta_records:
            if acc_a not in fasta_records:
                missing_accessions.append(acc_a)
            if acc_b not in fasta_records:
                missing_accessions.append(acc_b)
            continue

        seq_a = fasta_records[acc_a]["seq"]
        a_lengths.append(len(seq_a))

    if not a_lengths:
        raise ValueError("No valid A/B pairs found in FASTA for joining.")

    b_start = compute_b_start(max(a_lengths), STEP)

    # second pass: build joined sequences
    for row_a, row_b in pairs:
        acc_a = (row_a.get("Accession") or "").strip()
        acc_b = (row_b.get("Accession") or "").strip()

        if acc_a not in fasta_records or acc_b not in fasta_records:
            continue

        seq_a = fasta_records[acc_a]["seq"]
        seq_b = fasta_records[acc_b]["seq"]

        gap_len = b_start - len(seq_a)
        if gap_len < 0:
            raise ValueError(f"Negative gap for pair {acc_a} / {acc_b}. This should never happen.")

        joined_seq = seq_a + (PAD_BASE * gap_len) + seq_b

        common_title = normalize_common_title(
            row_a.get("GenBank_Title", ""),
            row_b.get("GenBank_Title", "")
        )

        header = (
            f"{common_title} | "
            f"A={acc_a} | B={acc_b} | "
            f"A_len={len(seq_a)} | B_start={b_start} | gap={gap_len} | B_len={len(seq_b)} | joined_len={len(joined_seq)}"
        )

        joined_entries.append((header, joined_seq))

    with open(output_fasta, "w", encoding="utf-8") as out:
        for header, seq in joined_entries:
            out.write(f">{header}\n")
            out.write(wrap_seq(seq) + "\n")

    print(f"Done.")
    print(f"Pairs in CSV: {len(pairs)}")
    print(f"Joined sequences written: {len(joined_entries)}")
    print(f"Common DNA-B start position: {b_start}")
    if missing_accessions:
        missing_unique = sorted(set(x for x in missing_accessions if x))
        print(f"Missing accessions in FASTA: {len(missing_unique)}")
        for acc in missing_unique[:20]:
            print(f"  - {acc}")
        if len(missing_unique) > 20:
            print(f"  ... and {len(missing_unique) - 20} more")


if __name__ == "__main__":
    main()