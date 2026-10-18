#!/usr/bin/env python3
"""Classify and pair DNA-A/DNA-B components from an NCBI Virus metadata CSV.

This script is for bipartite viruses (e.g. bipartite begomoviruses). It is run
twice by `pipeline2.sh prepare`:

1. Before downloading (no --fasta): every complete record is labelled DNA-A or
   DNA-B from its Segment field, with GenBank_Title as a fallback, and the
   records to download are written to the output directory.
2. After downloading (--fasta): records that are still unlabelled are labelled
   by k-mer similarity to the labelled records of the same dataset, RefSeq
   records identical to another record are removed, and all records are paired
   again.

Optional quality filters (off by default) exclude UNVERIFIED records, records
from unselected hosts, records with too many ambiguous bases and records of
unexpected length. Filters run before similarity labelling, so excluded records
are neither labelled nor used as references; each is listed with its reason in
excluded_records.csv.

Optional files in the folder of the input CSV, read in both passes:
  segments_options.txt   command-line options, one per line (e.g.
                         "--host=Manihot esculenta"), so that both passes use
                         the same settings and the selection is recorded.
  manual_components.csv  columns Accession,Component (A or B); overrides labels.
  manual_pairs.csv       columns Accession_A,Accession_B; pairs these records
                         before automatic pairing (e.g. after reviewing
                         ambiguous.csv).
"""
import argparse
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = ["Accession", "Segment", "GenBank_Title"]

EXPECTED_COLUMNS = [
    "Accession",
    "BioSample",
    "Organism_Name",
    "Species",
    "Isolate",
    "Segment",
    "GenBank_Title",
    "Length",
    "Nuc_Completeness",
    "Geo_Location",
    "Country",
    "Host",
    "Collection_Date",
]

PAIR_FIELDS = [
    "BioSample",
    "Organism_Name",
    "Species",
    "Nuc_Completeness",
    "Geo_Location",
    "Country",
    "Host",
    "Collection_Date",
]

# Similarity labelling. Thresholds were set from seven cassava mosaic virus
# datasets: the nearest labelled record of the same component scored >= 0.30,
# the nearest record of the other component <= 0.20 (shared common region),
# and shuffled sequences <= 0.13.
K = 8
DEFAULT_MIN_SIMILARITY = 0.30
DEFAULT_MIN_MARGIN = 0.10

# Segment field: the whole value must name the component ("A", "DNA-B",
# "segment A", "DNA1", ...). Bare digits are not accepted, since "1"/"2" are
# also used by segmented RNA viruses.
SEGMENT_LETTER = re.compile(r"(?:(?:segment|component)[\s_-]*)?(?:dna[\s_-]*)?([ab])", re.I)
SEGMENT_NUMBER = re.compile(r"(?:(?:segment|component)[\s_-]*)?dna[\s_-]*([12])", re.I)
# GenBank_Title: explicit component wording anywhere in the title
# ("segment A", "segment: A", "component B", "DNA-A", "DNA A", "segment DNA1").
TITLE_COMPONENT = re.compile(
    r"\b(?:(?:segment|component)[\s_:-]+(?:dna[\s_-]?)?([ab])|dna[\s_-]?([ab12]))\b", re.I)
TITLE_COMPLETENESS = re.compile(r"\bcomplete\s+(?:sequence|genome)\b", re.I)
SATELLITE = re.compile(r"satellite", re.I)
NUMBER_TO_COMPONENT = {"1": "A", "2": "B", "A": "A", "B": "B"}

COMPLEMENT = str.maketrans("ACGT", "TGCA")
BASE4 = str.maketrans("ACGT", "0123")
UNAMBIGUOUS_RUN = re.compile(r"[ACGT]+")


def fail(message):
    raise SystemExit(f"ERROR: {message}")


def clean_text(x: str) -> str:
    x = str(x).strip().lower()
    x = re.sub(r"\s+", " ", x)
    return x


def normalize_accession(accession: str) -> str:
    return re.sub(r"\.\d+$", "", str(accession or "").strip())


def segment_component(segment: str):
    """Return "A", "B", or None when the Segment value does not name a component."""
    value = segment.strip()
    match = SEGMENT_LETTER.fullmatch(value) or SEGMENT_NUMBER.fullmatch(value)
    return NUMBER_TO_COMPONENT[match.group(1).upper()] if match else None


def title_components(title: str) -> set:
    return {NUMBER_TO_COMPONENT[(m.group(1) or m.group(2)).upper()]
            for m in TITLE_COMPONENT.finditer(title)}


def classify_component(segment: str, title: str):
    """Return (component, source, note) from the metadata of one record.

    Segment is the primary source and GenBank_Title the fallback. Records whose
    two sources disagree are left UNKNOWN, so that the similarity pass decides.
    """
    segment, title = str(segment).strip(), str(title).strip()

    if SATELLITE.search(f"{segment} {title}"):
        return "UNKNOWN", "", "satellite record"

    from_title = title_components(title)

    if segment:
        from_segment = segment_component(segment)
        if from_segment is None:
            return "UNKNOWN", "", f"unrecognised Segment value '{segment}'"
        if from_title and from_title != {from_segment}:
            return ("UNKNOWN", "",
                    f"Segment says {from_segment} but GenBank_Title says "
                    f"{'/'.join(sorted(from_title))}")
        return from_segment, "Segment", ""

    if len(from_title) == 1:
        return from_title.pop(), "GenBank_Title", ""
    if from_title:
        return "UNKNOWN", "", "GenBank_Title names both components"
    return "UNKNOWN", "", "no component in Segment or GenBank_Title"


def normalize_isolate_for_pairing(isolate: str) -> str:
    x = clean_text(isolate)
    if not x:
        return ""

    x = x.replace("_", "-")
    x = re.sub(r"\s+", " ", x)

    # ACMV-A-CI008-01_A_M -> ACMV-CI008-01-M
    # ACMV-B-CI008-01_B_M -> ACMV-CI008-01-M
    x = re.sub(r"^([a-z0-9]+)-[ab]-", r"\1-", x)
    x = re.sub(r"(?<=-)[ab](?=-)", "", x)

    # remove trailing A/B after sample number: CMD03A -> CMD03
    x = re.sub(r"(?<=\d)[ab]$", "", x)

    # remove terminal -a / -b
    x = re.sub(r"([a-z0-9])-[ab]$", r"\1", x)

    # normalize leftovers
    x = re.sub(r"-+", "-", x).strip("- ")
    return x


def title_pair_core(title: str) -> str:
    """Title without component and completeness wording, used when Isolate is empty.

    "SLCMV CaKT2 DNA, segment: A, complete sequence" and its DNA-B title give
    the same core, while CaKT2 and CaKT6 stay distinct.
    """
    x = TITLE_COMPONENT.sub(" ", str(title))
    x = TITLE_COMPLETENESS.sub(" ", x)
    x = re.sub(r"[\s,;]+", " ", x).strip(" ,;")
    return normalize_isolate_for_pairing(x)


def isolate_pair_core(isolate: str, title: str) -> str:
    return normalize_isolate_for_pairing(isolate) or title_pair_core(title)


def build_pair_key(row) -> str:
    parts = [clean_text(row[col]) for col in PAIR_FIELDS if col in row.index]

    # use isolate core if available
    parts.append(row["Isolate_Pair_Core"])

    return " || ".join(parts)


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def load_metadata(input_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(input_csv, dtype=str).fillna("")

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        fail(f"{input_csv} lacks required columns {missing}. "
             "Export the metadata from NCBI Virus with these columns.")
    optional_missing = [c for c in EXPECTED_COLUMNS if c not in df.columns and c not in REQUIRED_COLUMNS]
    if optional_missing:
        print("Warning: missing optional columns (not used for pairing):", optional_missing)

    df["Accession_Key"] = df["Accession"].map(normalize_accession)
    if (df["Accession_Key"] == "").any():
        fail(f"{input_csv} has rows with an empty Accession.")
    duplicated = sorted(df.loc[df["Accession_Key"].duplicated(), "Accession"])
    if duplicated:
        fail(f"{input_csv} lists these accessions more than once: {duplicated}")

    if "Nuc_Completeness" in df.columns:
        df = df[df["Nuc_Completeness"].str.lower() == "complete"].copy()
    else:
        print("Warning: no Nuc_Completeness column; all records are treated as complete.")
    return df.reset_index(drop=True)


def read_fasta(path: Path) -> dict:
    records, accession, chunks = {}, None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if accession is not None:
                    records[accession] = "".join(chunks).upper()
                accession = normalize_accession(line[1:].split()[0]) if line[1:].split() else ""
                chunks = []
            elif line:
                chunks.append(line.replace(" ", ""))
    if accession is not None:
        records[accession] = "".join(chunks).upper()
    return records


def read_manual_csv(path: Path, columns: list) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str).fillna("")
    missing = [c for c in columns if c not in table.columns]
    if missing:
        fail(f"{path} must have the columns {columns}; missing {missing}.")
    return table


def resolve_manual_path(option, default: Path):
    """An explicit path must exist; the default file is optional."""
    if option is not None:
        if not option.is_file():
            fail(f"manual file not found: {option}")
        return option
    return default if default.is_file() else None


# ---------------------------------------------------------------------------
# Component labels
# ---------------------------------------------------------------------------

def label_from_metadata(df: pd.DataFrame) -> None:
    labels = [classify_component(s, t) for s, t in zip(df["Segment"], df["GenBank_Title"])]
    df["Component"] = [c for c, _, _ in labels]
    df["Component_Source"] = [s for _, s, _ in labels]
    df["Component_Note"] = [n for _, _, n in labels]


def index_of(df: pd.DataFrame, accession: str, source: Path) -> int:
    matches = df.index[df["Accession_Key"] == normalize_accession(accession)]
    if len(matches) == 0:
        fail(f"{source}: accession '{accession}' is not among the complete records "
             "of the input CSV.")
    return matches[0]


def apply_manual_components(df: pd.DataFrame, path: Path) -> int:
    table = read_manual_csv(path, ["Accession", "Component"])
    seen = set()
    for _, entry in table.iterrows():
        idx = index_of(df, entry["Accession"], path)
        if idx in seen:
            fail(f"{path}: accession '{entry['Accession']}' is listed more than once.")
        seen.add(idx)
        component = segment_component(entry["Component"])
        if component is None:
            fail(f"{path}: Component for '{entry['Accession']}' must be A or B, "
                 f"not '{entry['Component']}'.")
        df.loc[idx, ["Component", "Component_Source", "Component_Note"]] = \
            [component, "manual", "set in manual_components.csv"]
    return len(seen)


def read_manual_pairs(df: pd.DataFrame, path: Path) -> list:
    """Return [(index_a, index_b)] and set both components accordingly."""
    table = read_manual_csv(path, ["Accession_A", "Accession_B"])
    pairs, seen = [], set()
    for _, entry in table.iterrows():
        idx_a = index_of(df, entry["Accession_A"], path)
        idx_b = index_of(df, entry["Accession_B"], path)
        for idx in (idx_a, idx_b):
            if idx in seen or idx_a == idx_b:
                fail(f"{path}: accession '{df.at[idx, 'Accession']}' appears in more "
                     "than one manual pair.")
            seen.add(idx)
        for idx, component in ((idx_a, "A"), (idx_b, "B")):
            if df.at[idx, "Component"] not in (component, "UNKNOWN"):
                print(f"Warning: {path.name} pairs {df.at[idx, 'Accession']} as DNA-{component}, "
                      f"but it was labelled DNA-{df.at[idx, 'Component']} "
                      f"({df.at[idx, 'Component_Source']}); the manual pair is used.")
            df.loc[idx, ["Component", "Component_Source", "Component_Note"]] = \
                [component, "manual", "set by manual_pairs.csv"]
        pairs.append((idx_a, idx_b))
    return pairs


def find_refseq_copies(df: pd.DataFrame, sequences: dict) -> dict:
    """Map each RefSeq (NC_) accession to a non-RefSeq record with the identical sequence."""
    by_sequence = defaultdict(list)
    for accession in df["Accession_Key"]:
        if sequences.get(accession):
            by_sequence[sequences[accession]].append(accession)
    copies = {}
    for accessions in by_sequence.values():
        originals = sorted(a for a in accessions if not a.startswith("NC_"))
        for accession in accessions:
            if originals and accession.startswith("NC_"):
                copies[accession] = originals[0]
    return copies


def popcount(x: int) -> int:
    return x.bit_count() if hasattr(x, "bit_count") else bin(x).count("1")


def kmer_set(sequence: str, k: int = K):
    """Canonical k-mers of a circular sequence as (bit set, size).

    Windows with ambiguous bases are skipped. Bit i is set when the k-mer with
    base-4 code i occurs, so set intersection is a bitwise AND.
    """
    n = len(sequence)
    if n < k:
        return 0, 0
    circular = sequence + sequence[:k - 1]
    reverse = circular.translate(COMPLEMENT)[::-1]
    total = len(circular)
    bits = bytearray(4 ** k // 8)
    for run in UNAMBIGUOUS_RUN.finditer(circular):
        for i in range(run.start(), min(run.end() - k, n - 1) + 1):
            forward = circular[i:i + k]
            backward = reverse[total - i - k:total - i]
            code = int(min(forward, backward).translate(BASE4), 4)
            bits[code >> 3] |= 1 << (code & 7)
    as_int = int.from_bytes(bits, "little")
    return as_int, popcount(as_int)


def nearest(query, candidates: list, kmers: dict, exclude: str):
    """Return (accession, cosine similarity) of the most similar candidate, or ("", nan)."""
    query_bits, query_size = query
    best_accession, best_similarity = "", float("nan")
    for accession in candidates:
        bits, size = kmers[accession]
        if accession == exclude or not size or not query_size:
            continue
        similarity = popcount(query_bits & bits) / math.sqrt(query_size * size)
        if best_accession == "" or similarity > best_similarity:
            best_accession, best_similarity = accession, similarity
    return best_accession, best_similarity


def label_by_similarity(df: pd.DataFrame, sequences: dict, keep: pd.Series,
                        min_similarity: float, min_margin: float) -> list:
    """Label UNKNOWN records from their nearest labelled records and check the labelled ones.

    References are the records labelled before this step (metadata or manual), so
    the result does not depend on record order. Returns one report row per record.
    """
    rows = df[keep]
    kmers = {a: kmer_set(sequences[a]) for a in rows["Accession_Key"] if sequences.get(a)}
    references = {c: sorted(a for a, comp in zip(rows["Accession_Key"], rows["Component"])
                            if comp == c and a in kmers) for c in "AB"}
    report = []

    for idx, row in rows.iterrows():
        accession, label = row["Accession_Key"], row["Component"]
        entry = {"Accession": row["Accession"], "Length": len(sequences.get(accession, "")),
                 "Label_Before": label, "Label_Source": row["Component_Source"],
                 "Nearest_A": "", "Similarity_A": "", "Nearest_B": "", "Similarity_B": "",
                 "Label_After": label}
        if accession not in kmers:
            entry["Note"] = "sequence not in FASTA"
            report.append(entry)
            continue

        best = {c: nearest(kmers[accession], references[c], kmers, accession) for c in "AB"}
        for c in "AB":
            entry[f"Nearest_{c}"] = best[c][0]
            entry[f"Similarity_{c}"] = "" if not best[c][0] else f"{best[c][1]:.3f}"
        lacking = [c for c in "AB" if not best[c][0]]

        if label == "UNKNOWN":
            if lacking:
                entry["Note"] = (f"not labelled: no labelled DNA-{lacking[0]} in the dataset "
                                 "to compare with")
            else:
                winner = "A" if best["A"][1] >= best["B"][1] else "B"
                loser = "B" if winner == "A" else "A"
                score, margin = best[winner][1], best[winner][1] - best[loser][1]
                if score >= min_similarity and margin >= min_margin:
                    entry["Label_After"] = winner
                    entry["Note"] = f"labelled DNA-{winner} by similarity"
                    df.loc[idx, ["Component", "Component_Source", "Component_Note"]] = \
                        [winner, "similarity",
                         f"nearest {best[winner][0]} ({score:.3f}); was: {row['Component_Note']}"]
                else:
                    entry["Note"] = (f"not labelled: nearest DNA-{winner} scores {score:.3f} "
                                     f"(min {min_similarity}), margin {margin:.3f} (min {min_margin})")
        elif lacking:
            entry["Note"] = "label kept; not checked (no other labelled record to compare with)"
        else:
            other = "B" if label == "A" else "A"
            own, rival = best[label][1], best[other][1]
            if own < rival:
                entry["Note"] = f"label kept, but nearer to DNA-{other}: check this record"
            elif own < min_similarity or own - rival < min_margin:
                entry["Note"] = "label kept, but similarity is weak: check this record"
            else:
                entry["Note"] = "label agrees with similarity"
        report.append(entry)
    return report


# ---------------------------------------------------------------------------
# Quality filters (all off by default)
# ---------------------------------------------------------------------------
# Each reason starts with a category followed by ":" (or is only a category),
# which the summary uses to count exclusions.

def exclude(df: pd.DataFrame, keep: pd.Series, mask: pd.Series, reasons) -> None:
    """Exclude kept records selected by mask; a record keeps its first reason."""
    hit = keep & mask
    df.loc[hit, "Excluded_Reason"] = reasons[hit] if isinstance(reasons, pd.Series) else reasons
    keep[hit] = False


def filter_by_metadata(df: pd.DataFrame, keep: pd.Series, args) -> None:
    """Filters that need only the metadata, applied before download."""
    if args.exclude_unverified:
        unverified = df["GenBank_Title"].str.contains(r"\bunverified\b", case=False, regex=True)
        exclude(df, keep, unverified, "UNVERIFIED GenBank record")

    if args.host:
        if "Host" not in df.columns:
            fail("--host needs a Host column in the input CSV.")
        wanted = {clean_text(h) for h in args.host}
        host = df["Host"].map(clean_text)
        exclude(df, keep, host == "", "host not recorded")
        exclude(df, keep, ~host.isin(wanted), "host not selected: " + df["Host"])


def filter_by_sequence(df: pd.DataFrame, keep: pd.Series, sequences: dict, args) -> None:
    """Filters that need the downloaded sequences, applied before similarity labelling.

    Records without a downloaded sequence are left to the later steps, which
    report them as missing.
    """
    seqs = df["Accession_Key"].map(lambda a: sequences.get(a, ""))
    lengths = seqs.str.len()
    has_sequence = lengths > 0

    if args.max_ambiguous is not None:
        fraction = seqs.map(lambda s: len(UNAMBIGUOUS_RUN.sub("", s)) / len(s) if s else 0.0)
        reasons = fraction.map(lambda f: f"ambiguous bases: {f:.1%} (max {args.max_ambiguous:.1%})")
        exclude(df, keep, has_sequence & (fraction > args.max_ambiguous), reasons)

    if args.max_length_deviation is not None:
        labelled = keep & has_sequence & df["Component"].isin(["A", "B"])
        if not labelled.any():
            return
        # Unlabelled records are compared with the median of all labelled records.
        medians = {c: lengths[labelled & (df["Component"] == c)].median() for c in "AB"}
        median = df["Component"].map(medians).fillna(lengths[labelled].median())
        deviation = (lengths - median).abs() / median
        reasons = [f"length: {n} bp, {d:.0%} from median {m:.0f} (max {args.max_length_deviation:.0%})"
                   for n, d, m in zip(lengths, deviation, median)]
        exclude(df, keep, has_sequence & (deviation > args.max_length_deviation),
                pd.Series(reasons, index=df.index))


# ---------------------------------------------------------------------------
# Pairing and output
# ---------------------------------------------------------------------------

def output_row(row, columns):
    return {col: row[col] for col in columns}


def run(args):
    input_csv = Path(args.input_csv)
    out_dir = Path(args.out_dir) if args.out_dir else input_csv.parent / "paired_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_metadata(input_csv)
    label_from_metadata(df)
    from_metadata = Counter(df["Component"])

    manual_components = resolve_manual_path(args.manual_components,
                                            input_csv.parent / "manual_components.csv")
    manual_pairs_path = resolve_manual_path(args.manual_pairs,
                                            input_csv.parent / "manual_pairs.csv")
    n_manual_components = apply_manual_components(df, manual_components) if manual_components else 0
    manual_pairs = read_manual_pairs(df, manual_pairs_path) if manual_pairs_path else []

    keep = pd.Series(True, index=df.index)
    df["Excluded_Reason"] = ""
    filter_by_metadata(df, keep, args)
    report = []
    if args.fasta:
        sequences = read_fasta(Path(args.fasta))
        # Only kept records count as originals, so a RefSeq copy of an excluded
        # record is judged on its own.
        copies = find_refseq_copies(df[keep], sequences)
        exclude(df, keep, df["Accession_Key"].isin(copies),
                "RefSeq copy: identical to " + df["Accession_Key"].map(copies).fillna(""))
        filter_by_sequence(df, keep, sequences, args)
        report = label_by_similarity(df, sequences, keep, args.min_similarity, args.min_margin)

    manual_indices = {i for pair in manual_pairs for i in pair} | set(
        df.index[df["Component_Source"] == "manual"])
    excluded_manual = [i for i in sorted(manual_indices) if not keep[i]]
    if excluded_manual:
        details = "; ".join(f"{df.at[i, 'Accession']} ({df.at[i, 'Excluded_Reason']})"
                            for i in excluded_manual)
        fail(f"records named in the manual files were excluded: {details}. Remove them from "
             "the manual files (for a RefSeq copy, name the original) or relax the filter.")

    known_components = set(df.loc[keep, "Component"]) & {"A", "B"}
    if not known_components:
        segments = Counter(df["Segment"]).most_common(5)
        fail("no record could be labelled DNA-A or DNA-B (most common Segment values: "
             f"{segments}). This script handles bipartite DNA-A/DNA-B viruses only; for "
             "non-segmented genomes use the 'complete' workflow, or label records in "
             "manual_components.csv.")

    df["Isolate_Clean"] = df["Isolate"].apply(clean_text) if "Isolate" in df.columns else ""
    isolates = df["Isolate"] if "Isolate" in df.columns else pd.Series("", index=df.index)
    df["Isolate_Pair_Core"] = [isolate_pair_core(i, t) for i, t in zip(isolates, df["GenBank_Title"])]
    df["Pair_Key"] = df.apply(build_pair_key, axis=1)

    # debug file
    df.drop(columns=["Accession_Key"]).to_csv(out_dir / "metadata_with_keys.csv", index=False)

    output_columns = [c for c in EXPECTED_COLUMNS if c in df.columns] + ["Component"]
    pairs_rows, ambiguous_rows, unpaired_rows = [], [], []

    for idx_a, idx_b in manual_pairs:
        pairs_rows.append(output_row(df.loc[idx_a], output_columns))
        pairs_rows.append(output_row(df.loc[idx_b], output_columns))

    in_manual_pair = df.index.isin([i for pair in manual_pairs for i in pair])
    known = df[keep & ~in_manual_pair & df["Component"].isin(["A", "B"])]

    ambiguous_groups = 0
    for _, group in known.groupby("Pair_Key", sort=True):
        a_rows = group[group["Component"] == "A"]
        b_rows = group[group["Component"] == "B"]

        if len(a_rows) == 1 and len(b_rows) == 1:
            pairs_rows.append(output_row(a_rows.iloc[0], output_columns))
            pairs_rows.append(output_row(b_rows.iloc[0], output_columns))
        elif len(a_rows) == 0 or len(b_rows) == 0:
            unpaired_rows.extend(output_row(row, output_columns) for _, row in group.iterrows())
        else:
            ambiguous_groups += 1
            ambiguous_rows.extend(output_row(row, output_columns + ["Pair_Key"])
                                  for _, row in group.iterrows())

    unknown = df[keep & (df["Component"] == "UNKNOWN")].rename(columns={"Component_Note": "Reason"})
    unknown_columns = [c for c in EXPECTED_COLUMNS if c in df.columns] + ["Reason"]

    pairs_df = pd.DataFrame(pairs_rows, columns=output_columns)
    ambiguous_df = pd.DataFrame(ambiguous_rows, columns=output_columns + ["Pair_Key"])
    unpaired_df = pd.DataFrame(unpaired_rows, columns=output_columns)
    unpaired_a_df = unpaired_df[unpaired_df["Component"] == "A"]
    unpaired_b_df = unpaired_df[unpaired_df["Component"] == "B"]

    # Every kept record must reach exactly one output file.
    written = (list(pairs_df["Accession"]) + list(ambiguous_df["Accession"])
               + list(unpaired_df["Accession"]) + list(unknown["Accession"]))
    if len(written) != len(set(written)) or set(written) != set(df.loc[keep, "Accession"]):
        raise RuntimeError("internal error: records were lost or duplicated while pairing")

    pairs_df.to_csv(out_dir / "pairs.csv", index=False)
    ambiguous_df.to_csv(out_dir / "ambiguous.csv", index=False)
    unpaired_df.to_csv(out_dir / "unpaired.csv", index=False)
    unpaired_a_df.to_csv(out_dir / "unpaired_A.csv", index=False)
    unpaired_b_df.to_csv(out_dir / "unpaired_B.csv", index=False)
    unknown[unknown_columns].to_csv(out_dir / "unknown_component.csv", index=False)

    excluded = df[~keep].rename(columns={"Excluded_Reason": "Reason"})
    excluded_columns = ["Accession", "Reason", "Component"] + [
        c for c in ("Host", "Length", "GenBank_Title") if c in df.columns]
    excluded[excluded_columns].to_csv(out_dir / "excluded_records.csv", index=False)

    similarity_path = out_dir / "similarity_report.csv"
    if args.fasta:
        pd.DataFrame(report, columns=[
            "Accession", "Length", "Label_Before", "Label_Source", "Nearest_A", "Similarity_A",
            "Nearest_B", "Similarity_B", "Label_After", "Note"]).to_csv(similarity_path, index=False)
    else:
        # A report from an earlier similarity pass would no longer match these outputs.
        similarity_path.unlink(missing_ok=True)

    sources = Counter(df.loc[keep & df["Component"].isin(["A", "B"]), "Component_Source"])
    print("Done.")
    if args.options_file:
        print(f"Options from {args.options_file}: {' '.join(args.file_options)}")
    print(f"Complete records: {len(df)}")
    print(f"Labelled from metadata: {from_metadata['A']} DNA-A, {from_metadata['B']} DNA-B, "
          f"{from_metadata['UNKNOWN']} unknown")
    if n_manual_components or manual_pairs:
        print(f"Manual: {n_manual_components} component labels, {len(manual_pairs)} pairs")
    if args.fasta:
        flagged = [r for r in report if r["Label_Before"] != "UNKNOWN" and "check" in r["Note"]]
        print(f"Labelled by similarity: {sources['similarity']} "
              f"(min similarity {args.min_similarity}, min margin {args.min_margin})")
        if flagged:
            print(f"Warning: {len(flagged)} labelled records disagree with similarity; "
                  f"see 'check' notes in {similarity_path}")
    else:
        print("Similarity labelling: not run (no --fasta)")
        if args.max_ambiguous is not None or args.max_length_deviation is not None:
            print("Sequence filters (ambiguous bases, length): applied after download (--fasta)")
    print(f"Excluded records: {len(excluded)} (see {out_dir / 'excluded_records.csv'})")
    for reason, count in Counter(r.split(":")[0] for r in excluded["Reason"]).most_common():
        print(f"  {count:4d}  {reason}")
    print(f"Pairs written: {len(pairs_df)} rows ({len(pairs_df) // 2} A/B pairs)")
    print(f"Unpaired rows: {len(unpaired_df)}")
    print(f"Unpaired A rows: {len(unpaired_a_df)}")
    print(f"Unpaired B rows: {len(unpaired_b_df)}")
    print(f"Unknown component rows: {len(unknown)}")
    for reason, count in Counter(unknown["Reason"]).most_common():
        print(f"  {count:4d}  {reason}")
    if ambiguous_groups:
        print(f"Warning: {len(ambiguous_df)} records in {ambiguous_groups} groups could not be "
              f"paired automatically and are not downloaded: {out_dir / 'ambiguous.csv'}")
        print("  To use them, list chosen pairs in manual_pairs.csv (Accession_A,Accession_B) "
              "and rerun.")
        for key, group in ambiguous_df.groupby("Pair_Key", sort=True):
            members = ", ".join(f"{a} ({c})" for a, c in zip(group["Accession"], group["Component"]))
            print(f"  - {members}")
    else:
        print("Ambiguous rows: 0")


OPTIONS_FILE = "segments_options.txt"


def read_options_file(path: Path) -> list:
    """One command-line option per line, e.g. "--host=Manihot esculenta"; "#" starts a comment."""
    lines = (line.split("#", 1)[0].strip() for line in path.read_text(encoding="utf-8").splitlines())
    return [line for line in lines if line]


def fraction(value: str) -> float:
    number = float(value)
    if not 0 <= number <= 1:
        raise argparse.ArgumentTypeError(f"must be between 0 and 1, not {value}")
    return number


def parse_args(argv=None):
    """Parse options; segments_options.txt next to the input CSV is read first if present.

    Options given on the command line are applied after the file, so they take
    precedence (or, for --host, add to it).
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    path = args.options_file or Path(args.input_csv).parent / OPTIONS_FILE
    if args.options_file and not path.is_file():
        fail(f"options file not found: {path}")
    file_options = read_options_file(path) if path.is_file() else []
    if file_options:
        args = parser.parse_args(file_options + argv)
    args.options_file = path if file_options else None
    args.file_options = file_options
    return args


def build_parser():
    parser = argparse.ArgumentParser(
        description="Create DNA-A/DNA-B pairs from a sequences metadata CSV.")
    parser.add_argument("input_csv", help="Input CSV file, e.g. sequences.csv")
    parser.add_argument("-o", "--out-dir", default=None,
                        help="Output directory. Default: paired_output next to the input CSV")
    parser.add_argument("--fasta", default=None,
                        help="Downloaded FASTA of the records; enables similarity labelling "
                             "and removal of RefSeq copies")
    parser.add_argument("--manual-components", type=Path, default=None,
                        help="CSV with Accession,Component. Default: manual_components.csv "
                             "next to the input CSV, if present")
    parser.add_argument("--manual-pairs", type=Path, default=None,
                        help="CSV with Accession_A,Accession_B. Default: manual_pairs.csv "
                             "next to the input CSV, if present")
    parser.add_argument("--min-similarity", type=float, default=DEFAULT_MIN_SIMILARITY,
                        help=f"Minimum k-mer similarity to label a record (default {DEFAULT_MIN_SIMILARITY})")
    parser.add_argument("--min-margin", type=float, default=DEFAULT_MIN_MARGIN,
                        help="Minimum similarity difference between the two components "
                             f"(default {DEFAULT_MIN_MARGIN})")
    quality = parser.add_argument_group(
        "quality filters (off by default; excluded records are listed in excluded_records.csv)")
    quality.add_argument("--exclude-unverified", action="store_true",
                         help="Exclude records whose GenBank_Title is marked UNVERIFIED")
    quality.add_argument("--host", action="append", default=None, metavar="NAME",
                         help="Keep only records from this host (exact Host value, case-insensitive). "
                              "Repeat for several hosts. Records without a host are excluded")
    quality.add_argument("--max-ambiguous", type=fraction, default=None, metavar="FRACTION",
                         help="Exclude records whose share of non-ACGT bases exceeds FRACTION "
                              "(e.g. 0.01). Needs --fasta")
    quality.add_argument("--max-length-deviation", type=fraction, default=None, metavar="FRACTION",
                         help="Exclude records whose length differs from the median of their "
                              "component by more than FRACTION (e.g. 0.10). Needs --fasta")
    parser.add_argument("--options-file", type=Path, default=None,
                        help=f"File of options, one per line. Default: {OPTIONS_FILE} next to the "
                             "input CSV, if present")
    return parser


if __name__ == "__main__":
    run(parse_args())