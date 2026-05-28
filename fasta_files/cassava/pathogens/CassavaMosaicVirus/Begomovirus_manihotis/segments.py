import pandas as pd
import re
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Create DNA-A/DNA-B pairs from a sequences metadata CSV."
)

parser.add_argument(
    "input_csv",
    help="Input CSV file, e.g. sequences.csv"
)

parser.add_argument(
    "-o", "--out-dir",
    default="paired_output",
    help="Output directory. Default: paired_output"
)

args = parser.parse_args()

INPUT_CSV = Path(args.input_csv)

OUT_DIR = INPUT_CSV.parent / "paired_output"
OUT_DIR.mkdir(exist_ok=True)

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

df = pd.read_csv(INPUT_CSV, dtype=str).fillna("")

missing = [c for c in EXPECTED_COLUMNS if c not in df.columns]
if missing:
    print("Warning: missing columns:", missing)

def get_col(row, col):
    return row[col] if col in row.index else ""

def clean_text(x: str) -> str:
    x = str(x).strip().lower()
    x = re.sub(r"\s+", " ", x)
    return x

def classify_component(segment: str, title: str) -> str:
    text = f"{segment} {title}".lower()

    if re.search(r"\b(dna[\s\-]?a|component a|segment a|dna1)\b", text):
        return "A"
    if re.search(r"\b(dna[\s\-]?b|component b|segment b|dna2)\b", text):
        return "B"
    return "UNKNOWN"

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

# keep only complete records
if "Nuc_Completeness" in df.columns:
    df = df[df["Nuc_Completeness"].str.lower() == "complete"].copy()

df["Component"] = df.apply(
    lambda r: classify_component(get_col(r, "Segment"), get_col(r, "GenBank_Title")),
    axis=1
)

df["Isolate_Clean"] = df["Isolate"].apply(clean_text) if "Isolate" in df.columns else ""
df["Isolate_Pair_Core"] = df["Isolate"].apply(normalize_isolate_for_pairing) if "Isolate" in df.columns else ""

def build_pair_key(row) -> str:
    parts = [clean_text(get_col(row, col)) for col in PAIR_FIELDS if col in row.index]

    # use isolate core if available
    parts.append(get_col(row, "Isolate_Pair_Core"))

    return " || ".join(parts)

df["Pair_Key"] = df.apply(build_pair_key, axis=1)

# debug file
df.to_csv(OUT_DIR / "metadata_with_keys.csv", index=False)

unknown = df[df["Component"] == "UNKNOWN"].copy()
unknown.drop(columns=["Component", "Isolate_Clean", "Isolate_Pair_Core", "Pair_Key"], errors="ignore") \
       .to_csv(OUT_DIR / "unknown_component.csv", index=False)

known = df[df["Component"].isin(["A", "B"])].copy()

pairs_rows = []
ambiguous_rows = []
unpaired_rows = []

output_columns = [c for c in EXPECTED_COLUMNS if c in df.columns]

for pair_key, group in known.groupby("Pair_Key", dropna=False):
    a_rows = group[group["Component"] == "A"]
    b_rows = group[group["Component"] == "B"]

    if len(a_rows) == 1 and len(b_rows) == 1:
        a = a_rows.iloc[0]
        b = b_rows.iloc[0]

        pairs_rows.append({col: get_col(a, col) for col in output_columns})
        pairs_rows.append({col: get_col(b, col) for col in output_columns})

    elif len(a_rows) == 0 or len(b_rows) == 0:
        for _, row in group.iterrows():
            unpaired_rows.append({col: get_col(row, col) for col in output_columns})
    else:
        for _, row in group.iterrows():
            ambiguous_rows.append({col: get_col(row, col) for col in output_columns})

pairs_df = pd.DataFrame(pairs_rows, columns=output_columns)
ambiguous_df = pd.DataFrame(ambiguous_rows, columns=output_columns)
unpaired_df = pd.DataFrame(unpaired_rows, columns=output_columns)

pairs_df.to_csv(OUT_DIR / "pairs.csv", index=False)
ambiguous_df.to_csv(OUT_DIR / "ambiguous.csv", index=False)
unpaired_df.to_csv(OUT_DIR / "unpaired.csv", index=False)

unpaired_a_df = unpaired_df[unpaired_df["Segment"].str.contains(
    r"(?i)\b(dna[\s\-]?a|component a|segment a|dna1)\b",
    na=False,
    regex=True
)].copy()

unpaired_b_df = unpaired_df[unpaired_df["Segment"].str.contains(
    r"(?i)\b(dna[\s\-]?b|component b|segment b|dna2)\b",
    na=False,
    regex=True
)].copy()

# fallback using title when segment is messy or empty
segment_mask = unpaired_df["Segment"].fillna("").str.strip() != ""

unpaired_a_title_df = unpaired_df[
    ~segment_mask &
    unpaired_df["GenBank_Title"].str.contains(
        r"(?i)\b(dna[\s\-]?a|component a|segment a|dna1)\b",
        na=False,
        regex=True
    )
].copy()

unpaired_b_title_df = unpaired_df[
    ~segment_mask &
    unpaired_df["GenBank_Title"].str.contains(
        r"(?i)\b(dna[\s\-]?b|component b|segment b|dna2)\b",
        na=False,
        regex=True
    )
].copy()

unpaired_a_df = pd.concat([unpaired_a_df, unpaired_a_title_df], ignore_index=True).drop_duplicates()
unpaired_b_df = pd.concat([unpaired_b_df, unpaired_b_title_df], ignore_index=True).drop_duplicates()

unpaired_a_df.to_csv(OUT_DIR / "unpaired_A.csv", index=False)
unpaired_b_df.to_csv(OUT_DIR / "unpaired_B.csv", index=False)

print("Done.")
print(f"Pairs written: {len(pairs_df)} rows ({len(pairs_df)//2} A/B pairs)")
print(f"Ambiguous rows: {len(ambiguous_df)}")
print(f"Unpaired rows: {len(unpaired_df)}")
print(f"Unpaired A rows: {len(unpaired_a_df)}")
print(f"Unpaired B rows: {len(unpaired_b_df)}")
print(f"Unknown component rows: {len(unknown)}")