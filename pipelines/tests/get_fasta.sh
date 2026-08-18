#!/usr/bin/env bash

set -u
set -o pipefail

usage() {
    echo "Usage: $0 <pathogen_directory> <output_fasta_name> <email> [api_key]" >&2
    echo "The directory may contain either sequences.csv (complete genomes) or" >&2
    echo "paired_output/{pairs,unpaired_A,unpaired_B}.csv (segmented genomes)." >&2
    echo "Example:" >&2
    echo "  $0 ../../fasta_files/cassava/pathogens/CassavaMosaicVirus/EastAfKenya EA_CMV.fasta user@example.com" >&2
}

if [[ "$#" -lt 3 || "$#" -gt 4 ]]; then
    usage
    exit 1
fi

pathogen_dir_input="$1"
output_name="$2"
email="$3"
api_key="${4:-}"

# ---------------------------------------------------------------------------
# Validate command-line arguments
# ---------------------------------------------------------------------------

if [[ ! -d "$pathogen_dir_input" ]]; then
    echo "ERROR: pathogen directory not found: $pathogen_dir_input" >&2
    exit 1
fi

if [[ "$output_name" == */* ]]; then
    echo "ERROR: output_fasta_name must be a filename, not a path." >&2
    echo "Example: EA_CMV.fasta" >&2
    exit 1
fi

case "$output_name" in
    *.fasta|*.fa|*.fna)
        ;;
    *)
        echo "ERROR: output filename must end in .fasta, .fa, or .fna." >&2
        exit 1
        ;;
esac

if [[ -z "$email" ]]; then
    echo "ERROR: email cannot be empty." >&2
    exit 1
fi

if ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: curl is not installed or not available in PATH." >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is not installed or not available in PATH." >&2
    exit 1
fi

# Resolve the pathogen directory to an absolute path.
pathogen_dir="$(cd "$pathogen_dir_input" && pwd -P)"

paired_dir="$pathogen_dir/paired_output"

pairs_csv="$paired_dir/pairs.csv"
unpaired_a_csv="$paired_dir/unpaired_A.csv"
unpaired_b_csv="$paired_dir/unpaired_B.csv"
sequences_csv="$pathogen_dir/sequences.csv"

output_file="$pathogen_dir/$output_name"
output_base="${output_name%.*}"

failed_file="$pathogen_dir/${output_base}_failed_accessions.txt"
log_file="$pathogen_dir/${output_base}_download.log"

base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# ---------------------------------------------------------------------------
# Validate input CSV files
# ---------------------------------------------------------------------------

check_csv_exists() {
    local input_csv="$1"

    if [[ ! -f "$input_csv" ]]; then
        echo "ERROR: CSV file not found: $input_csv" >&2
        exit 1
    fi

    if [[ ! -r "$input_csv" ]]; then
        echo "ERROR: CSV file is not readable: $input_csv" >&2
        exit 1
    fi
}

input_mode=""
input_csvs=()

if [[ -f "$pairs_csv" && -f "$unpaired_a_csv" && -f "$unpaired_b_csv" ]]; then
    input_mode="segmented"
    input_csvs=("$pairs_csv" "$unpaired_a_csv" "$unpaired_b_csv")
elif [[ -f "$sequences_csv" ]]; then
    input_mode="complete"
    input_csvs=("$sequences_csv")
else
    echo "ERROR: no supported CSV input was found in: $pathogen_dir" >&2
    echo "Expected either $sequences_csv or all three CSV files in $paired_dir." >&2
    exit 1
fi

for input_csv in "${input_csvs[@]}"; do
    check_csv_exists "$input_csv"
done

# ---------------------------------------------------------------------------
# Create temporary working files
#
# The temporary directory is created inside the pathogen directory. Therefore,
# the completed temporary FASTA and the final FASTA are normally on the same
# filesystem, allowing the final mv operation to be atomic.
# ---------------------------------------------------------------------------

tmp_dir="$(mktemp -d "$pathogen_dir/.get_fasta.XXXXXX")"

tmp_output="$tmp_dir/downloaded.fasta"
manifest_file="$tmp_dir/accessions.txt"
counts_file="$tmp_dir/counts.txt"

cleanup() {
    rm -rf "$tmp_dir"
}

trap cleanup EXIT

: > "$tmp_output"

# Remove old diagnostic files, but do not remove the existing FASTA.
# The existing FASTA is only replaced after a fully successful download.
rm -f "$failed_file" "$log_file"

# ---------------------------------------------------------------------------
# Build an accession manifest
#
# Python's csv module is used so that:
#   - the Accession field is identified by column name;
#   - quoted CSV fields are handled correctly;
#   - duplicate accessions can be detected before downloading.
# ---------------------------------------------------------------------------

if ! python3 - \
    "$manifest_file" \
    "$counts_file" \
    "${input_csvs[@]}" \
    2> >(tee -a "$log_file" >&2) <<'PY'
import csv
import sys
from pathlib import Path

manifest_file = Path(sys.argv[1])
counts_file = Path(sys.argv[2])
csv_paths = [Path(value) for value in sys.argv[3:]]

sources = [(csv_path.stem, csv_path) for csv_path in csv_paths]

seen = {}
duplicates = []
counts = {label: 0 for label, _ in sources}
ordered_accessions = []

for label, csv_path in sources:
    with csv_path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames or "Accession" not in reader.fieldnames:
            raise SystemExit(
                f"ERROR: missing 'Accession' column in {csv_path}"
            )

        for line_number, row in enumerate(reader, start=2):
            accession = (row.get("Accession") or "").strip()

            if not accession:
                continue

            counts[label] += 1

            if accession in seen:
                previous_path, previous_line = seen[accession]

                duplicates.append(
                    f"{accession}: "
                    f"{previous_path}:{previous_line} and "
                    f"{csv_path}:{line_number}"
                )
                continue

            seen[accession] = (csv_path, line_number)
            ordered_accessions.append(accession)

if duplicates:
    print(
        "ERROR: duplicate accessions found in the input CSV files:",
        file=sys.stderr
    )

    for duplicate in duplicates:
        print(f"  - {duplicate}", file=sys.stderr)

    raise SystemExit(1)

with manifest_file.open(
    "w",
    encoding="utf-8",
    newline="\n"
) as handle:
    for accession in ordered_accessions:
        handle.write(accession + "\n")

with counts_file.open(
    "w",
    encoding="utf-8",
    newline="\n"
) as handle:
    for label, _ in sources:
        handle.write(f"{label}={counts[label]}\n")

    handle.write(f"total={len(ordered_accessions)}\n")
PY
then
    echo "ERROR: could not create a valid accession manifest." \
        | tee -a "$log_file" >&2
    exit 1
fi

# The counts file contains simple shell assignments created internally by
# the Python code above.
# shellcheck disable=SC1090
source "$counts_file"

expected_count="${total:-0}"

if [[ "$expected_count" -eq 0 ]]; then
    echo "ERROR: no accessions found in the input CSV files." \
        | tee -a "$log_file" >&2
    exit 1
fi

# Use a shorter delay when an NCBI API key is supplied.
if [[ -n "$api_key" ]]; then
    request_delay="0.12"
else
    request_delay="0.40"
fi

# ---------------------------------------------------------------------------
# FASTA response validation
# ---------------------------------------------------------------------------

validate_fasta() {
    local fasta_file="$1"
    local accession="$2"
    local header_count
    local sequence

    header_count="$(grep -c '^>' "$fasta_file" || true)"

    if [[ "$header_count" -ne 1 ]]; then
        echo "WARNING: $accession returned $header_count FASTA headers." \
            | tee -a "$log_file"
        return 1
    fi

    sequence="$(
        grep -v '^>' "$fasta_file" \
            | tr -d '[:space:]' \
            | tr '[:lower:]' '[:upper:]'
    )"

    if [[ -z "$sequence" ]]; then
        echo "WARNING: $accession returned an empty sequence." \
            | tee -a "$log_file"
        return 1
    fi

    # Accepted DNA IUPAC symbols:
    # A C G T R Y S W K M B D H V N
    if [[ ! "$sequence" =~ ^[ACGTRYSWKMBDHVN]+$ ]]; then
        echo "WARNING: $accession returned invalid nucleotide symbols." \
            | tee -a "$log_file"
        return 1
    fi

    return 0
}

# ---------------------------------------------------------------------------
# Download one accession
# ---------------------------------------------------------------------------

download_accession() {
    local id="$1"
    local tmp_fasta="$tmp_dir/response.fasta"
    local tmp_err="$tmp_dir/curl.err"
    local success=0
    local attempt
    local -a curl_args

    echo "Downloading $id ..." | tee -a "$log_file"

    for attempt in 1 2 3 4 5; do
        rm -f "$tmp_fasta" "$tmp_err"

        curl_args=(
            --http1.1
            --silent
            --show-error
            --location
            --fail-with-body
            --connect-timeout 20
            --max-time 120
            --retry 2
            --retry-delay 2
            --retry-all-errors
            --get "$base_url"
            --data-urlencode "db=nuccore"
            --data-urlencode "id=$id"
            --data-urlencode "rettype=fasta"
            --data-urlencode "retmode=text"
            --data-urlencode "tool=altair-thesis"
            --data-urlencode "email=$email"
        )

        if [[ -n "$api_key" ]]; then
            curl_args+=(
                --data-urlencode "api_key=$api_key"
            )
        fi

        if curl "${curl_args[@]}" \
            -o "$tmp_fasta" \
            2> "$tmp_err"
        then
            if validate_fasta "$tmp_fasta" "$id"; then
                cat "$tmp_fasta" >> "$tmp_output"
                printf '\n' >> "$tmp_output"

                success=1

                echo "OK: $id" >> "$log_file"
                break
            fi

            echo "Response preview for $id (attempt $attempt):" \
                >> "$log_file"

            head -n 20 "$tmp_fasta" \
                >> "$log_file" \
                2>/dev/null || true

            printf '\n' >> "$log_file"
        else
            echo "WARNING: curl failed for $id on attempt $attempt." \
                | tee -a "$log_file"

            if [[ -s "$tmp_err" ]]; then
                tee -a "$log_file" < "$tmp_err"
            fi
        fi

        if [[ "$attempt" -lt 5 ]]; then
            sleep $((attempt * 2))
        fi
    done

    if [[ "$success" -eq 0 ]]; then
        echo "$id" >> "$failed_file"
        echo "FAILED: $id" | tee -a "$log_file"
        return 1
    fi

    return 0
}

# ---------------------------------------------------------------------------
# Log initial information
# ---------------------------------------------------------------------------

echo "Download started: $(date)" | tee -a "$log_file"
echo "Pathogen directory: $pathogen_dir" | tee -a "$log_file"
echo "Input mode: $input_mode" | tee -a "$log_file"
for input_csv in "${input_csvs[@]}"; do
    echo "Input CSV: $input_csv" | tee -a "$log_file"
done
echo "Output FASTA: $output_file" | tee -a "$log_file"
echo "curl version: $(curl --version | head -n 1)" | tee -a "$log_file"
echo "" | tee -a "$log_file"

echo "Expected unique total records: $expected_count" \
    | tee -a "$log_file"

echo "" | tee -a "$log_file"

# ---------------------------------------------------------------------------
# Download all unique accessions
# ---------------------------------------------------------------------------

while IFS= read -r accession || [[ -n "$accession" ]]; do
    [[ -z "$accession" ]] && continue

    # The script continues after an individual failure so that every failed
    # accession can be reported in the failure file.
    download_accession "$accession" || true

    sleep "$request_delay"
done < "$manifest_file"

actual_count="$(grep -c '^>' "$tmp_output" || true)"

# ---------------------------------------------------------------------------
# Validate the complete download
# ---------------------------------------------------------------------------

echo "" | tee -a "$log_file"
echo "Download finished: $(date)" | tee -a "$log_file"
echo "Expected records: $expected_count" | tee -a "$log_file"
echo "Downloaded records: $actual_count" | tee -a "$log_file"

if [[ "$actual_count" -ne "$expected_count" ]]; then
    echo "ERROR: incomplete download." \
        | tee -a "$log_file" >&2

    echo "Expected $expected_count records but downloaded $actual_count." \
        | tee -a "$log_file" >&2

    if [[ -s "$failed_file" ]]; then
        echo "Failed accessions saved in: $failed_file" \
            | tee -a "$log_file" >&2
    fi

    echo "The existing output FASTA, if present, was left unchanged." \
        | tee -a "$log_file" >&2

    exit 1
fi

if [[ -s "$failed_file" ]]; then
    echo "ERROR: some accessions failed." \
        | tee -a "$log_file" >&2

    echo "Failed accessions saved in: $failed_file" \
        | tee -a "$log_file" >&2

    echo "The existing output FASTA, if present, was left unchanged." \
        | tee -a "$log_file" >&2

    exit 1
fi

# ---------------------------------------------------------------------------
# Publish the completed FASTA
#
# This happens only after every accession has been successfully downloaded
# and validated. Until this point, an existing output FASTA is untouched.
# ---------------------------------------------------------------------------

mv -f "$tmp_output" "$output_file"

echo "SUCCESS: all records downloaded and validated correctly." \
    | tee -a "$log_file"

echo "Output FASTA: $output_file" \
    | tee -a "$log_file"
