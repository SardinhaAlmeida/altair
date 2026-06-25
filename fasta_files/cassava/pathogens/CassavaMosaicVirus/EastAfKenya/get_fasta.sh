#!/usr/bin/env bash

set -u
set -o pipefail

if [[ "$#" -lt 5 || "$#" -gt 6 ]]; then
    echo "Usage: $0 <pairs.csv> <unpaired_A.csv> <unpaired_B.csv> <output.fasta> <email> [api_key]" >&2
    exit 1
fi

pairs_csv="$1"
unpaired_a_csv="$2"
unpaired_b_csv="$3"
output_file="$4"
email="$5"
api_key="${6:-}"

base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

output_stem="${output_file%.*}"
failed_file="${output_stem}_failed_accessions.txt"
log_file="${output_stem}_download.log"
tmp_dir="$(mktemp -d)"

trap 'rm -rf "$tmp_dir"' EXIT

rm -f "$output_file" "$failed_file" "$log_file"

if ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: curl is not installed or not available in PATH." >&2
    exit 1
fi

check_csv_exists() {
    local input_csv="$1"

    if [[ ! -f "$input_csv" ]]; then
        echo "ERROR: CSV file not found: $input_csv" >&2
        exit 1
    fi
}

count_accessions_in_csv() {
    local input_csv="$1"

    awk -F',' '
        NR > 1 {
            gsub(/\r/, "", $1)
            gsub(/^[ \t]+|[ \t]+$/, "", $1)
            if ($1 != "") count++
        }
        END { print count + 0 }
    ' "$input_csv"
}

download_accession() {
    local id="$1"
    local tmp_fasta="$tmp_dir/${id}.fasta"
    local tmp_err="$tmp_dir/${id}.err"

    echo "Downloading $id ..." | tee -a "$log_file"

    local success=0

    for attempt in 1 2 3 4 5 6 7 8; do
        rm -f "$tmp_fasta" "$tmp_err"

        curl_args=(
            --http1.1
            --silent
            --show-error
            --location
            --fail-with-body
            --connect-timeout 20
            --max-time 120
            --retry 5
            --retry-delay 3
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
            curl_args+=(--data-urlencode "api_key=$api_key")
        fi

        if curl "${curl_args[@]}" -o "$tmp_fasta" 2> "$tmp_err"; then
            if grep -q "^>" "$tmp_fasta"; then
                header_count=$(grep -c "^>" "$tmp_fasta")

                if [[ "$header_count" -eq 1 ]]; then
                    cat "$tmp_fasta" >> "$output_file"
                    echo "" >> "$output_file"
                    success=1
                    echo "OK: $id" >> "$log_file"
                    break
                else
                    echo "WARNING: $id returned $header_count FASTA headers." | tee -a "$log_file"
                fi
            else
                echo "WARNING: no FASTA header for $id on attempt $attempt." | tee -a "$log_file"
                cat "$tmp_fasta" >> "$log_file"
                echo "" >> "$log_file"
            fi
        else
            echo "WARNING: curl failed for $id on attempt $attempt." | tee -a "$log_file"
            cat "$tmp_err" | tee -a "$log_file"
        fi

        sleep $((attempt * 3))
    done

    if [[ "$success" -eq 0 ]]; then
        echo "$id" >> "$failed_file"
        echo "FAILED: $id" | tee -a "$log_file"
        return 1
    fi

    return 0
}

download_from_csv() {
    local input_csv="$1"
    local label="$2"

    echo "" | tee -a "$log_file"
    echo "=== Downloading accessions from $label: $input_csv ===" | tee -a "$log_file"

    while IFS=, read -r accession _ || [[ -n "$accession" ]]; do
        accession="${accession//$'\r'/}"
        accession="$(echo "$accession" | xargs)"

        [[ -z "$accession" ]] && continue

        download_accession "$accession"

        # NCBI recommends limiting E-utilities requests without an API key.
        sleep 0.4

    done < <(tail -n +2 "$input_csv")
}

check_csv_exists "$pairs_csv"
check_csv_exists "$unpaired_a_csv"
check_csv_exists "$unpaired_b_csv"

expected_pairs=$(count_accessions_in_csv "$pairs_csv")
expected_unpaired_a=$(count_accessions_in_csv "$unpaired_a_csv")
expected_unpaired_b=$(count_accessions_in_csv "$unpaired_b_csv")
expected_count=$((expected_pairs + expected_unpaired_a + expected_unpaired_b))

echo "Download started: $(date)" | tee -a "$log_file"
echo "Pairs CSV: $pairs_csv" | tee -a "$log_file"
echo "Unpaired DNA-A CSV: $unpaired_a_csv" | tee -a "$log_file"
echo "Unpaired DNA-B CSV: $unpaired_b_csv" | tee -a "$log_file"
echo "Output FASTA: $output_file" | tee -a "$log_file"
echo "curl version: $(curl --version | head -1)" | tee -a "$log_file"
echo "" | tee -a "$log_file"

echo "Expected paired records: $expected_pairs" | tee -a "$log_file"
echo "Expected unpaired DNA-A records: $expected_unpaired_a" | tee -a "$log_file"
echo "Expected unpaired DNA-B records: $expected_unpaired_b" | tee -a "$log_file"
echo "Expected total records: $expected_count" | tee -a "$log_file"

if [[ "$expected_count" -eq 0 ]]; then
    echo "ERROR: no accessions found in the input CSV files." | tee -a "$log_file"
    exit 1
fi

download_from_csv "$pairs_csv" "paired DNA-A/DNA-B sequences"
download_from_csv "$unpaired_a_csv" "unpaired DNA-A sequences"
download_from_csv "$unpaired_b_csv" "unpaired DNA-B sequences"

actual_count=0

if [[ -f "$output_file" ]]; then
    actual_count=$(grep -c "^>" "$output_file")
fi

echo "" | tee -a "$log_file"
echo "Download finished: $(date)" | tee -a "$log_file"
echo "Expected records: $expected_count" | tee -a "$log_file"
echo "Downloaded records: $actual_count" | tee -a "$log_file"

if [[ "$actual_count" -ne "$expected_count" ]]; then
    echo "ERROR: incomplete download." | tee -a "$log_file"
    echo "Expected $expected_count records but downloaded $actual_count." | tee -a "$log_file"

    if [[ -s "$failed_file" ]]; then
        echo "Failed accessions saved in: $failed_file" | tee -a "$log_file"
    fi

    exit 1
fi

if [[ -s "$failed_file" ]]; then
    echo "ERROR: some accessions failed." | tee -a "$log_file"
    echo "Failed accessions saved in: $failed_file" | tee -a "$log_file"
    exit 1
fi

if grep -qiE "error|invalid|not found|cannot retrieve" "$output_file"; then
    echo "ERROR: output FASTA may contain an error message." | tee -a "$log_file"
    exit 1
fi

echo "SUCCESS: all records downloaded correctly." | tee -a "$log_file"
echo "Output FASTA: $output_file" | tee -a "$log_file"