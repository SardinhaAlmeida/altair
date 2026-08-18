#!/usr/bin/env bash

set -euo pipefail

accession="GCF_000203835.1"
output_file="Sc.fasta"
zip_file="Sc_NCBI.zip"
tmp_dir="Sc_NCBI"

rm -f "$output_file" "$zip_file"
rm -rf "$tmp_dir"

echo "Downloading ${accession} from NCBI..."

datasets download genome accession "$accession" \
    --include genome \
    --filename "$zip_file"

echo "Extracting genome..."

unzip -q "$zip_file" -d "$tmp_dir"

genome_file=$(find "$tmp_dir/ncbi_dataset/data/$accession" \
    -name "*_genomic.fna" \
    -type f | head -n 1)

if [[ -z "$genome_file" || ! -s "$genome_file" ]]; then
    echo "ERROR: NCBI did not return a genomic FASTA." >&2
    exit 1
fi

cp "$genome_file" "$output_file"

if ! grep -q '^>' "$output_file"; then
    echo "ERROR: Invalid FASTA file." >&2
    exit 1
fi

rm -rf "$tmp_dir" "$zip_file"

echo "Downloaded ${accession} to ${output_file}"