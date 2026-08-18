#!/usr/bin/env bash

set -euo pipefail

accession="JN555585.1"
output_file="HP.fasta"
email="YOUR_EMAIL@example.com"

curl \
    --fail \
    --silent \
    --show-error \
    --location \
    --retry 3 \
    --get "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
    --data-urlencode "db=nuccore" \
    --data-urlencode "id=${accession}" \
    --data-urlencode "rettype=fasta" \
    --data-urlencode "retmode=text" \
    --data-urlencode "tool=altair-thesis" \
    --data-urlencode "email=${email}" \
    --output "${output_file}"

if [[ ! -s "${output_file}" ]] || ! grep -q '^>' "${output_file}"; then
    echo "ERROR: NCBI did not return a valid FASTA file." >&2
    exit 1
fi

echo "Downloaded ${accession} to ${output_file}"
