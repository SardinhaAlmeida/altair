# #!/bin/bash

# input_csv="paired_output/pairs.csv"
# output_fasta="ACMV.fasta"

# rm -f "$output_fasta"

# tail -n +2 "$input_csv" | while IFS=, read -r accession _; do
#   accession="${accession//$'\r'/}"
#   [[ -z "$accession" ]] && continue
#   echo "Downloading $accession ..."
#   efetch -db nucleotide -format fasta -id "$accession" >> "$output_fasta"
# done

pairs_csv="paired_output/pairs.csv"
unpaired_a_csv="paired_output/unpaired_A.csv"
unpaired_b_csv="paired_output/unpaired_B.csv"

output_fasta="ACMV.fasta"

rm -f "$output_fasta"

download_from_csv() {
  local input_csv="$1"
  local label="$2"

  if [[ ! -f "$input_csv" ]]; then
    echo "[WARN] File not found: $input_csv"
    return
  fi

  echo
  echo "=== Downloading accessions from $label: $input_csv ==="

  tail -n +2 "$input_csv" | while IFS=, read -r accession _; do
    accession="${accession//$'\r'/}"
    accession="$(echo "$accession" | xargs)"

    [[ -z "$accession" ]] && continue

    echo "Downloading $accession ..."
    efetch -db nucleotide -format fasta -id "$accession" >> "$output_fasta"
  done
}

download_from_csv "$pairs_csv" "paired A/B sequences"
download_from_csv "$unpaired_a_csv" "unpaired DNA-A sequences"
download_from_csv "$unpaired_b_csv" "unpaired DNA-B sequences"

echo
echo "Done. All sequences saved to: $output_fasta"