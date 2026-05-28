#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Usage:
#   ./run_all_variants.sh HOST_FASTA PATHOGEN_ROOT [KMIN] [KMAX]
#
# Example:
#   ./run_all_variants.sh \
#     ../../fasta_files/cassava/hosts/tme7_1/tme7_phase1.fasta \
#     ../../fasta_files/cassava/pathogens/CassavaMosaicVirus \
#     11 13

HOST="$1"          # e.g. ../../fasta_files/.../hosts/tme7_1/tme7_phase1.fasta
PROOT="$2"         # e.g. ../../fasta_files/.../pathogens/CassavaMosaicVirus
KMIN="${3:-10}"
KMAX="${4:-13}"

# ----------------------------
# 1) Run run_raw.sh for ALL .fasta in ALL variant folders
# ----------------------------
echo "[1/2] Running AltaiR for all pathogen FASTA files..."

for variant_dir in "$PROOT"/*; do
  # Only folders inside CassavaMosaicVirus (EastAf, Indian, ...)
  [ -d "$variant_dir" ] || continue

  echo "  - Variant folder: $(basename "$variant_dir")"

  # All FASTA-like files inside this variant folder
  # for fasta in "$variant_dir"/*_joined.fasta "$variant_dir"/*_joined.fa "$variant_dir"/*_joined.fna; do
  #   [ -e "$fasta" ] || continue
  #   echo "      > Processing pathogen FASTA: $(basename "$fasta")"
  #   ./run_raw.sh "$HOST" "$fasta" "$KMIN" "$KMAX"
  # done
  for fasta in "$variant_dir"/*_joined.fasta "$variant_dir"/*_joined.fa "$variant_dir"/*_joined.fna; do
    [ -e "$fasta" ] || continue

    # Skip empty files (0 bytes)
    if [ ! -s "$fasta" ]; then
      echo "      ! Skipping empty joined FASTA: $(basename "$fasta")"
      continue
    fi

    # Skip files that have no FASTA records (no '>' headers)
    if ! grep -q '^>' "$fasta"; then
      echo "      ! Skipping joined FASTA with no records: $(basename "$fasta")"
      continue
    fi

    echo "      > Processing pathogen FASTA: $(basename "$fasta")"
    ./run_raw.sh "$HOST" "$fasta" "$KMIN" "$KMAX"
  done

done

echo "[1/2] AltaiR runs finished."

# ----------------------------
# 2) Run GCContent.py on all *.fasta-k*.eg under host_dir/tme7_1/*
# ----------------------------
echo "[2/2] Running GCContent on generated *.fasta-k*.eg files..."

host_dir="$(dirname "$HOST")"   # .../hosts/tme7_1

for variant_dir in "$host_dir"/*; do
  [ -d "$variant_dir" ] || continue

  eg_files=( "$variant_dir"/*_joined.fasta-k*.eg )
  if [ "${#eg_files[@]}" -eq 0 ]; then
    continue
  fi

  echo "  - Variant host folder: $(basename "$variant_dir")"
  echo "      > Found ${#eg_files[@]} _joined.fasta-k*.eg file(s)"

  # Run GCContent.py:
  #   -i with all .eg files in this variant folder
  #   -p using the same folder for plots/output
  python GCcontent.py -i "${eg_files[@]}" -p "$variant_dir"/results
done

echo "Complete"
