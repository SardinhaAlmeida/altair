#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

HOST=$1         # e.g. ../../fasta_files/.../hosts/tme7_1/tme7_phase1.fasta
PATHOGEN=$2     # e.g. ../../fasta_files/.../pathogens/.../EastAf/EA_CMV.fasta
KMIN=${3:-10}
KMAX=${4:-13}

# Names and dirs
host_dir=$(dirname "$HOST")
hbase=$(basename "$HOST");      hbase="${hbase%.*}"

pdir=$(dirname "$PATHOGEN")     # .../EastAf
pgroup=$(basename "$pdir")      # EastAf
pbase=$(basename "$PATHOGEN");  pbase="${pbase%.*}"   # EA_CMV

# 1) Run AltaiR (writes outputs next to PATHOGEN)
AltaiR raw -min "$KMIN" -max "$KMAX" -f -p "$HOST" "$PATHOGEN"

# 2) Destination beside host, under the pathogen’s parent dir name
dest="$host_dir/$pgroup"
mkdir -p "$dest"

# delete previous .eg files in the destination
echo "[i] Cleaning previous files inside: $dest"
find "$dest" -type f -delete

# 3) Move outputs over (keep filenames, including .fasta part if present)
src_dir="$pdir"
moved_any=false
for f in "$src_dir/$pbase.fasta"-k*.eg "$src_dir/$pbase"-k*.eg; do
  [ -e "$f" ] || continue
  mv "$f" "$dest/"
  moved_any=true
done

if ! $moved_any; then
  echo "No output files found to move from: $src_dir"
  exit 1
fi

echo "Results in: $dest"
