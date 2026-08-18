#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

HOST=$1         # e.g. ../../fasta_files/.../hosts/tme7_1/tme7_phase1.fasta
PATHOGEN=$2     # e.g. ../../fasta_files/.../pathogens/.../EastAf/EA_CMV.fasta
KMIN=${3:-10}
KMAX=${4:-13}
RESULTS_DIR=${5:-}
ANALYSIS_NAME=${6:-}

if ! command -v gnuplot >/dev/null 2>&1; then
  echo "gnuplot is required to generate the AltaiR CG plot." >&2
  exit 1
fi

# AltaiR creates some outputs in its current working directory. Resolve the
# inputs before changing directory so callers can continue to use relative
# paths.
HOST=$(realpath "$HOST")
PATHOGEN=$(realpath "$PATHOGEN")

# Names and dirs
host_dir=$(dirname "$HOST")
hbase=$(basename "$HOST");      hbase="${hbase%.*}"

pdir=$(dirname "$PATHOGEN")     # .../EastAf
pgroup=$(basename "$pdir")      # EastAf
pbase=$(basename "$PATHOGEN");  pbase="${pbase%.*}"   # EA_CMV
ANALYSIS_NAME=${ANALYSIS_NAME:-$hbase-$pbase}

# 1) Run AltaiR in an isolated directory. The RAW .eg files are written next
# to PATHOGEN; fixed-name plot files are written to the current directory.
work_dir=$(mktemp -d)
cleanup() {
  rm -rf -- "$work_dir"
}
trap cleanup EXIT

(
  cd "$work_dir"
  AltaiR raw -min "$KMIN" -max "$KMAX" -f -p "$HOST" "$PATHOGEN"
)

# 2) Destination beside host, under the pathogen’s parent dir name
dest="$host_dir/$pgroup"
mkdir -p "$dest"

# Delete previous files produced directly in the RAW destination, without
# deleting files stored in its results subdirectory.
echo "[i] Cleaning previous AltaiR files inside: $dest"
find "$dest" -maxdepth 1 -type f -delete

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

# 4) Keep AltaiR's supporting data and scripts with the RAW files.
for generated_file in CG-data.eg CGplot.sh AVG-data.eg AVGplot.sh; do
  if [ -e "$work_dir/$generated_file" ]; then
    mv "$work_dir/$generated_file" "$dest/"
  fi
done

if [ ! -s "$dest/CG-data.eg" ]; then
  echo "AltaiR did not generate a non-empty CG-data.eg file." >&2
  exit 1
fi

if [ ! -s "$dest/CGplot.sh" ]; then
  echo "AltaiR did not generate a non-empty CGplot.sh file." >&2
  exit 1
fi

# Use <RAW directory>/results when the caller does not provide a destination.
if [ -z "$RESULTS_DIR" ]; then
  RESULTS_DIR="$dest/results"
else
  mkdir -p "$RESULTS_DIR"
  RESULTS_DIR=$(cd "$RESULTS_DIR" && pwd -P)
fi
mkdir -p "$RESULTS_DIR"

# CGplot.sh expects CG-data.eg in its working directory and creates
# CGplot.pdf there. Move that PDF to its single final location afterwards.
echo "[i] Generating AltaiR CG plot..."
(
  cd "$dest"
  sh ./CGplot.sh
)

if [ ! -s "$dest/CGplot.pdf" ]; then
  echo "CGplot.sh did not generate a non-empty CGplot.pdf." >&2
  exit 1
fi

plot_path="$RESULTS_DIR/${ANALYSIS_NAME}_CGplot.pdf"
mv "$dest/CGplot.pdf" "$plot_path"

echo "RAW results in: $dest"
echo "CG plot in: $plot_path"
