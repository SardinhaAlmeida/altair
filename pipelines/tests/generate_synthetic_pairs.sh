#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GTO_BIN="${SCRIPT_DIR}/../../../gto/bin"
ALTAIR="${SCRIPT_DIR}/AltaiR"
OUTPUT_DIR="${SCRIPT_DIR}/synthetic_pairs"
HOST_SIZE=50000000
PATHOGEN_SIZE=10000
PATHOGEN_REPLICATES=3
MIN_K=10
MAX_K=15
RUN_ALTAIR=true
declare -A K_RANGES=()

host_percentages=(80 60 50 40 20 60 40)
pathogen_percentages=(20 40 50 60 80 60 40)

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Generate H80/P20, H60/P40, H50/P50, H40/P60, and H20/P80
A+T composition pairs. Each pair has one host and three independently
seeded pathogens. The seeds are recorded in seeds.tsv. Each pathogen has a
p1, p2, or p3 directory containing its seed record, raw results, and CG plot.

Options:
  --gto-bin DIR       GTO binary directory (default: ${GTO_BIN})
  --altair PATH       AltaiR executable (default: ${ALTAIR})
  --output-dir DIR    Output directory (default: ${OUTPUT_DIR})
  --host-size N       Host length (default: ${HOST_SIZE})
  --pathogen-size N   Pathogen length (default: ${PATHOGEN_SIZE})
  --pathogens N       Pathogens per pair (default: ${PATHOGEN_REPLICATES})
  --min-k N           AltaiR minimum context (default: ${MIN_K})
  --max-k N           AltaiR maximum context (default: ${MAX_K})
  --k-range SPEC      Override one pair as PAIR:MIN:MAX (repeatable)
  --generate-only     Generate FASTA files without running AltaiR
  -h, --help          Show this help

Pairs without a --k-range override use --min-k and --max-k.
EOF
}

while (($#)); do
  case "$1" in
    --gto-bin) GTO_BIN="$2"; shift 2 ;;
    --altair) ALTAIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --host-size) HOST_SIZE="$2"; shift 2 ;;
    --pathogen-size) PATHOGEN_SIZE="$2"; shift 2 ;;
    --pathogens) PATHOGEN_REPLICATES="$2"; shift 2 ;;
    --min-k) MIN_K="$2"; shift 2 ;;
    --max-k) MAX_K="$2"; shift 2 ;;
    --k-range)
      (($# >= 2)) || { echo "--k-range requires PAIR:MIN:MAX" >&2; exit 2; }
      IFS=: read -r range_pair range_min range_max range_extra <<< "$2"
      if [[ -z "$range_pair" || -z "$range_min" || -z "$range_max" ||
            -n "$range_extra" ]]; then
        echo "Invalid --k-range '$2'; expected PAIR:MIN:MAX" >&2
        exit 2
      fi
      if [[ -n "${K_RANGES[$range_pair]+x}" ]]; then
        echo "Duplicate --k-range for ${range_pair}" >&2
        exit 2
      fi
      K_RANGES["$range_pair"]="${range_min}:${range_max}"
      shift 2
      ;;
    --generate-only) RUN_ALTAIR=false; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

is_supported_pair() {
  local requested_pair="$1" pair_index supported_pair
  for pair_index in "${!host_percentages[@]}"; do
    supported_pair="H${host_percentages[$pair_index]}_P${pathogen_percentages[$pair_index]}"
    [[ "$requested_pair" == "$supported_pair" ]] && return 0
  done
  return 1
}

# Validate every override before generating any large FASTA file.
for range_pair in "${!K_RANGES[@]}"; do
  is_supported_pair "$range_pair" || {
    echo "Unknown pair in --k-range: ${range_pair}" >&2
    exit 2
  }
  IFS=: read -r range_min range_max <<< "${K_RANGES[$range_pair]}"
  for value in "$range_min" "$range_max"; do
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || {
      echo "Invalid k range for ${range_pair}: ${K_RANGES[$range_pair]}" >&2
      exit 2
    }
  done
  ((range_min <= range_max)) || {
    echo "Minimum k exceeds maximum k for ${range_pair}: ${K_RANGES[$range_pair]}" >&2
    exit 2
  }
done

for value in "$HOST_SIZE" "$PATHOGEN_SIZE" "$PATHOGEN_REPLICATES" "$MIN_K" "$MAX_K"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] ||
    { echo "All numeric options must be positive integers." >&2; exit 2; }
done
((MIN_K <= MAX_K)) || {
  echo "Minimum k exceeds maximum k: ${MIN_K}:${MAX_K}" >&2
  exit 2
}

GENERATOR="${GTO_BIN}/gto_genomic_gen_random_dna"
TO_FASTA="${GTO_BIN}/gto_fasta_from_seq"
for executable in "$GENERATOR" "$TO_FASTA"; do
  [[ -x "$executable" ]] || { echo "Executable not found: $executable" >&2; exit 1; }
done
if "$RUN_ALTAIR" && [[ ! -x "$ALTAIR" ]]; then
  echo "AltaiR executable not found: $ALTAIR" >&2
  exit 1
fi
if "$RUN_ALTAIR" && ! command -v gnuplot >/dev/null 2>&1; then
  echo "gnuplot is required to build the AltaiR CG plots." >&2
  echo "Install gnuplot, or use --generate-only to create only the FASTA files." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
GENERATOR="$(realpath "$GENERATOR")"
TO_FASTA="$(realpath "$TO_FASTA")"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
if "$RUN_ALTAIR"; then
  ALTAIR="$(realpath "$ALTAIR")"
fi
MANIFEST="${OUTPUT_DIR}/seeds.tsv"
printf "pair\trole\treplicate\tat_percent\tseed\tlength\tfile\n" > "$MANIFEST"

frequency_for_at() {
  case "$1" in
    80) echo "0.40,0.10,0.10,0.40" ;;
    60) echo "0.30,0.20,0.20,0.30" ;;
    50) echo "0.25,0.25,0.25,0.25" ;;
    40) echo "0.20,0.30,0.30,0.20" ;;
    20) echo "0.10,0.40,0.40,0.10" ;;
    *) echo "Unsupported A+T percentage: $1" >&2; return 1 ;;
  esac
}

generate_fasta() {
  local pair="$1" role="$2" replicate="$3" at="$4"
  local size="$5" seed="$6" output="$7" frequency
  frequency="$(frequency_for_at "$at")"

  echo "Generating ${pair} ${role}${replicate}: length=${size}, A+T=${at}%, seed=${seed}"
  "$GENERATOR" -s "$seed" -n "$size" -f "$frequency" |
    "$TO_FASTA" -n "${pair}_${role}${replicate}_seed${seed}" > "$output"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$pair" "$role" "$replicate" "$at" "$seed" "$size" "$output" >> "$MANIFEST"
}

for pair_index in "${!host_percentages[@]}"; do
  host_at="${host_percentages[$pair_index]}"
  pathogen_at="${pathogen_percentages[$pair_index]}"
  pair_name="H${host_at}_P${pathogen_at}"
  pair_min_k="$MIN_K"
  pair_max_k="$MAX_K"
  if [[ -n "${K_RANGES[$pair_name]+x}" ]]; then
    IFS=: read -r pair_min_k pair_max_k <<< "${K_RANGES[$pair_name]}"
  fi
  pair_dir="${OUTPUT_DIR}/${pair_name}"
  host_seed="$((11 + pair_index))"
  host_file="${pair_dir}/syn_H_${host_at}_seed${host_seed}.fasta"
  mkdir -p "$pair_dir"

  generate_fasta "$pair_name" host 1 "$host_at" \
    "$HOST_SIZE" "$host_seed" "$host_file"

  for ((replicate = 1; replicate <= PATHOGEN_REPLICATES; replicate++)); do
    pathogen_seed="$((101 + pair_index * PATHOGEN_REPLICATES + replicate - 1))"
    result_dir="${pair_dir}/p${replicate}"
    pathogen_file="${result_dir}/syn_P_${pathogen_at}_r${replicate}_seed${pathogen_seed}.fasta"
    seed_file="${result_dir}/seeds.tsv"
    mkdir -p "$result_dir"

    generate_fasta "$pair_name" pathogen "$replicate" "$pathogen_at" \
      "$PATHOGEN_SIZE" "$pathogen_seed" "$pathogen_file"

    printf "role\tat_percent\tseed\tlength\tfile\n" > "$seed_file"
    printf "host\t%s\t%s\t%s\t%s\n" \
      "$host_at" "$host_seed" "$HOST_SIZE" "$host_file" >> "$seed_file"
    printf "pathogen\t%s\t%s\t%s\t%s\n" \
      "$pathogen_at" "$pathogen_seed" "$PATHOGEN_SIZE" "$pathogen_file" >> "$seed_file"

    if "$RUN_ALTAIR"; then
      echo "Running AltaiR: ${pair_name}, pathogen replicate ${replicate}, k=${pair_min_k}:${pair_max_k}"
      (
        cd "$result_dir"
        "$ALTAIR" raw -v -p -min "$pair_min_k" -max "$pair_max_k" \
          "$host_file" "$pathogen_file"

        # AltaiR currently creates CGplot.sh; accept GCplot.sh as well in case
        # the executable's naming is corrected in a future version.
        if [[ -f CGplot.sh ]]; then
          echo "Building plot: ${result_dir}/CGplot.pdf"
          sh CGplot.sh
          [[ -s CGplot.pdf ]] || {
            echo "CGplot.sh did not produce a non-empty CGplot.pdf in: ${result_dir}" >&2
            exit 1
          }
        elif [[ -f GCplot.sh ]]; then
          echo "Building plot: ${result_dir}/GCplot.pdf"
          sh GCplot.sh
          [[ -s GCplot.pdf ]] || {
            echo "GCplot.sh did not produce a non-empty GCplot.pdf in: ${result_dir}" >&2
            exit 1
          }
        else
          echo "AltaiR did not create CGplot.sh or GCplot.sh in: ${result_dir}" >&2
          exit 1
        fi
      )
    fi
  done
done

echo "Done. Results: ${OUTPUT_DIR}"
echo "Seed manifest: ${MANIFEST}"
