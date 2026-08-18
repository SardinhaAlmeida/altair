#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

usage() {
    cat >&2 <<USAGE
Usage:
  $0 prepare <pathogen_directory> <output_prefix> <email> [api_key]

  $0 complete <pathogen_directory> <host_fasta> <analysis_name> <email> \
[kmin] [kmax] [api_key]

  $0 analyse <pathogen_fasta_or_directory> <host_fasta> <analysis_name> \
[kmin] [kmax]

  $0 <pathogen_directory> <host_fasta> <output_prefix> <email> \
[kmin] [kmax] [api_key]

Workflows:
  prepare
      (This command was made for Bipartite Begomoviruses to handle DNA-A and DNA-B components)
      Classify the metadata, download nucleotide sequences,
      and organise paired and unpaired DNA components.

  complete
      Download every accession from sequences.csv as an independent complete
      genome, then run the RAW analysis. This workflow does not classify or
      join DNA-A/DNA-B components.

  analyse
      Analyse an existing pathogen FASTA. The first argument may be:
        - a FASTA file; or
        - a directory containing the FASTA.

      If a directory is supplied, the script first looks for files named from
      <analysis_name>, including <analysis_name>_joined.fasta and
      <analysis_name>.fasta. If neither exists, the directory must contain
      exactly one non-empty .fasta, .fa, or .fna file.

  no workflow argument
      Run the complete preparation and analysis pipeline for a segmented
      pathogen dataset.

Examples:

  $0 complete \
     ../../fasta_files/virus_complete_genomes \
     ../../fasta_files/human/human.fna \
     VIRUS \
     user@example.com \
     10 \
     13

  $0 prepare \
     ../../fasta_files/cassava/pathogens/CassavaMosaicVirus/Begomovirus_manihotis \
     ACMV \
     user@example.com

  $0 analyse \
     ../../fasta_files/synthetic/synthetic_pathogens/60AT \
     ../../fasta_files/synthetic/synthetic_hosts/40AT/synthetic_h_40at.fasta \
     P60 \
     11 \
     13

  $0 analyse \
     ../../fasta_files/cassava/pathogens/CassavaMosaicVirus/Begomovirus_manihotis/ACMV_joined.fasta \
     ../../fasta_files/cassava/hosts/tmeb117/TMEB117.fa \
     ACMV \
     11 \
     13

  $0 \
     ../../fasta_files/cassava/pathogens/CassavaMosaicVirus/Begomovirus_manihotis \
     ../../fasta_files/cassava/hosts/tmeb117/TMEB117.fa \
     ACMV \
     user@example.com \
     11 \
     13
USAGE
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    local command_name="$1"

    if ! command -v "$command_name" >/dev/null 2>&1; then
        fail "required command is not available: $command_name"
    fi
}

require_script() {
    local script_path="$1"

    if [[ ! -f "$script_path" ]]; then
        fail "required pipeline script not found: $script_path"
    fi
}

resolve_directory() {
    local directory="$1"

    if [[ ! -d "$directory" ]]; then
        fail "directory not found: $directory"
    fi

    (
        cd "$directory"
        pwd -P
    )
}

resolve_file() {
    local file_path="$1"
    local file_dir
    local file_name

    if [[ ! -f "$file_path" ]]; then
        fail "file not found: $file_path"
    fi

    if [[ ! -s "$file_path" ]]; then
        fail "file is empty: $file_path"
    fi

    file_dir="$(
        cd "$(dirname "$file_path")"
        pwd -P
    )"

    file_name="$(basename "$file_path")"

    printf '%s/%s\n' "$file_dir" "$file_name"
}

resolve_pathogen_fasta() {
    local pathogen_input="$1"
    local analysis_name="$2"
    local pathogen_dir
    local candidate
    local -a candidates=()
    local -a nonempty_candidates=()

    # Direct FASTA input: use it immediately.
    if [[ -f "$pathogen_input" ]]; then
        resolve_file "$pathogen_input"
        return
    fi

    if [[ ! -d "$pathogen_input" ]]; then
        fail "pathogen FASTA or directory not found: $pathogen_input"
    fi

    pathogen_dir="$(resolve_directory "$pathogen_input")"

    # Prefer a file that explicitly matches the supplied analysis name.
    for candidate in \
        "$pathogen_dir/${analysis_name}_joined.fasta" \
        "$pathogen_dir/${analysis_name}_joined.fa" \
        "$pathogen_dir/${analysis_name}_joined.fna" \
        "$pathogen_dir/${analysis_name}.fasta" \
        "$pathogen_dir/${analysis_name}.fa" \
        "$pathogen_dir/${analysis_name}.fna"
    do
        if [[ -s "$candidate" ]]; then
            resolve_file "$candidate"
            return
        fi
    done

    # If no filename matches the analysis name, accept the directory only when
    # it contains exactly one non-empty FASTA file.
    shopt -s nullglob
    candidates=(
        "$pathogen_dir"/*.fasta
        "$pathogen_dir"/*.fa
        "$pathogen_dir"/*.fna
    )
    shopt -u nullglob

    for candidate in "${candidates[@]}"; do
        if [[ -s "$candidate" ]]; then
            nonempty_candidates+=("$candidate")
        fi
    done

    case "${#nonempty_candidates[@]}" in
        0)
            fail "no non-empty FASTA file found in pathogen directory: $pathogen_dir"
            ;;
        1)
            resolve_file "${nonempty_candidates[0]}"
            ;;
        *)
            echo "ERROR: several pathogen FASTA files were found in:" >&2
            echo "  $pathogen_dir" >&2
            echo "Pass the intended FASTA file directly to the analyse command:" >&2

            for candidate in "${nonempty_candidates[@]}"; do
                echo "  - $candidate" >&2
            done

            exit 1
            ;;
    esac
}

validate_output_prefix() {
    local output_prefix="$1"

    if [[ -z "$output_prefix" || "$output_prefix" == */* ]]; then
        fail "output prefix or analysis name must be a simple name such as ACMV, P60, or H40_P60"
    fi

    if [[ ! "$output_prefix" =~ ^[A-Za-z0-9._-]+$ ]]; then
        fail "output prefix or analysis name may only contain letters, numbers, dots, underscores, and hyphens"
    fi
}

validate_kmer_range() {
    local kmin="$1"
    local kmax="$2"

    if ! [[ "$kmin" =~ ^[1-9][0-9]*$ && "$kmax" =~ ^[1-9][0-9]*$ ]]; then
        fail "kmin and kmax must be positive integers"
    fi

    if (( kmin > kmax )); then
        fail "kmin cannot be greater than kmax"
    fi
}

print_separator() {
    echo "============================================================"
}

prepare_workflow() {
    local pathogen_dir="$1"
    local output_prefix="$2"
    local email="$3"
    local api_key="${4:-}"

    local metadata_csv="$pathogen_dir/sequences.csv"
    local paired_dir="$pathogen_dir/paired_output"
    local downloaded_fasta="$pathogen_dir/${output_prefix}.fasta"
    local joined_fasta="$pathogen_dir/${output_prefix}_joined.fasta"
    local required_csv

    if [[ -z "$email" ]]; then
        fail "email cannot be empty"
    fi

    if [[ ! -f "$metadata_csv" ]]; then
        fail "metadata CSV not found: $metadata_csv"
    fi

    require_command python3
    require_command curl

    require_script "$script_dir/segments.py"
    require_script "$script_dir/get_fasta.sh"
    require_script "$script_dir/JoinSegments2.py"

    print_separator
    echo "Pathogen data-preparation workflow"
    print_separator
    echo "Pathogen directory : $pathogen_dir"
    echo "Output prefix      : $output_prefix"
    echo "Metadata CSV       : $metadata_csv"
    echo "Downloaded FASTA   : $downloaded_fasta"
    echo "Organised FASTA    : $joined_fasta"
    print_separator
    echo

    echo "[1/3] Classifying viral genome components..."

    python3 "$script_dir/segments.py" \
        "$metadata_csv"

    echo "[1/3] Completed."
    echo

    echo "[2/3] Downloading and validating nucleotide sequences..."

    if [[ -n "$api_key" ]]; then
        bash "$script_dir/get_fasta.sh" \
            "$pathogen_dir" \
            "${output_prefix}.fasta" \
            "$email" \
            "$api_key"
    else
        bash "$script_dir/get_fasta.sh" \
            "$pathogen_dir" \
            "${output_prefix}.fasta" \
            "$email"
    fi

    echo "[2/3] Completed."
    echo

    for required_csv in \
        "$paired_dir/pairs.csv" \
        "$paired_dir/unpaired_A.csv" \
        "$paired_dir/unpaired_B.csv"
    do
        if [[ ! -f "$required_csv" ]]; then
            fail "expected classification output not found: $required_csv"
        fi
    done

    if [[ ! -s "$downloaded_fasta" ]]; then
        fail "downloaded FASTA not found or empty: $downloaded_fasta"
    fi

    echo "[3/3] Organising paired and unpaired DNA components..."

    python3 "$script_dir/JoinSegments2.py" \
        "$paired_dir/pairs.csv" \
        "$paired_dir/unpaired_A.csv" \
        "$paired_dir/unpaired_B.csv" \
        "$downloaded_fasta" \
        "$joined_fasta"

    if [[ ! -s "$joined_fasta" ]]; then
        fail "organised FASTA was not created or is empty: $joined_fasta"
    fi

    echo "[3/3] Completed."
    echo

    print_separator
    echo "Preparation workflow completed successfully"
    print_separator
    echo "Downloaded FASTA : $downloaded_fasta"
    echo "Organised FASTA  : $joined_fasta"
    print_separator
}

complete_workflow() {
    local pathogen_dir="$1"
    local host_fasta="$2"
    local analysis_name="$3"
    local email="$4"
    local kmin="$5"
    local kmax="$6"
    local api_key="${7:-}"
    local metadata_csv="$pathogen_dir/sequences.csv"
    local downloaded_fasta="$pathogen_dir/${analysis_name}.fasta"

    if [[ ! -f "$metadata_csv" ]]; then
        fail "metadata CSV not found: $metadata_csv"
    fi

    require_command python3
    require_command curl
    require_script "$script_dir/get_fasta.sh"

    print_separator
    echo "Complete-genome download workflow"
    print_separator
    echo "Pathogen directory : $pathogen_dir"
    echo "Metadata CSV       : $metadata_csv"
    echo "Downloaded FASTA   : $downloaded_fasta"
    print_separator
    echo

    # Ensure this mode reads sequences.csv rather than stale segmented output.
    if [[ -d "$pathogen_dir/paired_output" ]]; then
        fail "complete workflow requires a non-segmented dataset, but paired_output exists: $pathogen_dir/paired_output"
    fi

    if [[ -n "$api_key" ]]; then
        bash "$script_dir/get_fasta.sh" \
            "$pathogen_dir" "${analysis_name}.fasta" "$email" "$api_key"
    else
        bash "$script_dir/get_fasta.sh" \
            "$pathogen_dir" "${analysis_name}.fasta" "$email"
    fi

    analyse_workflow \
        "$downloaded_fasta" \
        "$host_fasta" \
        "$analysis_name" \
        "$kmin" \
        "$kmax"
}

analyse_workflow() {
    local pathogen_fasta="$1"
    local host_fasta="$2"
    local analysis_name="$3"
    local kmin="$4"
    local kmax="$5"

    local pathogen_filename
    local pathogen_stem
    local pathogen_group
    local host_dir
    local raw_dir
    local results_dir
    local -a raw_files

    pathogen_filename="$(basename "$pathogen_fasta")"
    pathogen_stem="${pathogen_filename%.*}"
    pathogen_group="$(basename "$(dirname "$pathogen_fasta")")"
    host_dir="$(dirname "$host_fasta")"

    raw_dir="$host_dir/$pathogen_group"
    results_dir="$raw_dir/results/$analysis_name"

    require_command python3
    require_command AltaiR

    require_script "$script_dir/run_raw.sh"
    require_script "$script_dir/GCcontent.py"
    require_script "$script_dir/Plot_RAWs.py"

    print_separator
    echo "Host-pathogen analysis workflow"
    print_separator
    echo "Pathogen FASTA     : $pathogen_fasta"
    echo "Host FASTA         : $host_fasta"
    echo "Analysis name      : $analysis_name"
    echo "k-mer range        : $kmin-$kmax"
    echo "RAW directory      : $raw_dir"
    echo "Results directory  : $results_dir"
    print_separator
    echo

    echo "[1/2] Running AltaiR RAW analysis..."

    bash "$script_dir/run_raw.sh" \
        "$host_fasta" \
        "$pathogen_fasta" \
        "$kmin" \
        "$kmax"

    echo "[1/2] Completed."
    echo

    shopt -s nullglob
    raw_files=(
        "$raw_dir/${pathogen_filename}"-k*.eg
        "$raw_dir/${pathogen_stem}"-k*.eg
    )
    shopt -u nullglob

    if [[ "${#raw_files[@]}" -eq 0 ]]; then
        fail "no AltaiR RAW files found for $pathogen_filename in: $raw_dir"
    fi

    mkdir -p "$results_dir"

    echo "[2/2] Calculating GC content and generating plots..."

    python3 "$script_dir/GCcontent.py" \
        -i "${raw_files[@]}" \
        -p "$results_dir"

    echo "[2/2] Completed."
    echo

    print_separator
    echo "Analysis workflow completed successfully"
    print_separator
    echo "Pathogen FASTA : $pathogen_fasta"
    echo "RAW files      : $raw_dir"
    echo "Plots and CSVs : $results_dir"
    print_separator
}

if [[ "$#" -eq 0 ]]; then
    usage
    exit 1
fi

case "$1" in
    -h|--help|help)
        usage
        exit 0
        ;;

    prepare)
        if [[ "$#" -lt 4 || "$#" -gt 5 ]]; then
            usage
            exit 1
        fi

        pathogen_dir_input="$2"
        output_prefix="$3"
        email="$4"
        api_key="${5:-}"

        validate_output_prefix "$output_prefix"
        pathogen_dir="$(resolve_directory "$pathogen_dir_input")"

        prepare_workflow \
            "$pathogen_dir" \
            "$output_prefix" \
            "$email" \
            "$api_key"
        ;;

    complete)
        if [[ "$#" -lt 5 || "$#" -gt 8 ]]; then
            usage
            exit 1
        fi

        pathogen_dir_input="$2"
        host_fasta_input="$3"
        analysis_name="$4"
        email="$5"
        kmin="${6:-10}"
        kmax="${7:-13}"
        api_key="${8:-}"

        validate_output_prefix "$analysis_name"
        validate_kmer_range "$kmin" "$kmax"

        pathogen_dir="$(resolve_directory "$pathogen_dir_input")"
        host_fasta="$(resolve_file "$host_fasta_input")"

        complete_workflow \
            "$pathogen_dir" \
            "$host_fasta" \
            "$analysis_name" \
            "$email" \
            "$kmin" \
            "$kmax" \
            "$api_key"
        ;;

    analyse)
        if [[ "$#" -lt 4 || "$#" -gt 6 ]]; then
            usage
            exit 1
        fi

        pathogen_input="$2"
        host_fasta_input="$3"
        analysis_name="$4"
        kmin="${5:-10}"
        kmax="${6:-13}"

        validate_output_prefix "$analysis_name"
        validate_kmer_range "$kmin" "$kmax"

        pathogen_fasta="$(resolve_pathogen_fasta "$pathogen_input" "$analysis_name")"
        host_fasta="$(resolve_file "$host_fasta_input")"

        analyse_workflow \
            "$pathogen_fasta" \
            "$host_fasta" \
            "$analysis_name" \
            "$kmin" \
            "$kmax"
        ;;

    *)
        if [[ "$#" -lt 4 || "$#" -gt 7 ]]; then
            usage
            exit 1
        fi

        pathogen_dir_input="$1"
        host_fasta_input="$2"
        output_prefix="$3"
        email="$4"
        kmin="${5:-10}"
        kmax="${6:-13}"
        api_key="${7:-}"

        validate_output_prefix "$output_prefix"
        validate_kmer_range "$kmin" "$kmax"

        pathogen_dir="$(resolve_directory "$pathogen_dir_input")"
        host_fasta="$(resolve_file "$host_fasta_input")"

        prepare_workflow \
            "$pathogen_dir" \
            "$output_prefix" \
            "$email" \
            "$api_key"

        echo

        pathogen_fasta="$pathogen_dir/${output_prefix}_joined.fasta"

        if [[ ! -s "$pathogen_fasta" ]]; then
            fail "prepared pathogen FASTA not found or empty: $pathogen_fasta"
        fi

        analyse_workflow \
            "$pathogen_fasta" \
            "$host_fasta" \
            "$output_prefix" \
            "$kmin" \
            "$kmax"
        ;;
esac
