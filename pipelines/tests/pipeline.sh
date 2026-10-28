#!/usr/bin/env bash

set -euo pipefail

pipeline_started_at=$SECONDS
pipeline_started_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

usage() {
    cat >&2 <<USAGE
Usage:
  $0 prepare <pathogen_directory> <output_prefix> <email> [api_key]

  $0 complete <host_fasta> <pathogen_directory> <analysis_name> <email> \
[kmin] [kmax] [api_key]

  $0 analyse <host_fasta> <pathogen_fasta_or_directory> <analysis_name> \
[kmin] [kmax]

  $0 <host_fasta> <pathogen_directory> <output_prefix> <email> \
[kmin] [kmax] [api_key]

Workflows:
  Run from the altair repository folder with AltaiR available on PATH:
      export PATH="\$PWD/src:\$PATH"
      bash pipelines/tests/pipeline.sh --help

  prepare
      (This command was made for Bipartite Begomoviruses to handle DNA-A and DNA-B components)
      Classify the metadata, download nucleotide sequences, label records
      without a DNA-A/DNA-B label by sequence similarity, and write unpadded
      DNA components plus a .layout.tsv for plotting. Optional
      manual_components.csv and manual_pairs.csv files in the pathogen
      directory are applied, as are the quality filters listed in an optional
      segments_options.txt (see GetSegments.py).

  complete
      Download every accession from sequences.csv as an independent complete
      genome, then run the RAW analysis. This workflow does not classify or
      join DNA-A/DNA-B components.

  analyse
      Analyse an existing pathogen FASTA. After the host FASTA, supply:
        - a FASTA file; or
        - a directory containing the FASTA.

      If a directory is supplied, the script first looks for files named from
      <analysis_name>, preferring <analysis_name>_segments.fasta over
      legacy <analysis_name>_joined.fasta and
      <analysis_name>.fasta. If neither exists, the directory must contain
      exactly one non-empty .fasta, .fa, or .fna file.

      All analysis outputs are stored under the host directory in
      <pathogen_group>/results/<analysis_name>/.
      Use a different analysis name to preserve each run.
      A run_summary.json records inputs, settings, timing, and output counts.
      Existing output files can be overwritten when an analysis name is reused.

  no workflow argument
      Run the complete preparation and analysis pipeline for a segmented
      pathogen dataset.

Examples:

  $0 complete \
     fasta_files/human/human.fna \
     fasta_files/virus_complete_genomes \
     VIRUS \
     user@example.com \
     10 \
     13

  $0 prepare \
     fasta_files/cassava/pathogens/CassavaMosaicVirus/Begomovirus_manihotis \
     ACMV \
     user@example.com

  $0 analyse \
     fasta_files/cassava/hosts/tmeb117/TMEB117.fa \
     fasta_files/cassava/pathogens/CassavaMosaicVirus/Begomovirus_manihotis/ACMV_segments.fasta \
     ACMV \
     11 \
     13

  $0 \
     fasta_files/cassava/hosts/tmeb117/TMEB117.fa \
     fasta_files/cassava/pathogens/CassavaMosaicVirus/Begomovirus_manihotis \
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

finish_run_summary() {
    local status=$?
    trap - EXIT
    if ! python3 "$script_dir/run_summary.py" finish "$run_summary_path" \
        --exit-code "$status" --elapsed-seconds "$((SECONDS - pipeline_started_at))"; then
        echo "ERROR: could not finalize run summary: $run_summary_path" >&2
        if [[ "$status" -eq 0 ]]; then
            status=1
        fi
    fi
    exit "$status"
}

start_run_summary() {
    local workflow="$1"
    local pathogen_fasta="$2"
    local host_fasta="${3:-}"
    local dataset_dir
    local -a summary_args=()

    dataset_dir="$(dirname "$pathogen_fasta")"
    if [[ "$workflow" == prepare ]]; then
        run_summary_path="$dataset_dir/paired_output/${output_prefix}_run_summary.json"
        summary_args=(--name "$output_prefix")
    else
        run_summary_path="$(dirname "$host_fasta")/$(basename "$dataset_dir")/results/$analysis_name/run_summary.json"
        summary_args=(--name "$analysis_name" --host "$host_fasta" --kmin "$kmin" --kmax "$kmax")
    fi
    require_command python3
    require_script "$script_dir/run_summary.py"
    python3 "$script_dir/run_summary.py" start "$run_summary_path" \
        --workflow "$workflow" --pathogen "$pathogen_fasta" \
        --started-at "$pipeline_started_utc" "${summary_args[@]}"
    trap finish_run_summary EXIT
    echo "Run summary: $run_summary_path"
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
        "$pathogen_dir/${analysis_name}_segments.fasta" \
        "$pathogen_dir/${analysis_name}_segments.fa" \
        "$pathogen_dir/${analysis_name}_segments.fna" \
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
    local components_fasta="$pathogen_dir/${output_prefix}_segments.fasta"
    local layout_metadata="$pathogen_dir/${output_prefix}_segments.layout.tsv"
    local required_csv

    if [[ -z "$email" ]]; then
        fail "email cannot be empty"
    fi

    if [[ ! -f "$metadata_csv" ]]; then
        fail "metadata CSV not found: $metadata_csv"
    fi

    require_command python3
    require_command curl

    require_script "$script_dir/GetSegments.py"
    require_script "$script_dir/get_fasta.sh"
    require_script "$script_dir/JoinSegments.py"

    print_separator
    echo "Pathogen data-preparation workflow"
    print_separator
    echo "Pathogen directory : $pathogen_dir"
    echo "Output prefix      : $output_prefix"
    echo "Metadata CSV       : $metadata_csv"
    echo "Downloaded FASTA   : $downloaded_fasta"
    echo "Component FASTA    : $components_fasta"
    echo "Display layout     : $layout_metadata"
    print_separator
    echo

    echo "[1/4] Classifying viral genome components from metadata..."

    python3 "$script_dir/GetSegments.py" \
        "$metadata_csv"

    echo "[1/4] Completed."
    echo

    echo "[2/4] Downloading and validating nucleotide sequences..."

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

    echo "[2/4] Completed."
    echo

    if [[ ! -s "$downloaded_fasta" ]]; then
        fail "downloaded FASTA not found or empty: $downloaded_fasta"
    fi

    echo "[3/4] Labelling unknown components by similarity and pairing again..."

    python3 "$script_dir/GetSegments.py" \
        "$metadata_csv" \
        --fasta "$downloaded_fasta" 

    echo "[3/4] Completed."
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

    echo "[4/4] Organising paired and unpaired DNA components..."

    python3 "$script_dir/JoinSegments.py" \
        "$paired_dir/pairs.csv" \
        "$paired_dir/unpaired_A.csv" \
        "$paired_dir/unpaired_B.csv" \
        "$downloaded_fasta" \
        "$components_fasta" \
        --metadata "$layout_metadata"

    if [[ ! -s "$components_fasta" || ! -s "$layout_metadata" ]]; then
        fail "component FASTA or display layout was not created: $components_fasta"
    fi

    echo "[4/4] Completed."
    echo

    print_separator
    echo "Preparation workflow completed successfully"
    print_separator
    echo "Downloaded FASTA : $downloaded_fasta"
    echo "Component FASTA  : $components_fasta"
    echo "Display layout   : $layout_metadata"
    print_separator
}

complete_workflow() {
    local host_fasta="$1"
    local pathogen_dir="$2"
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
        "$host_fasta" \
        "$downloaded_fasta" \
        "$analysis_name" \
        "$kmin" \
        "$kmax"
}

analyse_workflow() {
    local host_fasta="$1"
    local pathogen_fasta="$2"
    local analysis_name="$3"
    local kmin="$4"
    local kmax="$5"

    local pathogen_filename
    local pathogen_stem
    local pathogen_group
    local host_dir
    local raw_dir
    local results_dir
    local k
    local -a raw_files
    local layout_metadata="${pathogen_fasta%.*}.layout.tsv"
    local -a metadata_args=()

    if [[ -f "$layout_metadata" ]]; then
        [[ -s "$layout_metadata" ]] || fail "display layout is empty: $layout_metadata"
        metadata_args=(--metadata "$layout_metadata")
    elif [[ "$(basename "$pathogen_fasta")" == *_segments.* ]] || \
        grep -qE '^>[^[:space:]]+ \| genome_idx=[0-9]+ \| component=[AB] \| status=' "$pathogen_fasta"; then
        fail "unpadded components require their display layout: $layout_metadata"
    fi

    pathogen_filename="$(basename "$pathogen_fasta")"
    pathogen_stem="${pathogen_filename%.*}"
    pathogen_group="$(basename "$(dirname "$pathogen_fasta")")"
    host_dir="$(dirname "$host_fasta")"

    results_dir="$host_dir/$pathogen_group/results/$analysis_name"
    raw_dir="$results_dir"

    require_command python3
    require_command AltaiR

    require_script "$script_dir/run_raw.sh"
    require_script "$script_dir/GCcontent.py"
    require_script "$script_dir/Filter_mRAWs.py"
    require_script "$script_dir/Plot_Raws.py"

    print_separator
    echo "Host-pathogen analysis workflow"
    print_separator
    echo "Host FASTA         : $host_fasta"
    echo "Pathogen FASTA     : $pathogen_fasta"
    echo "Analysis name      : $analysis_name"
    echo "k-mer range        : $kmin-$kmax"
    echo "RAW directory      : $raw_dir"
    echo "Results directory  : $results_dir"
    if [[ "${#metadata_args[@]}" -gt 0 ]]; then
        echo "Display layout     : $layout_metadata"
    fi
    print_separator
    echo

    echo "[1/2] Running AltaiR RAW analysis..."

    bash "$script_dir/run_raw.sh" \
        "$host_fasta" \
        "$pathogen_fasta" \
        "$kmin" \
        "$kmax" \
        "$results_dir" \
        "$analysis_name"

    echo "[1/2] Completed."
    echo

    shopt -s nullglob
    raw_files=()
    for ((k = kmin; k <= kmax; k++)); do
        if [[ -f "$raw_dir/${pathogen_filename}-k${k}.eg" ]]; then
            raw_files+=("$raw_dir/${pathogen_filename}-k${k}.eg")
        elif [[ -f "$raw_dir/${pathogen_stem}-k${k}.eg" ]]; then
            raw_files+=("$raw_dir/${pathogen_stem}-k${k}.eg")
        fi
    done
    shopt -u nullglob

    if [[ "${#raw_files[@]}" -eq 0 ]]; then
        fail "no AltaiR RAW files found for $pathogen_filename in: $raw_dir"
    fi

    mkdir -p "$results_dir"
    if [[ "${#metadata_args[@]}" -gt 0 ]]; then
        if [[ "$layout_metadata" != "$results_dir/${pathogen_stem}.layout.tsv" ]]; then
            cp "$layout_metadata" "$results_dir/${pathogen_stem}.layout.tsv"
        fi
        metadata_args=(--metadata "$results_dir/${pathogen_stem}.layout.tsv")
    fi

    echo "[2/2] Calculating GC content and generating plots..."

    python3 "$script_dir/GCcontent.py" \
        -i "${raw_files[@]}"

    local host_name
    host_name="$(basename "$host_fasta")"
    host_name="${host_name%.*}"

    python3 "$script_dir/Filter_mRAWs.py" \
        --altair "$(command -v AltaiR)" \
        --host "$host_fasta" \
        --input "$raw_dir/merged_gc.csv" \
        --output "$raw_dir/merged_mraw_gc.csv"

    python3 "$script_dir/Plot_Raws.py" \
        "$raw_dir/merged_gc.csv" "$results_dir/raws" \
        --pathogen-name "$pathogen_group" --host-name "$host_name" \
        --word-label RAWs "${metadata_args[@]}"

    if [[ -s "$raw_dir/merged_mraw_gc.csv" ]]; then
        MPLBACKEND=Agg python3 "$script_dir/Plot_Raws.py" \
            "$raw_dir/merged_mraw_gc.csv" \
            "$results_dir/mraws" \
            --pathogen-name "$pathogen_group" --host-name "$host_name" \
            --word-label mRAWs "${metadata_args[@]}"
    else
        echo "No mRAWs retained; skipping plots."
    fi

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
        if [[ "${2:-}" == -h || "${2:-}" == --help ]]; then
            usage
            exit 0
        fi
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

        start_run_summary prepare "$pathogen_dir/${output_prefix}_segments.fasta"

        prepare_workflow \
            "$pathogen_dir" \
            "$output_prefix" \
            "$email" \
            "$api_key"
        ;;

    complete)
        if [[ "${2:-}" == -h || "${2:-}" == --help ]]; then
            usage
            exit 0
        fi
        if [[ "$#" -lt 5 || "$#" -gt 8 ]]; then
            usage
            exit 1
        fi

        host_fasta_input="$2"
        pathogen_dir_input="$3"
        analysis_name="$4"
        email="$5"
        kmin="${6:-10}"
        kmax="${7:-13}"
        api_key="${8:-}"

        validate_output_prefix "$analysis_name"
        validate_kmer_range "$kmin" "$kmax"

        pathogen_dir="$(resolve_directory "$pathogen_dir_input")"
        host_fasta="$(resolve_file "$host_fasta_input")"

        start_run_summary complete "$pathogen_dir/${analysis_name}.fasta" "$host_fasta"

        complete_workflow \
            "$host_fasta" \
            "$pathogen_dir" \
            "$analysis_name" \
            "$email" \
            "$kmin" \
            "$kmax" \
            "$api_key"
        ;;

    analyse)
        if [[ "${2:-}" == -h || "${2:-}" == --help ]]; then
            usage
            exit 0
        fi
        if [[ "$#" -lt 4 || "$#" -gt 6 ]]; then
            usage
            exit 1
        fi

        host_fasta_input="$2"
        pathogen_input="$3"
        analysis_name="$4"
        kmin="${5:-10}"
        kmax="${6:-13}"

        validate_output_prefix "$analysis_name"
        validate_kmer_range "$kmin" "$kmax"

        pathogen_fasta="$(resolve_pathogen_fasta "$pathogen_input" "$analysis_name")"
        host_fasta="$(resolve_file "$host_fasta_input")"

        start_run_summary analyse "$pathogen_fasta" "$host_fasta"

        analyse_workflow \
            "$host_fasta" \
            "$pathogen_fasta" \
            "$analysis_name" \
            "$kmin" \
            "$kmax"
        ;;

    *)
        if [[ "$#" -lt 4 || "$#" -gt 7 ]]; then
            usage
            exit 1
        fi

        host_fasta_input="$1"
        pathogen_dir_input="$2"
        output_prefix="$3"
        email="$4"
        kmin="${5:-10}"
        kmax="${6:-13}"
        api_key="${7:-}"

        validate_output_prefix "$output_prefix"
        validate_kmer_range "$kmin" "$kmax"

        pathogen_dir="$(resolve_directory "$pathogen_dir_input")"
        host_fasta="$(resolve_file "$host_fasta_input")"

        analysis_name="$output_prefix"
        start_run_summary segmented "$pathogen_dir/${output_prefix}_segments.fasta" "$host_fasta"

        prepare_workflow \
            "$pathogen_dir" \
            "$output_prefix" \
            "$email" \
            "$api_key"

        echo

        pathogen_fasta="$pathogen_dir/${output_prefix}_segments.fasta"

        if [[ ! -s "$pathogen_fasta" ]]; then
            fail "prepared pathogen FASTA not found or empty: $pathogen_fasta"
        fi

        analyse_workflow \
            "$host_fasta" \
            "$pathogen_fasta" \
            "$output_prefix" \
            "$kmin" \
            "$kmax"
        ;;
esac

pipeline_elapsed=$((SECONDS - pipeline_started_at))
printf "Total pipeline runtime: %02dh %02dm %02ds (%d seconds)\n" \
    "$((pipeline_elapsed / 3600))" "$((pipeline_elapsed / 60 % 60))" \
    "$((pipeline_elapsed % 60))" "$pipeline_elapsed"
