# Filtering RAW occurrences with AltaiR's membership models

`Filter_mRAWs.py` uses only the Python standard library and an AltaiR executable.
It accepts the **unfiltered, headerless** `merged_gc.csv` from `GCcontent.py`:

```text
genome_idx,position,raw,gc_pct,kmer_len
```

That line describes the columns; do not add it as a header. DNA candidates must
contain uppercase A/C/G/T, with lengths 2–20. All occurrences of retained words
are preserved, including their ordering, original coordinates and GC values.

## Run on the existing TME204 CSV

From the repository root, in the same environment as your original RAW analysis:

```bash
DATA="fasta_files/cassava/hosts/tme204/Begomovirus_manihotisafricaense"
python3 pipelines/tests/Filter_mRAWs.py \
  --host fasta_files/cassava/hosts/tme204/TME204.fa \
  --input "$DATA/merged_gc.csv" \
  --output "$DATA/merged_mraw_gc.csv"
```

By default, `AltaiR` is resolved from PATH. If necessary, use
`--altair /absolute/path/to/AltaiR` to select **the executable used for the original
RAW run**, not an arbitrary other build. Add `--ignore-ir` only when the original
run passed `-i` to AltaiR. The current `run_raw.sh` does not pass `-i`, so omit it.

Outputs are `merged_mraw_gc.csv` and `merged_mraw_gc.summary.json`. The JSON records
the executable, host path, strand flag, occurrence counts and unique-word counts,
including rejection reasons per k. Rerunning replaces these output files but does
not change the input. AltaiR's probe FASTAs, `.eg` outputs and logs are temporary.
Use `--temp-dir` to choose a different parent directory for temporary files.

## How membership is checked

For each distinct candidate `w`, retain it when AltaiR reports:

```text
w absent AND w[:-1] present AND w[1:] present
```

The script queries complete words as well as both factors, so a candidate that
AltaiR now considers present is rejected and reported. These rejections can point
to mismatched input data, host, executable or flags; inspect the summary.

For each required length, the script runs AltaiR once on a temporary FASTA of all
distinct queries of that length. The host is scanned once per required length;
model memory is bounded by running lengths sequentially. A 10–13 candidate set
can require host models for lengths 9–13. There is no need to generate a full
pathogen RAW set for k=9. This does add host-model construction time.

The current AltaiR target loop reads a following base before emitting a k-mer.
Each probe record therefore contains `word + A`, giving exactly one query at
position 1. Returned indices, words and positions are validated; unexpected
formats or failed subprocesses abort before replacing the CSV.

This is **minimality relative to the selected AltaiR models**, not an independent
correction of their matching semantics. Ambiguous-base handling, strand handling
and sequence-boundary behaviour are inherited, including implementation defects.
The filter does not establish completeness of the original RAW set. Artificial
padding and padding-boundary overlaps in joined pathogen sequences are not
removed. Those require a separate coordinate-aware filter.

## Changes to existing scripts (not applied automatically)

### `pipeline2.sh`

In `analyse_workflow`, add these checks alongside the existing `require_script`
calls:

```bash
require_script "$script_dir/Filter_mRAWs.py"
require_script "$script_dir/Plot_RAWs2.py"
```

Replace the existing `python3 "$script_dir/GCcontent.py" ... -p "$results_dir"`
command with:

```bash
# Generate all RAW GC data without plotting it.
python3 "$script_dir/GCcontent.py" -i "${raw_files[@]}"

# Reuse the AltaiR found in this pipeline's environment.
python3 "$script_dir/Filter_mRAWs.py" \
    --altair "$(command -v AltaiR)" \
    --host "$host_fasta" \
    --input "$raw_dir/merged_gc.csv" \
    --output "$raw_dir/merged_mraw_gc.csv"

# Plot separately from the existing RAW figures.
if [[ -s "$raw_dir/merged_mraw_gc.csv" ]]; then
    MPLBACKEND=Agg python3 "$script_dir/Plot_RAWs2.py" \
        "$raw_dir/merged_mraw_gc.csv" \
        "$results_dir/mraws" \
        "$pathogen_group $analysis_name — mRAWs only"
else
    echo "No mRAWs retained; no mRAW plots generated."
fi
```

If a future change adds `-i` to `run_raw.sh`, add `--ignore-ir` to this filter
command as well. The existing pipeline's `set -e` stops execution if filtering
fails. Its progress labels can be updated to describe GC calculation, filtering
and plotting separately.

The existing `run_raw.sh` removes files directly in `raw_dir` on a new run, so it
also removes these filtered CSV/report files before they are regenerated. Its
AltaiR CG plot still describes **all RAWs**, not filtered mRAWs.

### `GCcontent.py` and `run_raw.sh`

No changes are necessary. GC content is a property of each word, so it can be
calculated before minimality filtering. `GCcontent.py` already supports omitting
`-p`. Its current import of `Plot_RAWs2` still requires the plotting dependencies
even when only calculating GC.

### `Plot_RAWs2.py`

No functional change is required to read the new CSV; its schema is unchanged.
The dataset title above identifies the mRAW-only input. Generic plot titles,
legends and filenames still say RAWs. For publication figures, add a configurable
word label and use `mRAWs` consistently, rather than hardcoding it for every run.

Use a fresh output directory for a new comparison. The plotting script does not
remove old per-k or GC-filtered figures when a later run has no corresponding
records. The empty-file guard above likewise does not remove old plots.

## Validation

```bash
ALTAIR_TEST_BINARY=/absolute/path/to/AltaiR \
  python3 -m unittest discover -s pipelines/tests -p test_filter_mraws.py -v
```

Tests cover both-factor minimality, input validation, independent small-DNA
membership, duplicate occurrence preservation, direct AltaiR agreement with both
strand settings and ambiguous/multiple FASTA records, empty results, and failure
without replacing existing output. Without `ALTAIR_TEST_BINARY`, only the pure
Python tests run; integration tests are explicitly skipped.
