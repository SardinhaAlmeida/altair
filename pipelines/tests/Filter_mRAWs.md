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

## Integration with the pipeline

`pipeline.sh` already runs `GCcontent.py`, `Filter_mRAWs.py`, and `Plot_Raws.py`
in sequence. Run it from the repository root with the host first:

```bash
export PATH="$PWD/src:$PATH"
bash pipelines/tests/pipeline.sh analyse HOST_FASTA PATHOGEN_FASTA NAME 11 13
```

Results are stored in `HOST_DIR/PATHOGEN_GROUP/results/NAME/`. The pipeline
keeps `merged_gc.csv` and `merged_mraw_gc.csv`, plots RAWs in `raws/`, and plots
retained mRAWs in `mraws/`. GC-filtered exports use `raws_gc_` and `mraws_gc_`
prefixes respectively. Empty mRAW outputs skip plotting. AltaiR's CG plot still
describes all RAWs. The pipeline records provenance and counts in
`run_summary.json`; the filter also writes `merged_mraw_gc.summary.json`.

`run_raw.sh` does not pass `-i` to AltaiR, so the filter uses the same default
strand handling. If a future change adds `-i`, also add `--ignore-ir` to the
filter command and update the run summary's recorded parameters.

`GCcontent.py` imports `Plot_Raws`, so the plotting dependencies must be
installed even when only calculating GC. See the [workflow documentation](../../README.md#dissertation-workflow)
for dependencies and help commands.
