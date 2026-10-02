#!/usr/bin/env python3
"""Filter headerless merged_gc.csv using the selected AltaiR's host models.

No Python reimplementation of host membership is used. For every needed length,
AltaiR checks a temporary FASTA of distinct candidates/prefixes/suffixes. Its DNA
target loop emits the preceding k bases after reading a following base, so each
query is written as word + 'A'. Only the word at position 1 is inspected.

Minimality here is relative to AltaiR's models, including their boundary and
ambiguous-base behaviour. This does not correct bugs in those models, establish
completeness of the input RAWs, or remove artificial joined-sequence padding.
Use the SAME AltaiR executable, host and --ignore-ir setting as the original run.
"""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict


def read_candidates(path):
    """Validate the GCcontent.py schema without changing occurrence coordinates."""
    rows = []
    with path.open(newline="") as handle:
        for number, row in enumerate(csv.reader(handle), 1):
            try:
                if len(row) != 5:
                    raise ValueError("expected five columns and no header")
                genome, position, word, gc, length = row
                k = int(length)
                if int(genome) < 1 or int(position) < 1:
                    raise ValueError("genome index and position must be positive integers")
                if not 2 <= k <= 20 or len(word) != k:
                    raise ValueError("word length must match k, with 2 <= k <= 20")
                if set(word) - set("ACGT"):
                    raise ValueError("only uppercase A/C/G/T DNA candidates are supported")
                if not math.isfinite(float(gc)) or not 0 <= float(gc) <= 100:
                    raise ValueError("GC percentage must be finite and between 0 and 100")
            except ValueError as error:
                raise ValueError(f"{path}:{number}: {error}") from error
            rows.append(row)
    return rows


def absent_queries(words, host, executable, ignore_ir, workdir):
    """Query one length with AltaiR; return exactly the words it reports absent."""
    words = sorted(words)
    k = len(words[0])
    query = workdir / f"queries_k{k}.fasta"
    with query.open("w") as handle:
        for index, word in enumerate(words, 1):
            # Long headers avoid the repository's progress counter dividing by
            # zero on query files shorter than 100 bytes, even for one query.
            handle.write(f">query_{index}_{'x' * 100}\n{word}A\n")
    command = [executable, "raw", "-f", "-min", str(k), "-max", str(k)]
    if ignore_ir:
        command.append("-i")
    command.extend([str(host), str(query)])
    log = workdir / f"altair_k{k}.log"
    with log.open("w") as handle:
        result = subprocess.run(command, cwd=workdir, stdout=handle,
                                stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        with log.open("rb") as handle:
            handle.seek(max(0, log.stat().st_size - 4000))
            tail = handle.read().decode(errors="replace")
        raise RuntimeError(f"AltaiR failed for k={k} (exit {result.returncode}).\n{tail}")
    output = Path(f"{query}-k{k}.eg")
    if not output.is_file():
        raise RuntimeError(f"AltaiR did not create its expected k={k} output")
    absent = set()
    with output.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise RuntimeError(f"Unexpected AltaiR output: {line!r}")
            index, position, word = fields
            index, position = int(index), int(position)
            if not 1 <= index <= len(words):
                raise RuntimeError("AltaiR returned an invalid query record index")
            if position != 1:
                raise RuntimeError("Unexpected query coordinates; this AltaiR version "
                                   "does not match the supported query protocol")
            if word != words[index - 1]:
                raise RuntimeError("AltaiR returned a different word than the query; "
                                   "refusing to infer membership")
            absent.add(word)
    return absent


def classify(word, absent):
    if word not in absent:
        return "rejected_present_in_host"
    prefix_absent = word[:-1] in absent
    suffix_absent = word[1:] in absent
    if prefix_absent and suffix_absent:
        return "rejected_both_absent"
    if prefix_absent:
        return "rejected_prefix_absent"
    if suffix_absent:
        return "rejected_suffix_absent"
    return "retained"


def atomic_write(path, writer):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", newline="", dir=path.parent,
                                         delete=False) as handle:
            temporary = Path(handle.name)
            writer(handle)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def filter_candidates(input_path, output_path, host, executable, ignore_ir=False,
                      report_path=None, temp_dir=None):
    input_path, output_path, host = (
        Path(path).resolve() for path in (input_path, output_path, host)
    )
    report_path = (Path(report_path).resolve() if report_path else
                   output_path.with_suffix(".summary.json"))
    if output_path in (input_path, host) or report_path in (input_path, host, output_path):
        raise ValueError("Input, host, output and report paths must be distinct")
    if not host.is_file() or host.stat().st_size == 0:
        raise ValueError(f"Host FASTA missing or empty: {host}")
    rows = read_candidates(input_path)
    resolved_executable = shutil.which(str(executable))
    if resolved_executable is None:
        raise ValueError(f"AltaiR executable not found: {executable}. "
                         "Use --altair with the executable from the original RAW run.")
    resolved_executable = str(Path(resolved_executable).resolve())
    if Path(resolved_executable) in (output_path, report_path):
        raise ValueError("Output/report must not overwrite the AltaiR executable")
    queries = defaultdict(set)
    candidates = {row[2] for row in rows}
    for word in candidates:
        queries[len(word)].add(word)
        queries[len(word) - 1].update((word[:-1], word[1:]))
    absent = set()
    # The selected binary's defaults are preserved; one k per process bounds
    # peak model memory and scans the host once per required length.
    with tempfile.TemporaryDirectory(prefix="altair-mraw-", dir=temp_dir) as directory:
        for k, words in sorted(queries.items()):
            candidate_count = len(words & candidates)
            factor_only_count = len(words - candidates)
            print(f"Checking {len(words):,} distinct membership queries at k={k} "
                  f"with AltaiR ({candidate_count:,} RAW candidate words; "
                  f"{factor_only_count:,} additional prefix/suffix words)...",
                  flush=True)
            absent.update(absent_queries(words, host, resolved_executable,
                                         ignore_ir, Path(directory)))
    decisions = {word: classify(word, absent) for word in candidates}
    retained = [row for row in rows if decisions[row[2]] == "retained"]
    by_k = {}
    for k in sorted({len(word) for word in candidates}):
        by_k[str(k)] = {
            "occurrences": dict(Counter(decisions[row[2]] for row in rows
                                        if len(row[2]) == k)),
            "unique_words": dict(Counter(decisions[word] for word in candidates
                                         if len(word) == k)),
        }
    report = {
        "input": str(input_path), "host": str(host), "output": str(output_path),
        "altair": resolved_executable, "ignore_ir": ignore_ir,
        "membership": "Selected AltaiR DNA host models; no independent matching",
        "input_occurrences": len(rows), "retained_occurrences": len(retained),
        "input_unique_words": len(candidates),
        "retained_unique_words": sum(value == "retained" for value in decisions.values()),
        "by_k": by_k,
    }
    atomic_write(output_path, lambda handle: csv.writer(handle).writerows(retained))
    atomic_write(report_path, lambda handle: json.dump(report, handle, indent=2))
    print(f"Retained {len(retained):,}/{len(rows):,} occurrences: {output_path}")
    print(f"Summary: {report_path}")
    if not retained:
        print("No mRAWs retained; skip plotting this empty CSV.")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path,
                        help="Headerless merged_gc.csv from GCcontent.py")
    parser.add_argument("--output", required=True, type=Path,
                        help="Separate headerless CSV; preserves occurrence order and columns")
    parser.add_argument("--altair", default="AltaiR",
                        help="Same executable as the RAW run (default: AltaiR on PATH)")
    parser.add_argument("--ignore-ir", action="store_true",
                        help="Pass -i to AltaiR; use ONLY if the original run used -i")
    parser.add_argument("--report", type=Path,
                        help="Summary JSON (default: output with .summary.json suffix)")
    parser.add_argument("--temp-dir", type=Path,
                        help="Parent directory for temporary AltaiR queries and outputs")
    args = parser.parse_args()
    try:
        filter_candidates(args.input, args.output, args.host, args.altair,
                          args.ignore_ir, args.report, args.temp_dir)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
