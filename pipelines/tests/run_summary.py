#!/usr/bin/env python3
"""Write the pipeline's JSON run record using only the Python standard library."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def fingerprint(path):
    path = Path(path).resolve()
    record = {"path": str(path), "exists": path.is_file()}
    if record["exists"]:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        record.update(size_bytes=path.stat().st_size, sha256=digest.hexdigest())
    return record


def write_summary(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, encoding="utf-8",
                                         delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(data, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def start(args):
    pathogen = args.pathogen.resolve()
    inputs = {"pathogen_fasta": fingerprint(pathogen)}
    if args.host:
        inputs["host_fasta"] = fingerprint(args.host)
    if args.workflow != "analyse":
        for name in ("sequences.csv", "manual_components.csv", "manual_pairs.csv",
                     "segments_options.txt"):
            inputs[name] = fingerprint(pathogen.parent / name)
    inputs["layout"] = fingerprint(pathogen.with_suffix(".layout.tsv"))
    executable = shutil.which("AltaiR")
    altair = None
    if executable and args.host:
        altair = fingerprint(executable)
        try:
            result = subprocess.run([executable, "-V"], capture_output=True,
                                    text=True, timeout=10, check=False)
            altair["version"] = (result.stdout or result.stderr).strip() if result.returncode == 0 else None
        except (OSError, subprocess.TimeoutExpired):
            altair["version"] = None
    script_dir = Path(__file__).resolve().parent
    scripts = {name: fingerprint(script_dir / name) for name in (
        "pipeline.sh", "get_fasta.sh", "GetSegments.py", "JoinSegments.py",
        "run_raw.sh", "GCcontent.py", "Filter_mRAWs.py", "Plot_Raws.py", "run_summary.py")}
    write_summary(args.summary, {
        "schema_version": 1,
        "status": "running", "workflow": args.workflow, "analysis_name": args.name,
        "started_at_utc": args.started_at, "finished_at_utc": None,
        "elapsed_seconds": None, "exit_code": None,
        "working_directory": str(Path.cwd()), "inputs": inputs,
        "parameters": {"kmin": args.kmin, "kmax": args.kmax,
                       "ignore_ir": False if args.host else None},
        "altair": altair, "python_version": sys.version,
        "scripts": scripts, "output_directory": str(args.summary.parent.resolve()),
    })


def csv_count(path):
    with path.open(newline="") as handle:
        return sum(1 for row in csv.reader(handle) if row)


def finish(args):
    data = json.loads(args.summary.read_text())
    data.update(status="completed" if args.exit_code == 0 else "failed",
                exit_code=args.exit_code, finished_at_utc=utc_now(),
                elapsed_seconds=args.elapsed_seconds)
    # Preserve the pre-run fingerprints; record generated preparation artifacts
    # separately so a stale FASTA present before preparation is not mistaken for
    # the FASTA that was analysed.
    pathogen = Path(data["inputs"]["pathogen_fasta"]["path"])
    if data["workflow"] != "analyse":
        data["prepared_artifacts"] = {
            "pathogen_fasta": fingerprint(pathogen),
            "layout": fingerprint(pathogen.with_suffix(".layout.tsv")),
        }
    output_dir = args.summary.parent
    counts = {}
    kmin, kmax = data["parameters"]["kmin"], data["parameters"]["kmax"]
    if kmin is not None:
        raw_counts = {}
        for k in range(kmin, kmax + 1):
            for name in (f"{pathogen.name}-k{k}.eg", f"{pathogen.stem}-k{k}.eg"):
                path = output_dir / name
                if path.is_file():
                    with path.open() as handle:
                        raw_counts[str(k)] = sum(1 for line in handle if line.strip())
                    break
        counts["raw_occurrences_by_k"] = raw_counts
        counts["missing_raw_k"] = [k for k in range(kmin, kmax + 1) if str(k) not in raw_counts]
        for name, label in (("merged_gc.csv", "raw_occurrences"),
                            ("merged_mraw_gc.csv", "mraw_occurrences")):
            path = output_dir / name
            if path.is_file():
                counts[label] = csv_count(path)
    else:
        for name in ("pairs.csv", "unpaired_A.csv", "unpaired_B.csv",
                     "unknown_component.csv", "excluded_records.csv"):
            path = output_dir / name
            if path.is_file():
                counts[name] = max(0, csv_count(path) - 1)
    data["output_counts"] = counts
    data["output_files"] = sorted(str(path.relative_to(output_dir))
                                  for path in output_dir.rglob("*")
                                  if path.is_file() and path != args.summary)
    data["output_inventory_note"] = (
        "Files present at completion; an existing analysis directory may also contain earlier outputs."
    )
    write_summary(args.summary, data)
    print(f"Run summary: {args.summary}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    begin = commands.add_parser("start", help="Record inputs before starting work")
    begin.add_argument("summary", type=Path)
    begin.add_argument("--workflow", choices=("prepare", "complete", "analyse", "segmented"), required=True)
    begin.add_argument("--name", required=True)
    begin.add_argument("--pathogen", type=Path, required=True)
    begin.add_argument("--host", type=Path)
    begin.add_argument("--kmin", type=int)
    begin.add_argument("--kmax", type=int)
    begin.add_argument("--started-at", required=True)
    end = commands.add_parser("finish", help="Record completion or failure and count outputs")
    end.add_argument("summary", type=Path)
    end.add_argument("--exit-code", type=int, required=True)
    end.add_argument("--elapsed-seconds", type=int, required=True)
    args = parser.parse_args()
    args.summary = args.summary.resolve()
    if args.command == "start":
        start(args)
    else:
        finish(args)


if __name__ == "__main__":
    main()
