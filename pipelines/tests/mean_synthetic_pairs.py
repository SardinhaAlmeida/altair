#!/usr/bin/env python3

"""Average the three pathogen CG distributions for each synthetic pair."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path


REPLICATES = ("p1", "p2", "p3")


def read_cg_data(path: Path) -> dict[int, tuple[float, float]]:
    """Read k-mer, AT percentage, and CG percentage from CG-data.eg."""
    values: dict[int, tuple[float, float]] = {}

    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue

            fields = line.split()
            if len(fields) != 3:
                raise ValueError(
                    f"{path}:{line_number}: expected exactly three columns"
                )

            kmer = int(fields[0])
            values[kmer] = (float(fields[1]), float(fields[2]))

    return values


def calculate_means(
    replicate_data: list[dict[int, tuple[float, float]]],
) -> list[tuple[int, float, float]]:
    """Average only k-mers that have valid data in all three replicates."""
    common_kmers = set(replicate_data[0])
    for data in replicate_data[1:]:
        common_kmers.intersection_update(data)

    means: list[tuple[int, float, float]] = []
    for kmer in sorted(common_kmers):
        mean_at = sum(data[kmer][0] for data in replicate_data) / len(
            replicate_data
        )
        # Each source row is rounded, so derive CG from AT to preserve 100%.
        mean_cg = 100.0 - mean_at
        means.append((kmer, mean_at, mean_cg))

    return means


def write_mean_data(
    output_path: Path,
    means: list[tuple[int, float, float]],
) -> None:
    """Write the headerless three-column format expected by CGplot.sh."""
    with output_path.open("w", encoding="utf-8") as handle:
        for kmer, mean_at, mean_cg in means:
            handle.write(f"{kmer}\t{mean_at:.3f}\t{mean_cg:.3f}\n")


def create_plot(
    cgplot_script: Path,
    mean_data: Path,
    output_path: Path,
) -> None:
    """
    Run the existing CGplot.sh without changing it.

    CGplot.sh expects CG-data.eg and creates CGplot.pdf in its working
    directory, so use a temporary directory and then copy the resulting PDF
    to the requested pair-level filename.
    """
    with tempfile.TemporaryDirectory(prefix="altair-cg-mean-") as temp_name:
        temp_dir = Path(temp_name)
        shutil.copyfile(mean_data, temp_dir / "CG-data.eg")

        subprocess.run(
            ["sh", str(cgplot_script.resolve())],
            cwd=temp_dir,
            check=True,
        )

        generated_plot = temp_dir / "CGplot.pdf"
        if not generated_plot.is_file() or generated_plot.stat().st_size == 0:
            raise RuntimeError(
                f"{cgplot_script} did not create a non-empty CGplot.pdf"
            )

        shutil.copyfile(generated_plot, output_path)


def process_pair(pair_dir: Path) -> None:
    data_files = [
        pair_dir / replicate / "CG-data.eg" for replicate in REPLICATES
    ]
    missing = [path for path in data_files if not path.is_file()]
    if missing:
        missing_text = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(
            f"{pair_dir.name}: missing replicate data: {missing_text}"
        )

    replicate_data = [read_cg_data(path) for path in data_files]
    means = calculate_means(replicate_data)
    if not means:
        raise ValueError(
            f"{pair_dir.name}: no k-mer has valid data in all three replicates"
        )

    plot_scripts = [
        pair_dir / replicate / "CGplot.sh" for replicate in REPLICATES
    ]
    cgplot_script = next(
        (path for path in plot_scripts if path.is_file()),
        None,
    )
    if cgplot_script is None:
        raise FileNotFoundError(
            f"{pair_dir.name}: no replicate CGplot.sh was found"
        )

    mean_data_path = pair_dir / "CG_mean_plot.eg"
    mean_plot_path = pair_dir / "CG_mean_plot.pdf"

    write_mean_data(mean_data_path, means)
    create_plot(cgplot_script, mean_data_path, mean_plot_path)

    used_kmers = ", ".join(str(kmer) for kmer, _, _ in means)
    print(f"{pair_dir.name}: k={used_kmers}")
    print(f"  {mean_data_path}")
    print(f"  {mean_plot_path}")


def main() -> None:
    default_root = Path(__file__).resolve().parent / "synthetic_pairs"

    parser = argparse.ArgumentParser(
        description=(
            "Average p1, p2, and p3 CG-data.eg values and create one "
            "pair-level CG plot with the existing CGplot.sh."
        )
    )
    parser.add_argument(
        "synthetic_pairs",
        nargs="?",
        type=Path,
        default=default_root,
        help=(
            "synthetic_pairs directory "
            f"(default: {default_root})"
        ),
    )
    args = parser.parse_args()

    root = args.synthetic_pairs.resolve()
    if not root.is_dir():
        parser.error(f"directory does not exist: {root}")

    pair_dirs = sorted(
        path
        for path in root.iterdir()
        if path.is_dir() and path.name.startswith("H")
    )
    if not pair_dirs:
        parser.error(f"no pair directories found under: {root}")

    if shutil.which("gnuplot") is None:
        parser.error("gnuplot is required to create CG_mean_plot.pdf")

    for pair_dir in pair_dirs:
        process_pair(pair_dir)


if __name__ == "__main__":
    main()
