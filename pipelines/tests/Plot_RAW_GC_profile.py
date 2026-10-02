#!/usr/bin/env python3
"""Create the two selected positional RAW and GC profile plots."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.ticker import MultipleLocator


ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "fasta_files/cassava/hosts/tme204/Begomovirus_manihotis"


def read_profile(path):
    records = []
    with path.open(newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle), 1):
            if len(row) != 5:
                raise ValueError(f"{path}:{line_number}: expected five columns")
            try:
                position = int(row[1])
                gc = float(row[3])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid position or GC") from error
            if position < 1 or not 0 <= gc <= 100:
                raise ValueError(f"{path}:{line_number}: position or GC out of range")
            records.append((position, gc))

    if not records:
        raise ValueError(f"No RAW records in {path}")

    positions = np.arange(1, max(position for position, _ in records) + 1)
    counts = np.zeros(len(positions), dtype=float)
    gc_sums = np.zeros(len(positions), dtype=float)
    target_counts = np.zeros(len(positions), dtype=float)
    for position, gc in records:
        counts[position - 1] += 1
        gc_sums[position - 1] += gc
        if 48 <= gc <= 52:
            target_counts[position - 1] += 1
    return positions, counts, gc_sums, target_counts


def window_profiles(counts, gc_sums, window):
    kernel = np.ones(window, dtype=float)
    raw_counts = np.convolve(counts, kernel, mode="same")
    window_gc_sums = np.convolve(gc_sums, kernel, mode="same")
    mean_gc = np.divide(
        window_gc_sums, raw_counts,
        out=np.full(len(counts), np.nan), where=raw_counts > 0,
    )
    return raw_counts, mean_gc


def plot_gradient_profile(
    positions, raw_counts, mean_gc, target_counts, output_path, title, window
):
    fig, raw_ax = plt.subplots(figsize=(15, 6))
    gc_ax = raw_ax.twinx()
    raw_ax.plot(positions, raw_counts, color="#2563a6", linewidth=1.7,
                label="RAWs in window")
    gc_ax.plot(positions, mean_gc, color="#b91c1c", linewidth=1.7,
               label="GC of RAWs (smoothed)")
    gc_ax.axhspan(48, 52, color="#facc15", alpha=0.25, label="48–52% GC")
    in_target = np.isfinite(mean_gc) & (mean_gc >= 48) & (mean_gc <= 52)
    gc_ax.plot(positions, np.where(in_target, mean_gc, np.nan),
               color="#ef4444", linewidth=2.8, label="GC within 48–52%")

    exact_mask = target_counts > 0
    if np.any(exact_mask):
        exact_positions = positions[exact_mask]
        exact_counts = target_counts[exact_mask]
        green_map = LinearSegmentedColormap.from_list(
            "target_count_green",
            plt.get_cmap("Greens")(np.linspace(0.25, 1.0, 256)),
        )
        green_norm = PowerNorm(gamma=0.5, vmin=0, vmax=exact_counts.max())
        raw_ax.vlines(
            exact_positions, 0, 1,
            transform=raw_ax.get_xaxis_transform(),
            colors=green_map(green_norm(exact_counts)),
            linewidth=4, alpha=0.6, zorder=0,
        )
        colorbar_ax = fig.add_axes([0.93, 0.2, 0.018, 0.6])
        colorbar = fig.colorbar(
            ScalarMappable(norm=green_norm, cmap=green_map), cax=colorbar_ax
        )
        colorbar.set_label("48–52% GC RAWs at exact position")

    raw_ax.set(
        xlabel="Position in genome", ylabel="RAWs in window",
        xlim=(positions[0], positions[-1]), ylim=(0, None),
    )
    raw_ax.xaxis.set_major_locator(MultipleLocator(500))
    gc_ax.set(ylabel="Mean RAW GC content (%)", ylim=(100, 0))
    raw_ax.grid(True, alpha=0.2)
    fig.suptitle(f"{title}: positional RAW and GC profile")
    raw_ax.set_title(f"{window}-position moving window", fontsize=10)
    lines1, labels1 = raw_ax.get_legend_handles_labels()
    lines2, labels2 = gc_ax.get_legend_handles_labels()
    raw_ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right")
    fig.subplots_adjust(left=0.07, right=0.82, bottom=0.12, top=0.88)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_exact_target_counts(
    positions, raw_counts, target_counts, output_path, title, window
):
    fig, (raw_ax, target_ax) = plt.subplots(
        2, 1, sharex=True, figsize=(15, 7),
        gridspec_kw={"height_ratios": [3, 1.5], "hspace": 0.08},
    )
    raw_ax.plot(positions, raw_counts, color="#2563a6", linewidth=1.7)
    target_ax.plot(positions, target_counts, color="#b91c1c", linewidth=1.2)
    raw_ax.set_ylabel("All RAWs in window")
    target_ax.set_ylabel("48–52% GC RAWs\nat exact position")
    target_ax.set_xlabel("Position in genome")
    for axis in (raw_ax, target_ax):
        axis.set_ylim(bottom=0)
        axis.grid(True, alpha=0.2)
    target_ax.set_xlim(positions[0], positions[-1])
    target_ax.xaxis.set_major_locator(MultipleLocator(250))
    target_ax.yaxis.set_major_locator(MultipleLocator(5))
    fig.suptitle(f"{title}: RAW counts by position")
    raw_ax.set_title(
        f"{window}-position RAW window; exact count of 48–52% GC RAWs at each position",
        fontsize=10,
    )
    fig.subplots_adjust(top=0.88, bottom=0.1, left=0.09, right=0.98)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=CASE / "merged_gc.csv")
    parser.add_argument("--output-dir", type=Path, default=CASE / "results")
    parser.add_argument("--title", default="TME204 vs Begomovirus manihotis")
    parser.add_argument("--window", type=int, default=51)
    args = parser.parse_args()
    if args.window < 1 or args.window % 2 != 1:
        parser.error("--window must be a positive odd integer")

    positions, counts, gc_sums, target_counts = read_profile(args.input)
    raw_counts, mean_gc = window_profiles(counts, gc_sums, args.window)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    gradient_path = args.output_dir / "raw_gc_profile_window51_exact_target_count_gradient.png"
    exact_count_path = args.output_dir / "raw_window51_target_count_window1.png"
    plot_gradient_profile(
        positions, raw_counts, mean_gc, target_counts,
        gradient_path, args.title, args.window,
    )
    plot_exact_target_counts(
        positions, raw_counts, target_counts,
        exact_count_path, args.title, args.window,
    )
    print(f"Saved {gradient_path}")
    print(f"Saved {exact_count_path}")
    print(f"{int(counts.sum())} RAWs; {int(target_counts.sum())} with 48–52% GC")


if __name__ == "__main__":
    main()
