#!/usr/bin/env python3
"""Original RAW figures plus offline interactive positional HTML plots, for RAWs or mRAWs.

Input: headerless genome_idx,position,raw,gc_pct,kmer_len CSV.
Profiles pool occurrences at the same position across all genome indices and k.
GC is the occurrence-weighted GC of words, not the reference genome's GC.
Optional --metadata maps unpadded component records to genome rows and display positions.
"""
import argparse
import json
import html
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
from matplotlib.collections import LineCollection
from matplotlib.ticker import MaxNLocator, MultipleLocator
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.text import Text


def save_plot(word_label, fig, path, **kwargs):
    """Keep the original artwork while labeling mRAW figures consistently."""
    if word_label == "mRAWs":
        for text in fig.findobj(Text):
            text.set_text(text.get_text().replace("RAW", "mRAW"))
    fig.savefig(path, **kwargs)
    plt.close(fig)

import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex

# ---- helpers ---------------------------------------------------------------
def ensure_outdir(d):
    os.makedirs(d, exist_ok=True)

def read_segment_layout(path):
    """Read the record-order mapping written beside an unpadded component FASTA."""
    layout = pd.read_csv(path, sep="\t", keep_default_na=False)
    required = {"record_idx", "genome_idx", "component", "accession", "length", "b_offset"}
    if not required.issubset(layout.columns) or layout.empty:
        raise ValueError("Layout must contain component records and columns: " + ", ".join(sorted(required)))
    for column in ("record_idx", "genome_idx", "length", "b_offset"):
        values = pd.to_numeric(layout[column], errors="raise")
        if not np.isfinite(values).all() or (values % 1 != 0).any():
            raise ValueError(f"Layout {column} must contain finite integers")
        if (values < (0 if column == "b_offset" else 1)).any():
            raise ValueError(f"Invalid layout {column}")
        layout[column] = values.astype(int)
    if layout.record_idx.tolist() != list(range(1, len(layout) + 1)):
        raise ValueError("Layout record_idx must follow FASTA record order, starting at 1")
    if not layout.component.isin(("A", "B")).all() or (layout.accession == "").any():
        raise ValueError("Layout requires A/B components and nonempty accessions")
    if layout.duplicated(["genome_idx", "component"]).any():
        raise ValueError("Layout has duplicate components for the same genome")
    if sorted(layout.genome_idx.unique()) != list(range(1, int(layout.genome_idx.max()) + 1)):
        raise ValueError("Layout genome indices must be consecutive, starting at 1")
    max_a = int(layout.loc[layout.component == "A", "length"].max()) if (layout.component == "A").any() else 0
    if layout.b_offset.nunique() != 1 or int(layout.b_offset.iloc[0]) < max_a + 500:
        raise ValueError("Layout requires one common DNA-B offset with at least 500 empty positions")
    return layout


def apply_segment_layout(df, metadata_path):
    """Keep source coordinates as columns; transform a copy for display only."""
    layout = read_segment_layout(metadata_path)
    records = layout.set_index("record_idx")
    result = df.copy()
    result["record_idx"] = result.genome_idx.astype(int)
    if not result.record_idx.isin(records.index).all():
        raise ValueError("RAW record index is missing from layout metadata")
    result["segment_position"] = result.position.astype(int)
    lengths = result.record_idx.map(records.length)
    if (result.segment_position + result.kmer_len - 1 > lengths).any():
        raise ValueError("RAW occurrence extends beyond the component length in metadata")
    for column in ("genome_idx", "component", "accession"):
        result[column] = result.record_idx.map(records[column])
    offsets = np.where(result.component == "B", int(layout.b_offset.iloc[0]), 0)
    result["position"] = result.segment_position + offsets
    ends = layout.length + np.where(layout.component == "B", layout.b_offset, 0)
    result.attrs["segment_layout"] = {
        "records": layout.to_dict("records"),
        "position_max": int(ends.max()),
        "genome_max": int(layout.genome_idx.max()),
        "b_offset": int(layout.b_offset.iloc[0]),
    }
    return result


def display_extents(df):
    layout = df.attrs.get("segment_layout", {})
    end = layout.get("position_max", int(df.position.max()))
    last_genome = layout.get("genome_max", int(df.genome_idx.max()))
    return max(250, int(np.ceil(end / 250)) * 250), max(5, int(np.ceil(last_genome / 5)) * 5)


def kmer_styles(df):
    """Use the original Plot_RAWs.py colours, with a fallback for other lengths."""
    kmer_lengths = sorted(df["kmer_len"].unique())
    color_map = plt.get_cmap("tab10", len(kmer_lengths))
    original_colors = {11: "red", 12: "blue", 13: "darkseagreen", 14: "#e87ba4"}
    k_colors = {
        k: original_colors.get(k, to_hex(color_map(index)))
        for index, k in enumerate(kmer_lengths)
    }
    k_sizes = {k: 3 for k in kmer_lengths}
    return kmer_lengths, k_colors, k_sizes

def topdown_axes(ax, x, y, title, layout=None):
    """Shared per-k limits and a fine grid with readable, aligned tick labels."""
    position_max = max(250, int(np.ceil(x.max() / 250)) * 250)
    genome_max = max(5, int(np.ceil(y.max() / 5)) * 5)
    if layout:
        position_max = max(250, int(np.ceil(layout["position_max"] / 250)) * 250)
        genome_max = max(5, int(np.ceil(layout["genome_max"] / 5)) * 5)
    # Use multiples of the distribution views' 250 bp / 5 genome intervals.
    # For larger datasets, retain the fine grid but space the labels farther apart.
    position_step = 250 * max(1, int(np.ceil(position_max / (250 * 28))))
    genome_step = 5 * max(1, int(np.ceil(genome_max / (5 * 32))))
    ax.set_title(title, pad=18, fontsize=18)
    ax.set_xlabel("Display position" if layout else "Position in Genome", fontsize=18, labelpad=12)
    ax.set_ylabel("Genome Index", fontsize=18, labelpad=12)
    ax.set_xlim(0, position_max)
    ax.set_ylim(-0.5, genome_max + 0.5)
    ax.xaxis.set_major_locator(MultipleLocator(position_step))
    ax.yaxis.set_major_locator(MultipleLocator(genome_step))
    ax.xaxis.set_minor_locator(MultipleLocator(50 if position_step == 250 else 250))
    ax.yaxis.set_minor_locator(MultipleLocator(1 if genome_step == 5 else 5))
    ax.tick_params(axis="both", which="major", labelsize=16, length=5, pad=6)
    ax.tick_params(axis="x", labelrotation=90)
    ax.tick_params(axis="both", which="minor", length=2)
    ax.grid(which="major", color="0.65", linewidth=0.7)
    ax.grid(which="minor", color="0.87", linewidth=0.4)
    ax.set_axisbelow(True)

def distribution_kmer_styles(df):
    """Use the older red/blue/green/grey palette for the distributions."""
    return kmer_styles(df)

def scatter_occurrences(ax, df, k_colors, ordinary_size=1):
    """Emphasise the input dataset's shortest k, even in filtered views."""
    shortest_k = min(k_colors)
    shortest = df.kmer_len.eq(shortest_k)
    for selected, size, layer in (
        (~shortest, ordinary_size, 2),
        (shortest, 3, 3),
    ):
        rows = df.loc[selected]
        if rows.empty:
            continue
        coordinates = (rows.genome_idx, rows.position, rows.kmer_len) if ax.name == "3d" else (rows.position, rows.genome_idx)
        options = {"depthshade": False} if ax.name == "3d" else {}
        ax.scatter(*coordinates, c=rows.kmer_len.map(k_colors), s=size,
                   alpha=0.9, zorder=layer, **options)


def distribution_topdown_axes(ax, x, y, title, layout=None):
    heading, _, context = title.partition(" [")
    heading = heading.replace(" (GC", "\n(GC")
    ax.set_title(heading + "\n" + context.rstrip("]"), fontsize=14, y=0.95, pad=0)
    ax.set_xlabel("Genome Index")
    ax.set_ylabel("Display position" if layout else "Position in Genome", labelpad=15)
    ax.set_xlim3d(0.5 if layout else x.min() - 0.5,
                  layout["genome_max"] + 0.5 if layout else x.max() + 0.5)
    position_max = max(250, int(np.ceil(y.max() / 250)) * 250)
    if layout:
        position_max = max(250, int(np.ceil(layout["position_max"] / 250)) * 250)
    ax.set_yticks(np.arange(0, position_max + 1, 250 * max(1, int(np.ceil(position_max / (250 * 28))))))
    ax.set_ylim3d(position_max, 0)  # azim=180 projects increasing y to the left
    if len(x): ax.set_xticks(np.arange(0, layout["genome_max"] + 1 if layout else x.max(), 5))
    ax.set_zticks([])
    ax.set_zlabel("")
    # Exclude the empty axis from tight-bbox calculations as well as drawing.
    ax.zaxis.set_visible(False)
    ax.invert_zaxis()
    ax.view_init(elev=90, azim=180)  # top-down
    ax.set_proj_type('ortho')
    ax.set_box_aspect([7.5, 8.1, 1], zoom=1.08)
    ax.yaxis.set_label_coords(-0.05, 0.5)
    ax.tick_params(axis='y', labelrotation=90)

def save_topdown_plot(df, fig, path, **kwargs):
    """Crop to projected data and visible annotations, not the empty 3D axes box."""
    from itertools import product
    from matplotlib.transforms import Bbox
    from mpl_toolkits.mplot3d import proj3d

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax = fig.axes[0]
    corners = np.array(list(product(ax.get_xlim(), ax.get_ylim(), ax.get_zlim())))
    px, py, _ = proj3d.proj_transform(*corners.T, ax.get_proj())
    pixels = ax.transData.transform(np.column_stack([px, py]))
    # Align the title and legend in one compact row above the projected plot.
    legend = ax.get_legend()
    title_box = ax.title.get_window_extent(renderer)
    legend_box = legend.get_window_extent(renderer) if legend is not None else None
    header_height = max(title_box.height, legend_box.height if legend_box else 0)
    gap = fig.dpi * 0.30
    header_top = pixels[:, 1].max() + gap + header_height
    inverse = ax.transAxes.inverted()
    ax.title.set_position(inverse.transform([
        pixels[:, 0].mean(), header_top - header_height / 2]))
    ax.title.set_verticalalignment("center")
    if legend is not None:
        legend.borderaxespad = 0
        legend.set_bbox_to_anchor(inverse.transform([pixels[:, 0].max(), header_top]))
    fig.canvas.draw()
    boxes = [Bbox.from_extents(*pixels.min(axis=0), *pixels.max(axis=0))]
    for artist in ax.findobj(Text):
        if artist.get_visible() and artist.get_text():
            boxes.append(artist.get_window_extent(renderer))
    if ax.get_legend() is not None:
        boxes.append(ax.get_legend().get_window_extent(renderer))
    crop = Bbox.union(boxes).transformed(fig.dpi_scale_trans.inverted()).padded(0.15)
    kwargs["bbox_inches"] = crop
    save_plot(df.attrs.get("word_label", "RAWs"), fig, path, **kwargs)


def legend_k(ax, k_colors):
    legend_items = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=6, label=f'k={k}')
        for k, c in k_colors.items()
    ]
    ax.legend(handles=legend_items, loc='upper right', title='k-mer length')

# ---- plotting -----------------------------------------------------
def hist_raws(df, out_png, va):
    plt.figure(figsize=(10,6))

    counts = df["kmer_len"].value_counts().sort_index()
    bars = plt.bar(counts.index.astype(str), counts.values, color="darkseagreen")

    # Add numbers above bars
    for bar, val in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width()/2, val + 0.2, str(val),
                ha='center', va='bottom', fontsize=8)
        
    plt.xlabel("k-mer length")
    plt.ylabel("Number of RAWs")
    plt.title(f"Number of RAWs per k-mer [{va}]")
    plt.tight_layout()
    save_plot(df.attrs.get("word_label", "RAWs"), plt.gcf(), f"{out_png}/number_raws.png",
            dpi=300, bbox_inches='tight')    
    plt.close()

def plot_3d_raws(df, out_png, va, k_colors, k_sizes):
    x = df["genome_idx"]
    y = df["position"]
    z = df["kmer_len"]

    colors = df["kmer_len"]
    fig = plt.figure(figsize=(13, 8))
    ax = fig.add_subplot(projection="3d", computed_zorder=False)

    scatter_occurrences(ax, df, k_colors)

    if df.attrs.get("segment_layout"):
        layout = df.attrs["segment_layout"]
        ax.set_xlim3d(0.5, layout["genome_max"] + 0.5)
        ax.set_ylim3d(0, display_extents(df)[0])
    ax.zaxis._axinfo['juggled'] = (1, 2, 0)  # moves the Z-axis pane to the left
    ax.set_xlabel("Genome Index", labelpad=10, fontsize=12)
    ax.set_ylabel("Display position" if df.attrs.get("segment_layout") else "Position in Genome", labelpad=14, fontsize=12)
    ax.set_zlabel("")
    ax.text2D(-0.12, 0.50, "k-mer Length", transform=ax.transAxes, fontsize=12, rotation=90, va="center")
    ax.tick_params(axis="both", labelsize=10, pad=3)
    ax.set_zticks(np.arange(z.min(), z.max() + 1, 1))

    ax.set_box_aspect([1.3, 3, 0.7], zoom=1.06)  # [X, Y, Z] relative scale

    ax.grid(True, linestyle=":", alpha=0.6)
    ax.view_init(elev=15, azim=-30)

    ax.invert_zaxis()

    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.04, top=0.94)
    legend_k(ax, k_colors)
    ax.get_legend().set_bbox_to_anchor((1.20, 0.86))
    ax.set_title(f"3D RAWs Distribution\n{va}", y=0.85, pad=8, fontsize=13)
    # Crop the unused vertical space of the low-elevation 3D axes.
    fig.canvas.draw()
    crop = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    from matplotlib.transforms import Bbox
    crop = Bbox.from_extents(crop.x0 - 1.3, crop.y0 + 0.06 * crop.height,
                            crop.x1 + 1.55, crop.y0 + 0.94 * crop.height)
    save_plot(df.attrs.get("word_label", "RAWs"), fig,
              f"{out_png}/3d_raws_distribution.png", dpi=300, bbox_inches=crop)
    plt.close()

def plot_3d_top_view(df, out_png, va, k_colors, k_sizes):
    x = df["genome_idx"]
    y = df["position"]
    z = df["kmer_len"]
    colors = df["kmer_len"]

    fig = plt.figure(figsize=(18.75, 18.75))
    ax = fig.add_subplot(projection="3d", computed_zorder=False)

    scatter_occurrences(ax, df, k_colors)

    distribution_topdown_axes(ax, x, y, f"Top-Down View of 3D RAWs Distribution [{va}]", df.attrs.get("segment_layout"))

    legend_k(ax, k_colors)
    save_topdown_plot(df, plt.gcf(), f"{out_png}/3d_raws_distribution_topdown.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_view_marked(df, out_png, va, k_colors, k_sizes):
    x = df["genome_idx"]
    y = df["position"]
    z = df["kmer_len"]
    colors = df["kmer_len"]

    mask_base = df["gc_pct"].between(48, 52)   # base GC% range (for dots)

    # --- Top-down (up) view ---
    fig2 = plt.figure(figsize=(18.75, 18.75))
    ax3 = fig2.add_subplot(projection="3d", computed_zorder=False)

    # Base layer (all points)
    scatter_occurrences(ax3, df, k_colors)

    # Overlay hollow rings (highlight only GC 48–52)
    ax3.scatter(
        x[mask_base], y[mask_base], z[mask_base],
        s=45,                     # slightly larger ring than points
        facecolors='none',        # hollow
        edgecolors='orange',         # highlight color
        linewidths=0.8,
        zorder=10
    )

    distribution_topdown_axes(ax3, x, y, f"Top-Down View of 3D RAWs Distribution (GC% 48-52 marked) [{va}]", df.attrs.get("segment_layout"))

    legend_k(ax3, k_colors)
    save_topdown_plot(df, plt.gcf(), f"{out_png}/3d_raws_distribution_topdown_gc48_52_marked.png",
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_view_45_55(df, out_png, va, k_colors, k_sizes):

    mask = df["gc_pct"].between(45, 55)

    x_sub = df.loc[mask, "genome_idx"].values
    y_sub = df.loc[mask, "position"].values
    z_sub = df.loc[mask, "kmer_len"].values

    colors_sub = pd.Series(z_sub).map(k_colors)
    sizes_sub  = pd.Series(z_sub).map(k_sizes)

    # --- New top-down (up) view with ONLY GC 45–55% ---
    fig2 = plt.figure(figsize=(18.75, 18.75))
    ax3 = fig2.add_subplot(projection="3d", computed_zorder=False)

    scatter_occurrences(ax3, df.loc[mask], k_colors)

    distribution_topdown_axes(ax3, x_sub, y_sub, f"Top-Down View of 3D RAWs Distribution (GC% 45-55 only) [{va}]", df.attrs.get("segment_layout"))
    legend_k(ax3, k_colors)
    save_topdown_plot(df, plt.gcf(), f"{out_png}/3d_raws_distribution_topdown_gc45_55.png",
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_view_highlight(df, out_png, va, k_colors, k_sizes):
    mask = df["gc_pct"].between(45, 55)
    x_sub = df.loc[mask, "genome_idx"].values
    y_sub = df.loc[mask, "position"].values
    z_sub = df.loc[mask, "kmer_len"].values

    colors_sub = pd.Series(z_sub).map(k_colors)
    sizes_sub  = pd.Series(z_sub).map(k_sizes)

    # --- NEW: mask for 48–52% ---
    mask_highlight = df["gc_pct"].between(48, 52)
    x_high = df.loc[mask_highlight, "genome_idx"].values
    y_high = df.loc[mask_highlight, "position"].values
    z_high = df.loc[mask_highlight, "kmer_len"].values

    fig2 = plt.figure(figsize=(18.75, 18.75))
    ax3 = fig2.add_subplot(projection="3d", computed_zorder=False)

    ax3.scatter(
        x_high, y_high, z_high,
        c="gold",
        marker="o",
        s=65,
        alpha=0.8,
        edgecolors="none",
        depthshade=False,
        zorder=1
    )

    scatter_occurrences(ax3, df.loc[mask], k_colors)

    distribution_topdown_axes(ax3, x_sub, y_sub, f"Top-Down View of 3D RAWs Distribution (GC% 45-55 only with 48-52% Highlight) [{va}]", df.attrs.get("segment_layout"))
    legend_k(ax3, k_colors)

    save_topdown_plot(df, plt.gcf(), f"{out_png}/3d_raws_distribution_topdown_gc45_55_highlight.png",
                dpi=300, bbox_inches='tight')
    plt.close()

def gc_vs_kmer_length(df, out_png, va, kmer_lengths):
    plt.figure(figsize=(10, 7))
    plt.scatter(pd.to_numeric(df["kmer_len"]), df["gc_pct"], s=5, alpha=0.7)

    plt.xlabel("k-mer Length")
    plt.ylabel("GC Percentage")

    plt.axhspan(45, 55, color='lightgreen', alpha=0.3)

    plt.xticks(kmer_lengths)
    plt.yticks(np.arange(0, 101, 5))  # from 0% to 100% every 5%
    plt.grid(True, linestyle=":", alpha=0.6)

    plt.title(f"GC% vs k-mer Length for RAWs [{va}]")
    save_plot(df.attrs.get("word_label", "RAWs"), plt.gcf(), f"{out_png}/gc%_vs_kmer.png", dpi=300, bbox_inches='tight')
    plt.close()

def results_csv(df, outputs_dir, va):
    prefix = "mraws" if df.attrs.get("word_label") == "mRAWs" else "raws"
    gc_filtered = df[(df["gc_pct"] >= 45) & (df["gc_pct"] <= 55)]
    gc_filtered.to_csv(f"{outputs_dir}/{prefix}_gc_45_55_{va}.csv", sep=",", index=False)

    gc_filtered2 = df[(df["gc_pct"] >= 48) & (df["gc_pct"] <= 52)]
    gc_filtered2.to_csv(f"{outputs_dir}/{prefix}_gc_48_52_{va}.csv", sep=",", index=False)



def clean_variant_name(variant_name):
    """
    Makes dataset names slightly more readable in figure titles.

    Example:
        Begomovirus_manihotis
    becomes:
        Begomovirus manihotis
    """
    return str(variant_name).replace("_", " ")

def make_title(main_title, variant_name, detail=None):
    """
    Standard title format used throughout the figures.
    """
    title = f"{main_title}\n{clean_variant_name(variant_name)}"

    if detail:
        title += f"\n{detail}"

    return title

def plot_top_view_per_k(df, output_dir, variant_name, k_colors, k_sizes, highlight_gc=True):
    """
    Generate one top-down positional distribution figure for each k-mer length.

    All plots deliberately use the FULL dataset x/y limits, not the subset
    limits. This makes the different k values directly comparable.
    """
    x_all = df['genome_idx']
    y_all = df['position']
    kmer_lengths = sorted(df['kmer_len'].unique())
    for k in kmer_lengths:
        subset = df[df['kmer_len'] == k].copy()
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot()
        scatter_occurrences(ax, subset, k_colors, ordinary_size=2)
        if highlight_gc:
            mask_gc = subset['gc_pct'].between(48, 52)
            ax.scatter(subset.loc[mask_gc, 'position'], subset.loc[mask_gc, 'genome_idx'],
                       s=45, facecolors='none', edgecolors='orange',
                       linewidths=0.8, zorder=10)
            ax.legend(handles=[Line2D([], [], marker='o', linestyle='None',
                      markerfacecolor='none', markeredgecolor='orange',
                      markeredgewidth=0.8, markersize=np.sqrt(45), label='GC 48–52%')],
                      loc='lower right', bbox_to_anchor=(1.0, 1.015),
                      borderaxespad=0, borderpad=0, frameon=False, fontsize=14)
        topdown_axes(ax, y_all, x_all, make_title('RAW positional distribution', variant_name, f'k = {int(k)}'), df.attrs.get('segment_layout'))
        ax.text(0.0, 1.015, f'n = {len(subset):,} RAWs', transform=ax.transAxes, ha='left', va='bottom', fontsize=14)
        save_plot(df.attrs.get('word_label', 'RAWs'), fig, os.path.join(output_dir, f'raw_position_topdown_k{int(k)}.png'), dpi=300, bbox_inches='tight')
        save_interactive_topdown(df, subset, output_dir, variant_name, k, k_colors[k], highlight_gc)

def save_interactive_topdown(df, subset, output_dir, variant_name, k, color, highlight_gc=True):
    """Fill the HTML template with the points for one k-mer length."""
    label = df.attrs.get("word_label", "RAWs")
    columns = ["position", "genome_idx", "raw", "gc_pct", "kmer_len"]
    extra_columns = ["component", "segment_position", "accession"] if df.attrs.get("segment_layout") else []
    xmax, ymax = display_extents(df)
    payload = {
        "rows": subset[columns + extra_columns].values.tolist(),
        "xmax": xmax, "ymax": ymax,
        "color": color, "highlight": highlight_gc, "label": label,
        "position_label": "Display position" if extra_columns else "Position in Genome",
        "extra_columns": ["Component", "Position in component", "Accession"] if extra_columns else [],
    }
    # Escape HTML script delimiters even if an input word contains unexpected text.
    data = json.dumps(payload, ensure_ascii=True).replace("<", "\\u003c")
    title = html.escape(f"{label} positional distribution — {variant_name} — k = {int(k)}")
    template_path = Path(__file__).with_name("raw_position_topdown.html")
    page = template_path.read_text(encoding="utf-8")
    page = page.replace("__TITLE__", title).replace("__DATA__", data)
    Path(output_dir, f"raw_position_topdown_k{int(k)}.html").write_text(page, encoding="utf-8")


def window_profiles(counts, gc_sums, window):
    kernel = np.ones(window, dtype=float)
    raw_counts = np.convolve(counts, kernel, mode="full")[window // 2:window // 2 + len(counts)]
    window_gc_sums = np.convolve(gc_sums, kernel, mode="full")[window // 2:window // 2 + len(counts)]
    mean_gc = np.divide(
        window_gc_sums, raw_counts,
        out=np.full(len(counts), np.nan), where=raw_counts > 0,
    )
    return raw_counts, mean_gc


def positional_profiles(df, window=51):
    """Pool profiles, clipping each record's window at its real component ends."""
    layout = df.attrs.get("segment_layout")
    end = layout["position_max"] if layout else int(df.position.max())
    positions = np.arange(1, end + 1)
    indices = df.position.to_numpy(dtype=int) - 1
    targets = np.bincount(indices[df.gc_pct.between(48, 52)], minlength=end).astype(float)
    if not layout:
        counts = np.bincount(indices, minlength=end).astype(float)
        gc_sums = np.bincount(indices, weights=df.gc_pct, minlength=end)
        window_counts, mean_gc = window_profiles(counts, gc_sums, window)
        return positions, window_counts, mean_gc, targets

    window_counts = np.zeros(end)
    window_gc = np.zeros(end)
    covered = np.zeros(end, dtype=bool)
    groups = {int(record): rows for record, rows in df.groupby("record_idx", sort=False)}
    for record in layout["records"]:
        offset = record["b_offset"] if record["component"] == "B" else 0
        length = record["length"]
        span = slice(offset, offset + length)
        covered[span] = True
        rows = groups.get(record["record_idx"])
        if rows is None:
            continue
        local_indices = rows.segment_position.to_numpy(dtype=int) - 1
        counts = np.bincount(local_indices, minlength=length).astype(float)
        gc_sums = np.bincount(local_indices, weights=rows.gc_pct, minlength=length)
        counts, means = window_profiles(counts, gc_sums, window)
        window_counts[span] += counts
        window_gc[span] += np.nan_to_num(means) * counts
    mean_gc = np.divide(window_gc, window_counts, out=np.full(end, np.nan), where=window_counts > 0)
    window_counts[~covered] = np.nan
    targets[~covered] = np.nan
    return positions, window_counts, mean_gc, targets


def plot_gradient_profile(
    positions, raw_counts, mean_gc, target_counts, output_path, title, window, word_label="RAWs",
    position_label="Position in genome",
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
        xlabel=position_label, ylabel="RAWs in window",
        xlim=(positions[0], positions[-1]), ylim=(0, None),
    )
    raw_ax.xaxis.set_major_locator(MultipleLocator(500 * max(1, int(np.ceil(positions[-1] / (500 * 15))))))
    gc_ax.set(ylabel="Mean RAW GC content (%)", ylim=(100, 0))
    raw_ax.grid(True, alpha=0.2)
    fig.suptitle(f"{title}: positional RAW and GC profile")
    raw_ax.set_title(f"{window}-position moving window", fontsize=10)
    lines1, labels1 = raw_ax.get_legend_handles_labels()
    lines2, labels2 = gc_ax.get_legend_handles_labels()
    raw_ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right")
    fig.subplots_adjust(left=0.07, right=0.82, bottom=0.12, top=0.88)
    save_plot(word_label, fig, output_path, dpi=200)
    plt.close(fig)

def plot_exact_target_counts(
    positions, raw_counts, target_counts, output_path, title, window, word_label="RAWs",
    position_label="Position in genome",
):
    fig, (raw_ax, target_ax) = plt.subplots(
        2, 1, sharex=True, figsize=(15, 7),
        gridspec_kw={"height_ratios": [3, 1.5], "hspace": 0.08},
    )
    raw_ax.plot(positions, raw_counts, color="#2563a6", linewidth=1.7)
    target_ax.plot(positions, target_counts, color="#b91c1c", linewidth=1.2)
    raw_ax.set_ylabel("All RAWs in window")
    target_ax.set_ylabel("48–52% GC RAWs\nat exact position")
    target_ax.set_xlabel(position_label)
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
    save_plot(word_label, fig, output_path, dpi=200)
    plt.close(fig)

def run_plots(input_csv, outputs_dir, pathogen_name, host_name, word_label="RAWs", metadata_path=None):
    if word_label not in ("RAWs", "mRAWs"):
        raise ValueError("word_label must be RAWs or mRAWs")
    ensure_outdir(outputs_dir)
    df = pd.read_csv(input_csv, header=None,
                     names=["genome_idx", "position", "raw", "gc_pct", "kmer_len"])
    if df.empty:
        print(f"No {word_label} records in {input_csv}; skipping plots.")
        return
    for column in ("genome_idx", "position", "gc_pct", "kmer_len"):
        df[column] = pd.to_numeric(df[column], errors="raise")
    if not np.isfinite(df[["genome_idx", "position", "gc_pct", "kmer_len"]]).all().all():
        raise ValueError("Numeric fields must be finite")
    for column in ("genome_idx", "position", "kmer_len"):
        if (df[column] % 1 != 0).any():
            raise ValueError(f"{column} must contain integers")
    if (df.position < 1).any() or (df.genome_idx < 0).any() or not df.gc_pct.between(0, 100).all():
        raise ValueError("Invalid coordinate or GC percentage")
    if not (df.raw.str.len() == df.kmer_len).all():
        raise ValueError("Word lengths do not match kmer_len")
    if metadata_path is not None:
        df = apply_segment_layout(df, metadata_path)
    df.attrs["word_label"] = word_label
    title = f"Pathogen: {pathogen_name.replace('_', ' ')} | Host: {host_name.replace('_', ' ')}"
    lengths, colors, sizes = kmer_styles(df)
    hist_raws(df, outputs_dir, title)
    _, distribution_colors, distribution_sizes = distribution_kmer_styles(df)
    for plot in (plot_3d_raws, plot_3d_top_view, plot_top_view_marked):
        plot(df, outputs_dir, title, distribution_colors, distribution_sizes)
    plot_top_view_per_k(df, outputs_dir, title, colors, sizes)
    if df.gc_pct.between(45, 55).any():
        for plot in (plot_top_view_45_55, plot_top_view_highlight):
            plot(df, outputs_dir, title, distribution_colors, distribution_sizes)
    gc_vs_kmer_length(df, outputs_dir, title, lengths)
    slug = "".join(c if c.isalnum() or c in "-_" else "_"
                   for c in f"{pathogen_name}_{host_name}")
    results_csv(df, outputs_dir, slug)
    positions, window_counts, mean_gc, targets = positional_profiles(df, 51)
    position_label = "Display position" if metadata_path is not None else "Position in genome"
    plot_gradient_profile(positions, window_counts, mean_gc, targets,
        Path(outputs_dir) / "raw_gc_profile_window51_exact_target_count_gradient.png",
        title, 51, word_label, position_label)
    plot_exact_target_counts(positions, window_counts, targets,
        Path(outputs_dir) / "raw_window51_target_count_window1.png", title, 51, word_label, position_label)
    print(f"Saved {word_label} plots to {outputs_dir}: {len(df):,} occurrences")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv")
    parser.add_argument("output_dir")
    parser.add_argument("--pathogen-name", "--virus-name", dest="pathogen_name", required=True)
    parser.add_argument("--host-name", required=True)
    parser.add_argument("--word-label", choices=("RAWs", "mRAWs"), default="RAWs")
    parser.add_argument("--metadata", type=Path, help="Component layout TSV from JoinSegments.py")
    args = parser.parse_args()
    run_plots(args.input_csv, args.output_dir, args.pathogen_name,
              args.host_name, args.word_label, args.metadata)


if __name__ == "__main__":
    main()