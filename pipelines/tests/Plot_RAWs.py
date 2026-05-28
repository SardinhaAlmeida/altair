import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# ---- helpers ---------------------------------------------------------------
k_colors = {11: "red", 12: "blue", 13: "darkseagreen", 14: "grey"}
k_sizes = {11: 3, 12: 1, 13: 1, 14: 1}

def ensure_outdir(d):
    os.makedirs(d, exist_ok=True)

def topdown_axes(ax, x, y, title):
    ax.set_title(title)
    ax.set_xlabel("Genome Index")
    ax.set_ylabel("Position", labelpad=15)
    ax.set_xlim3d(x.min(), x.max())
    ax.set_ylim3d(y.min(), y.max())
    if len(y): ax.set_yticks(np.arange(0, y.max() + 250, 250))
    if len(x): ax.set_xticks(np.arange(0, x.max(), 5))
    ax.invert_zaxis()
    ax.view_init(elev=90, azim=180)  # top-down
    ax.set_proj_type('ortho')
    ax.set_box_aspect([6, 6, 1])
    ax.yaxis.set_label_coords(-0.05, 0.5)
    ax.tick_params(axis='y', labelrotation=90)

def legend_k(ax):
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
    plt.savefig(f"{out_png}/number_raws.png",
            dpi=300, bbox_inches='tight')    
    plt.close()

def plot_3d_raws(df, out_png, va):
    x = df["genome_idx"]
    y = df["position"]
    z = df["kmer_len"]

    colors = df["kmer_len"]
    k_colors = {11: "red", 12: "blue", 13: "darkseagreen" , 14: "grey"}
    k_sizes = {11: 3, 12: 1, 13: 1, 14: 1}

    fig = plt.figure(figsize=(18,10))
    ax = fig.add_subplot(projection="3d")

    ax.scatter(
        x, y, z,
        c=colors.map(k_colors),
        s=colors.map(k_sizes),               # dot size
        alpha=0.9,
        depthshade=False
    )

    ax.zaxis._axinfo['juggled'] = (1, 2, 0)  # moves the Z-axis pane to the left
    ax.set_xlabel("Genome Index")
    ax.set_ylabel("Position in Genome", labelpad=15)
    ax.set_zlabel("k-mer Length")
    ax.set_zticks(np.arange(z.min(), z.max() + 1, 1))

    ax.set_box_aspect([1.3, 3, 0.5])  # [X, Y, Z] relative scale

    ax.grid(True, linestyle=":", alpha=0.6)
    ax.view_init(elev=15, azim=-30)

    ax.invert_zaxis()

    legend_k(ax)

    plt.title(f"3D RAWs Distribution [{va}]", pad=20)
    plt.savefig(f"{out_png}3d_raws_distribution.png", dpi=300)
    plt.close()

def plot_3d_top_view(df, out_png, va):
    x = df["genome_idx"]
    y = df["position"]
    z = df["kmer_len"]
    colors = df["kmer_len"]

    fig = plt.figure(figsize=(30, 20))
    ax = fig.add_subplot(projection="3d")

    ax.scatter(
        x, y, z,
        c=colors.map(k_colors),
        s=colors.map(k_sizes),               # dot size
        alpha=0.9
    )

    topdown_axes(ax, x, y, f"Top-Down View of 3D RAWs Distribution [{va}]")

    legend_k(ax)
    plt.savefig(f"{out_png}3d_raws_distribution_topdown.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_view_marked(df, out_png, va):
    x = df["genome_idx"]
    y = df["position"]
    z = df["kmer_len"]
    colors = df["kmer_len"]

    mask_base = df["gc_pct"].between(48, 52)   # base GC% range (for dots)

    # --- Top-down (up) view ---
    fig2 = plt.figure(figsize=(15, 35))
    ax3 = fig2.add_subplot(projection="3d")

    # Base layer (all points)
    ax3.scatter(
        x, y, z,
        c=colors.map(k_colors),
        s=colors.map(k_sizes),
        alpha=0.9
    )

    # Overlay hollow rings (highlight only GC 48–52)
    ax3.scatter(
        x[mask_base], y[mask_base], z[mask_base],
        s=45,                     # slightly larger ring than points
        facecolors='none',        # hollow
        edgecolors='black',         # highlight color
        linewidths=0.8,
        zorder=10
    )

    topdown_axes(ax3, x, y, f"Top-Down View of 3D RAWs Distribution (GC% 48-52 marked) [{va}]")

    legend_k(ax3)
    plt.savefig(f"{out_png}/3d_raws_distribution_topdown_gc48_52_marked.png",
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_view_45_55(df, out_png, va):

    mask = df["gc_pct"].between(45, 55)

    x_sub = df.loc[mask, "genome_idx"].values
    y_sub = df.loc[mask, "position"].values
    z_sub = df.loc[mask, "kmer_len"].values

    colors_sub = pd.Series(z_sub).map(k_colors)
    sizes_sub  = pd.Series(z_sub).map(k_sizes)

    # --- New top-down (up) view with ONLY GC 45–55% ---
    fig2 = plt.figure(figsize=(15, 15))
    ax3 = fig2.add_subplot(projection="3d")

    ax3.scatter(
        x_sub, y_sub, z_sub,
        c=colors_sub,
        s=sizes_sub,
        alpha=0.9
    )

    topdown_axes(ax3, x_sub, y_sub, f"Top-Down View of 3D RAWs Distribution (GC% 45-55 only) [{va}]")
    legend_k(ax3)
    plt.savefig(f"{out_png}/3d_raws_distribution_topdown_gc45_55.png",
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_top_view_highlight(df, out_png, va):
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

    fig2 = plt.figure(figsize=(15, 15))
    ax3 = fig2.add_subplot(projection="3d")

    ax3.scatter(
        x_high, y_high, z_high,
        c="gold",  # soft highlight color (you can change it)
        s=50,      # larger than normal dots
        alpha=0.25,
        edgecolors="none"
    )

    ax3.scatter(
        x_sub, y_sub, z_sub,
        c=colors_sub,
        s=sizes_sub,
        alpha=0.9
    )

    topdown_axes(ax3, x_sub, y_sub, f"Top-Down View of 3D RAWs Distribution (GC% 45-55 only with 48-52% Highlight) [{va}]")
    legend_k(ax3)

    plt.savefig(f"{out_png}/3d_raws_distribution_topdown_gc45_55_highlight.png",
                dpi=300, bbox_inches='tight')
    plt.close()

def gc_vs_kmer_length(df, out_png, va):
    plt.figure(figsize=(10, 7))
    plt.scatter(pd.to_numeric(df["kmer_len"]), df["gc_pct"], s=5, alpha=0.7)

    plt.xlabel("k-mer Length")
    plt.ylabel("GC Percentage")

    plt.axhspan(45, 55, color='lightgreen', alpha=0.3)

    plt.xticks(np.arange(11, 15, 1))  # from 0% to 100% every 5%
    plt.yticks(np.arange(0, 101, 5))  # from 0% to 100% every 5%
    plt.grid(True, linestyle=":", alpha=0.6)

    plt.title(f"GC% vs k-mer Length for RAWs [{va}]")
    plt.savefig(f"{out_png}gc%_vs_kmer.png", dpi=300, bbox_inches='tight')
    plt.close()

def results_csv(df, outputs_dir, va):
    gc_filtered = df[(df["gc_pct"] >= 45) & (df["gc_pct"] <= 55)]
    gc_filtered.to_csv(f"{outputs_dir}raws_gc_45_55_{va}.csv", sep=",", index=False)

    gc_filtered2 = df[(df["gc_pct"] >= 48) & (df["gc_pct"] <= 52)]
    gc_filtered2.to_csv(f"{outputs_dir}raws_gc_48_52_{va}.csv", sep=",", index=False)

# ---- running plots ------------------------------------------------------------------
def run_plots(input_csv: str, outputs_dir: str, variant_name: str):

    ensure_outdir(outputs_dir)

    df = pd.read_csv(
        input_csv,
        header=None,
        names=["genome_idx", "position", "raw", "gc_pct", "kmer_len"]
    )
    
    hist_raws(df, outputs_dir, variant_name)
    plot_3d_raws(df, outputs_dir, variant_name)
    plot_3d_top_view(df, outputs_dir, variant_name)
    plot_top_view_marked(df, outputs_dir, variant_name)
    plot_top_view_45_55(df, outputs_dir, variant_name)
    plot_top_view_highlight(df, outputs_dir, variant_name)
    gc_vs_kmer_length(df, outputs_dir, variant_name)
    results_csv(df, outputs_dir, variant_name)
