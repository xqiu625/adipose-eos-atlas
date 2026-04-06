#!/usr/bin/env python
"""
24c_pseudotime_violin.py
Violin plots for CytoTRACE + DPT pseudotime from three-way integration.
Groups: Gurtner BM, Gurtner Blood, Hu BAT/eWAT/scWAT, Qiu C0-C5

Input:  cellrank_scores.csv + metadata.csv from scripts 24/24b
        hu_adipose_eos.rds metadata for depot info
Output: 24c_cytotrace_violin.pdf, 24c_dpt_violin.pdf
SLURM:  batch, 16G, 4 CPUs, 30min, cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Publication figure settings ---------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# --- Color palettes -----------------------------------------------------------
EOS_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4",
}

EXTRA_COLORS = {
    "Gurtner\nBone Marrow": "#2E8B57",
    "Gurtner\nBlood": "#DAA520",
    "Hu eWAT": "#FF8C00",
    "Hu scWAT": "#8B4513",
    "Hu BAT": "#DC143C",
}

# --- Paths --------------------------------------------------------------------
BASE = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
THREE_WAY = os.path.join(BASE, "integration", "three_way")
META_PATH = os.path.join(THREE_WAY, "for_cellrank", "metadata.csv")
SCORES_PATH = os.path.join(THREE_WAY, "cellrank_scores.csv")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# --- Load data ----------------------------------------------------------------
print("Loading metadata and CytoTRACE/DPT scores ...")
meta = pd.read_csv(META_PATH, index_col=0)
scores = pd.read_csv(SCORES_PATH, index_col=0)

print(f"  Metadata: {len(meta)} cells, columns: {meta.columns.tolist()}")
print(f"  Scores: {len(scores)} cells, columns: {scores.columns.tolist()}")

# Merge scores into metadata
score_cols = [c for c in ["ct_score", "dpt_pseudotime"] if c in scores.columns]
df = meta.join(scores[score_cols], how="inner")
print(f"  Merged: {len(df)} cells")

has_dpt = "dpt_pseudotime" in df.columns
has_ct = "ct_score" in df.columns

# --- Load Hu depot info from original object ----------------------------------
# The three-way integration set hu$organ = "Adipose_Healthy" uniformly.
# Get depot (BAT/eWAT/scWAT) from the 'depot' column if it exists,
# or fall back to parsing from barcodes.
print("Resolving Hu depot info ...")

hu_mask = df["dataset"].astype(str).str.contains("Hu", na=False)
print(f"  Hu cells in three-way: {hu_mask.sum()}")

if "depot" in df.columns:
    print("  Using 'depot' column from metadata")
else:
    # Try to load depot mapping from original Hu metadata
    depot_loaded = False
    hu_meta_candidates = [
        os.path.join(BASE, "integration", "three_way", "for_cellrank", "metadata.csv"),
    ]
    # Check if organ column has depot detail for Hu cells
    hu_organs = df.loc[hu_mask, "organ"].unique()
    print(f"  Hu organ values: {hu_organs}")
    if any(x in str(hu_organs) for x in ["BAT", "eWAT", "scWAT"]):
        df["depot"] = df["organ"]
        depot_loaded = True
        print("  Depot info found in organ column")

    if not depot_loaded:
        # Parse from barcode: Hu cells have "Hu_<barcode>" format
        # Original Hu barcodes may have depot prefix: BAT_xxx, eWAT_xxx, scWAT_xxx
        print("  Attempting to parse depot from cell barcodes ...")
        def parse_depot(idx):
            # Remove "Hu_" prefix added by Seurat merge
            bc = str(idx)
            if bc.startswith("Hu_"):
                bc = bc[3:]
            parts = bc.split("_")
            if len(parts) >= 2:
                prefix = parts[0].upper()
                if prefix in ["BAT", "EWAT", "SCWAT"]:
                    return prefix.replace("EWAT", "eWAT").replace("SCWAT", "scWAT").replace("BAT", "BAT")
            return "Unknown"

        df.loc[hu_mask, "depot"] = [parse_depot(idx) for idx in df.index[hu_mask]]
        depot_vals = df.loc[hu_mask, "depot"].value_counts()
        print(f"  Parsed depots: {depot_vals.to_dict()}")

# --- Assign violin groups -----------------------------------------------------
print("Assigning violin groups ...")

def assign_group(row):
    ds = str(row.get("dataset", ""))
    organ = str(row.get("organ", ""))
    cluster = str(row.get("original_cluster", ""))
    depot = str(row.get("depot", "")) if "depot" in row.index else ""

    if "Gurtner" in ds:
        if "Bone Marrow" in organ:
            return "Gurtner\nBone Marrow"
        elif "Blood" in organ:
            return "Gurtner\nBlood"
        else:
            return None
    elif "Hu" in ds:
        if "eWAT" in depot:
            return "Hu eWAT"
        elif "scWAT" in depot:
            return "Hu scWAT"
        elif "BAT" in depot:
            return "Hu BAT"
        else:
            return "Hu Adipose"
    elif "Qiu" in ds:
        # Use cluster name without prefix
        cl = cluster.replace("Eos_", "")
        if cl in ["C0", "C1", "C2", "C3", "C4", "C5"]:
            return cl
        return None
    return None

df["violin_group"] = df.apply(assign_group, axis=1)
df = df.dropna(subset=["violin_group"])
print(f"  Cells with valid groups: {len(df)}")
print(f"  Groups: {df['violin_group'].value_counts().to_dict()}")


# --- Helper: make violin plot -------------------------------------------------
def make_violin(df, score_col, title, ylabel, out_path, ylim=None):
    """Create violin + boxplot for a score column, ordered by median."""
    plot_df = df.dropna(subset=[score_col]).copy()
    if score_col == "dpt_pseudotime":
        plot_df = plot_df[plot_df[score_col] < 1.0]  # remove inf

    group_medians = plot_df.groupby("violin_group")[score_col].median().sort_values(ascending=False)
    ordered_groups = group_medians.index.tolist()

    all_colors = {}
    all_colors.update(EXTRA_COLORS)
    all_colors.update(EOS_COLORS)
    colors = [all_colors.get(g, "#999999") for g in ordered_groups]

    fig, ax = plt.subplots(figsize=(max(12, len(ordered_groups) * 0.9 + 2), 6))

    data_list = [plot_df.loc[plot_df["violin_group"] == g, score_col].values
                 for g in ordered_groups]

    parts = ax.violinplot(data_list, positions=range(len(ordered_groups)),
                          showmeans=False, showmedians=False, showextrema=False)

    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor("black")
        pc.set_linewidth(0.8)
        pc.set_alpha(0.85)

    ax.boxplot(data_list, positions=range(len(ordered_groups)),
               widths=0.15, showfliers=False, patch_artist=True,
               medianprops=dict(color="white", linewidth=2),
               boxprops=dict(facecolor="black", linewidth=0.8),
               whiskerprops=dict(color="black", linewidth=0.8),
               capprops=dict(color="black", linewidth=0.8))

    ax.set_xticks(range(len(ordered_groups)))
    ax.set_xticklabels(ordered_groups, rotation=45, ha="right", fontsize=10,
                       fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=13, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    if ylim:
        ax.set_ylim(ylim)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path)
    plt.close()
    print(f"  Saved: {out_path}")

    # Summary stats
    summary = plot_df.groupby("violin_group")[score_col].agg(
        ["count", "mean", "median", "std"]).reindex(ordered_groups)
    return summary


# --- CytoTRACE violin --------------------------------------------------------
if has_ct:
    print("\nGenerating CytoTRACE violin ...")
    ct_summary = make_violin(
        df, "ct_score",
        title="CytoTRACE Score by Group\n(high = progenitor-like)",
        ylabel="CytoTRACE Score",
        out_path=os.path.join(FIG_DIR, "24c_cytotrace_violin.pdf"),
        ylim=(-0.02, 1.02)
    )
    ct_summary.to_csv(os.path.join(THREE_WAY, "cytotrace_violin_summary.csv"))
    print(ct_summary)

# --- DPT pseudotime violin ---------------------------------------------------
if has_dpt:
    print("\nGenerating DPT pseudotime violin ...")
    dpt_summary = make_violin(
        df, "dpt_pseudotime",
        title="Diffusion Pseudotime by Group\n(BM root → tissue maturation)",
        ylabel="Diffusion Pseudotime",
        out_path=os.path.join(FIG_DIR, "24c_dpt_violin.pdf")
    )
    dpt_summary.to_csv(os.path.join(THREE_WAY, "dpt_violin_summary.csv"))
    print(dpt_summary)

print("\nDone.")
