#!/usr/bin/env python
"""
24c_pseudotime_violin.py
Fig 2H-style violin plot: DPT pseudotime by group, ranked by median pseudotime.

Groups:
  - Gurtner Bone Marrow (most immature)
  - Gurtner Blood
  - Qiu eos clusters C0-C6
  - Hu adipose types: eWAT, scWAT, BAT (parsed from CellName prefix)

Input:  three-way CellRank scores + metadata from script 24b
Output: violin plot PDF
SLURM:  batch, 16G, 2 CPUs, 30min, cellrank_env
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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

# --- Color palettes (match 00_theme.R) ---------------------------------------
EOS_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4", "C6": "#2F3C7E",
}

EXTRA_COLORS = {
    "Gurtner\nBone Marrow": "#2E8B57",
    "Gurtner\nBlood": "#DAA520",
    "Hu eWAT": "#FF8C00",
    "Hu scWAT": "#8B4513",
    "Hu BAT": "#DC143C",
}

# --- Paths -------------------------------------------------------------------
BASE = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
THREE_WAY = os.path.join(BASE, "integration", "three_way")
META_PATH = os.path.join(THREE_WAY, "for_cellrank", "metadata.csv")
SCORES_PATH = os.path.join(THREE_WAY, "cellrank_scores.csv")
FIG_DIR = os.path.join(THREE_WAY, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# --- Load data ---------------------------------------------------------------
print("Loading metadata and CellRank scores ...")

meta = pd.read_csv(META_PATH, index_col=0)
scores = pd.read_csv(SCORES_PATH, index_col=0)

# Merge
df = meta.join(scores[["ct_score", "dpt_pseudotime"]], how="inner")
print(f"  Merged: {len(df)} cells")

# --- Assign violin groups ----------------------------------------------------
print("Assigning violin groups ...")

def assign_group(row):
    ds = row.get("dataset", "")
    organ = row.get("organ", "")
    cellname = str(row.get("CellName", ""))
    cluster = row.get("original_cluster", "")

    if ds == "Gurtner_2022_Blood_BM":
        if organ == "Bone Marrow":
            return "Gurtner\nBone Marrow"
        elif organ == "Blood":
            return "Gurtner\nBlood"
        else:
            return None  # exclude other Gurtner organs
    elif ds == "Hu_2026_Healthy":
        # Parse tissue type from CellName prefix (e.g., BAT_1000383 -> BAT)
        prefix = cellname.split("_")[0] if "_" in cellname else cellname
        if prefix in ("eWAT", "ewat"):
            return "Hu eWAT"
        elif prefix in ("scWAT", "scwat"):
            return "Hu scWAT"
        elif prefix in ("BAT", "bat"):
            return "Hu BAT"
        else:
            return f"Hu {prefix}"
    elif ds == "Qiu_2025_Adipose":
        return cluster  # C0, C1, ..., C6
    else:
        return None

df["violin_group"] = df.apply(assign_group, axis=1)
df = df.dropna(subset=["violin_group", "dpt_pseudotime"])
df = df[df["dpt_pseudotime"] < 1.0]  # remove inf/unreachable

print(f"  Cells with valid groups: {len(df)}")
print(f"  Groups: {df['violin_group'].value_counts().to_dict()}")

# --- Order by median pseudotime ----------------------------------------------
group_medians = df.groupby("violin_group")["dpt_pseudotime"].median().sort_values()
ordered_groups = group_medians.index.tolist()
print(f"\n  Group order by median DPT:")
for g in ordered_groups:
    n = (df["violin_group"] == g).sum()
    print(f"    {g:20s}  median={group_medians[g]:.3f}  n={n}")

# --- Build color map ---------------------------------------------------------
all_colors = {}
all_colors.update(EXTRA_COLORS)
all_colors.update(EOS_COLORS)
colors = [all_colors.get(g, "#999999") for g in ordered_groups]

# --- Violin plot (Fig 2H style) ----------------------------------------------
print("\nGenerating violin plot ...")

fig, ax = plt.subplots(figsize=(max(10, len(ordered_groups) * 0.9 + 2), 6))

data_list = [df.loc[df["violin_group"] == g, "dpt_pseudotime"].values
             for g in ordered_groups]

parts = ax.violinplot(data_list, positions=range(len(ordered_groups)),
                      showmeans=False, showmedians=False, showextrema=False)

# Color each violin
for i, pc in enumerate(parts["bodies"]):
    pc.set_facecolor(colors[i])
    pc.set_edgecolor("black")
    pc.set_linewidth(0.8)
    pc.set_alpha(0.85)

# Add boxplot inside (thin, white median line)
bp = ax.boxplot(data_list, positions=range(len(ordered_groups)),
                widths=0.15, showfliers=False, patch_artist=True,
                medianprops=dict(color="white", linewidth=2),
                boxprops=dict(facecolor="black", linewidth=0.8),
                whiskerprops=dict(color="black", linewidth=0.8),
                capprops=dict(color="black", linewidth=0.8))

# Style
ax.set_xticks(range(len(ordered_groups)))
ax.set_xticklabels(ordered_groups, rotation=45, ha="right", fontsize=10,
                   fontweight="bold")
ax.set_ylabel("Pseudotime", fontsize=13, fontweight="bold")
ax.set_title("Diffusion Pseudotime by Group", fontsize=14, fontweight="bold")
ax.set_ylim(-0.02, 1.02)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

out_path = os.path.join(FIG_DIR, "24c_pseudotime_violin.pdf")
fig.savefig(out_path)
plt.close()
print(f"\n  Saved: {out_path}")

# Also save group summary
summary = df.groupby("violin_group")["dpt_pseudotime"].agg(
    ["count", "mean", "median", "std", "min", "max"]
).reindex(ordered_groups)
summary.to_csv(os.path.join(THREE_WAY, "pseudotime_violin_summary.csv"))
print(f"  Saved: pseudotime_violin_summary.csv")

print("\nDone.")
