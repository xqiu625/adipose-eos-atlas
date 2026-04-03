#!/usr/bin/env python
"""
24c_pseudotime_violin.py
Fig 2H-style violin plot: CytoTRACE score by group, ranked by median score.
Higher CytoTRACE = more progenitor-like (immature).

Uses THREE-WAY integration (Qiu + Gurtner + Hu) CytoTRACE scores from 24b.

Groups:
  - Gurtner Bone Marrow (most immature, highest CytoTRACE)
  - Gurtner Blood
  - Hu adipose types: eWAT, scWAT, BAT
  - Qiu eos clusters C0-C6

Input:  three-way CellRank scores + metadata from script 24b
Output: violin plot PDF
SLURM:  batch, 128G, 4 CPUs, 30min, cellrank_env
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

# --- Paths (three-way integration from script 24b) ---------------------------
BASE = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
THREE_WAY = os.path.join(BASE, "integration", "three_way")
META_PATH = os.path.join(THREE_WAY, "for_cellrank", "metadata.csv")
SCORES_PATH = os.path.join(THREE_WAY, "cellrank_scores.csv")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# --- Load data ---------------------------------------------------------------
print("Loading three-way integration metadata and CytoTRACE scores ...")

meta = pd.read_csv(META_PATH, index_col=0)
scores = pd.read_csv(SCORES_PATH, index_col=0)

print(f"  Metadata: {len(meta)} cells")
print(f"  Metadata columns: {meta.columns.tolist()}")
print(f"  Scores: {len(scores)} cells")
print(f"  Scores columns: {scores.columns.tolist()}")

# Merge
if "ct_score" not in scores.columns:
    print("ERROR: ct_score column not found. Make sure script 24b completed.")
    raise SystemExit(1)

merge_cols = ["ct_score"]
df = meta.join(scores[merge_cols], how="inner")
print(f"  Merged: {len(df)} cells")

# --- Assign violin groups ----------------------------------------------------
print("Assigning violin groups ...")

def assign_group(row):
    ds = str(row.get("dataset", ""))
    organ = str(row.get("organ", ""))
    cellname = str(row.get("CellName", row.name))
    cluster = str(row.get("original_cluster", ""))

    if "Gurtner" in ds:
        if organ == "Bone Marrow":
            return "Gurtner\nBone Marrow"
        elif organ == "Blood":
            return "Gurtner\nBlood"
        else:
            return None  # exclude other Gurtner organs
    elif "Hu" in ds:
        # Parse tissue type from CellName prefix (e.g., BAT_1000383 -> BAT)
        prefix = cellname.split("_")[0] if "_" in cellname else cellname
        if prefix.upper() == "EWAT":
            return "Hu eWAT"
        elif prefix.upper() == "SCWAT":
            return "Hu scWAT"
        elif prefix.upper() == "BAT":
            return "Hu BAT"
        else:
            return f"Hu {prefix}"
    elif "Qiu" in ds:
        return cluster  # C0, C1, ..., C6
    else:
        return None

df["violin_group"] = df.apply(assign_group, axis=1)
df = df.dropna(subset=["violin_group", "ct_score"])

print(f"  Cells with valid groups: {len(df)}")
print(f"  Groups: {df['violin_group'].value_counts().to_dict()}")

# --- Order by median CytoTRACE (descending: most progenitor first) ----------
group_medians = df.groupby("violin_group")["ct_score"].median().sort_values(ascending=False)
ordered_groups = group_medians.index.tolist()
print(f"\n  Group order by median CytoTRACE (high = progenitor):")
for g in ordered_groups:
    n = (df["violin_group"] == g).sum()
    print(f"    {g:20s}  median={group_medians[g]:.4f}  n={n}")

# --- Build color map ---------------------------------------------------------
all_colors = {}
all_colors.update(EXTRA_COLORS)
all_colors.update(EOS_COLORS)
colors = [all_colors.get(g, "#999999") for g in ordered_groups]

# --- Violin plot (Fig 2H style) ----------------------------------------------
print("\nGenerating violin plot ...")

fig, ax = plt.subplots(figsize=(max(12, len(ordered_groups) * 0.9 + 2), 6))

data_list = [df.loc[df["violin_group"] == g, "ct_score"].values
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
ax.set_ylabel("CytoTRACE Score", fontsize=13, fontweight="bold")
ax.set_title("CytoTRACE Score by Group\n(high = progenitor-like)",
             fontsize=14, fontweight="bold")
ax.set_ylim(-0.02, 1.02)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

out_path = os.path.join(FIG_DIR, "24c_cytotrace_violin.pdf")
fig.savefig(out_path)
plt.close()
print(f"\n  Saved: {out_path}")

# Also save group summary
summary = df.groupby("violin_group")["ct_score"].agg(
    ["count", "mean", "median", "std", "min", "max"]
).reindex(ordered_groups)
summary.to_csv(os.path.join(THREE_WAY, "cytotrace_violin_summary.csv"))
print(f"  Saved: cytotrace_violin_summary.csv")

print("\nDone.")
