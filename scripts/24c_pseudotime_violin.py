#!/usr/bin/env python
"""
24c_pseudotime_violin.py
Fig 2H-style violin plot: Diffusion pseudotime by group, ranked by median DPT.
Uses TWO-WAY integration (Qiu + Gurtner) from script 13b to match original
manuscript Methods (CellRank2 + Scanpy DPT, root = Gurtner bone marrow).

Groups:
  - Gurtner Bone Marrow (most immature, lowest DPT)
  - Gurtner Blood
  - Qiu eos clusters C0-C6

Input:  13b CellRank scores + metadata from two-way integration
Output: violin plot PDF
SLURM:  batch, 16G, 2 CPUs, 30min, cellrank_env
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
}

# --- Paths (two-way integration from script 13b) ----------------------------
BASE = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
INTEG_DIR = os.path.join(BASE, "integration")
META_PATH = os.path.join(INTEG_DIR, "for_cellrank", "metadata.csv")
SCORES_PATH = os.path.join(INTEG_DIR, "cellrank_cytotrace_scores.csv")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# --- Load data ---------------------------------------------------------------
print("Loading two-way integration data (script 13b) ...")

meta = pd.read_csv(META_PATH, index_col=0)
scores = pd.read_csv(SCORES_PATH, index_col=0)

print(f"  Metadata: {len(meta)} cells")
print(f"  Scores:   {len(scores)} cells")
print(f"  Score columns: {scores.columns.tolist()}")

# Merge — use dpt_pseudotime as primary metric (matches original manuscript)
if "dpt_pseudotime" not in scores.columns:
    print("ERROR: dpt_pseudotime column not found in scores. "
          "Make sure script 13b completed successfully with DPT.")
    raise SystemExit(1)

df = meta.join(scores[["dpt_pseudotime", "ct_score"]], how="inner")
print(f"  Merged: {len(df)} cells")

# Remove cells with infinite DPT
n_inf = np.isinf(df["dpt_pseudotime"]).sum()
if n_inf > 0:
    print(f"  Removing {n_inf} cells with infinite DPT")
    df = df[~np.isinf(df["dpt_pseudotime"])]

# --- Assign violin groups (two-way: Qiu + Gurtner only) ---------------------
print("Assigning violin groups ...")

def assign_group(row):
    ds = row.get("dataset", "")
    organ = row.get("organ", "")
    cluster = row.get("original_cluster", "")

    if ds == "Gurtner_2022_Blood_BM" or ds == "Gurtner_2022":
        if organ == "Bone Marrow":
            return "Gurtner\nBone Marrow"
        elif organ == "Blood":
            return "Gurtner\nBlood"
        else:
            return None
    elif ds == "Qiu_2025_Adipose":
        return cluster  # C0, C1, ..., C6
    else:
        return None

df["violin_group"] = df.apply(assign_group, axis=1)
df = df.dropna(subset=["violin_group", "dpt_pseudotime"])

print(f"  Cells with valid groups: {len(df)}")
print(f"  Groups: {df['violin_group'].value_counts().to_dict()}")

# --- Order by median DPT (ascending: most immature / lowest DPT first) ------
group_medians = df.groupby("violin_group")["dpt_pseudotime"].median().sort_values(ascending=True)
ordered_groups = group_medians.index.tolist()
print(f"\n  Group order by median DPT (low = progenitor):")
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
ax.set_ylabel("Diffusion Pseudotime", fontsize=13, fontweight="bold")
ax.set_title("Diffusion Pseudotime by Group\n(root = Gurtner bone marrow)",
             fontsize=14, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

out_path = os.path.join(FIG_DIR, "24c_dpt_violin.pdf")
fig.savefig(out_path)
plt.close()
print(f"\n  Saved: {out_path}")

# Also save group summary
summary = df.groupby("violin_group")["dpt_pseudotime"].agg(
    ["count", "mean", "median", "std", "min", "max"]
).reindex(ordered_groups)
summary.to_csv(os.path.join(INTEG_DIR, "dpt_violin_summary.csv"))
print(f"  Saved: dpt_violin_summary.csv")

# --- Also generate CytoTRACE version for supplementary ----------------------
print("\nGenerating supplementary CytoTRACE violin ...")

ct_medians = df.groupby("violin_group")["ct_score"].median().sort_values(ascending=False)
ct_order = ct_medians.index.tolist()
ct_colors = [all_colors.get(g, "#999999") for g in ct_order]

fig2, ax2 = plt.subplots(figsize=(max(10, len(ct_order) * 0.9 + 2), 6))

ct_data = [df.loc[df["violin_group"] == g, "ct_score"].values for g in ct_order]

parts2 = ax2.violinplot(ct_data, positions=range(len(ct_order)),
                        showmeans=False, showmedians=False, showextrema=False)

for i, pc in enumerate(parts2["bodies"]):
    pc.set_facecolor(ct_colors[i])
    pc.set_edgecolor("black")
    pc.set_linewidth(0.8)
    pc.set_alpha(0.85)

ax2.boxplot(ct_data, positions=range(len(ct_order)),
            widths=0.15, showfliers=False, patch_artist=True,
            medianprops=dict(color="white", linewidth=2),
            boxprops=dict(facecolor="black", linewidth=0.8),
            whiskerprops=dict(color="black", linewidth=0.8),
            capprops=dict(color="black", linewidth=0.8))

ax2.set_xticks(range(len(ct_order)))
ax2.set_xticklabels(ct_order, rotation=45, ha="right", fontsize=10,
                    fontweight="bold")
ax2.set_ylabel("CytoTRACE Score", fontsize=13, fontweight="bold")
ax2.set_title("CytoTRACE Score by Group\n(high = progenitor-like)",
              fontsize=14, fontweight="bold")
ax2.set_ylim(-0.02, 1.02)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.tight_layout()
out_path2 = os.path.join(FIG_DIR, "24c_cytotrace_violin.pdf")
fig2.savefig(out_path2)
plt.close()
print(f"  Saved: {out_path2}")

print("\nDone.")
