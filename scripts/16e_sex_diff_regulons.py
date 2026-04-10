#!/usr/bin/env python3
"""
16e_sex_diff_regulons.py
Eosinophil Sex Differential Regulons heatmap (Figure 5C equivalent).

Identifies regulons enriched in WT Female vs WT Male eosinophils,
then shows z-scored AUCell means across WT♀, WT♂, KO♀, KO♂.

SLURM: --partition=epyc --mem=16G --cpus-per-task=2 --time=00:30:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"])

BASE      = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
AUC_FILE  = os.path.join(BASE, "scenic", "eos_scenic_aucell.csv")
META_FILE = os.path.join(BASE, "scenic", "eos_cell_metadata.csv")
FIG_DIR   = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

TOP_FEMALE  = 3      # top female-enriched regulons
TOP_MALE    = 2      # top male-enriched regulons
PVAL_THRESH = 0.05
ZSCORE_CAP  = 2.5

EXCLUDE = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

print("=" * 60)
print("16e_sex_diff_regulons.py")
print("=" * 60)

# =============================================================================
# Load data
# =============================================================================
print("\nLoading AUCell scores...")
auc = pd.read_csv(AUC_FILE, index_col=0)
print(f"  {auc.shape[0]} cells x {auc.shape[1]} regulons")

print("Loading metadata...")
meta = pd.read_csv(META_FILE)
# set barcode as index
if "barcode" in meta.columns:
    meta = meta.set_index("barcode")
print(f"  {meta.shape[0]} cells")
print(f"  Sex values: {meta['sex'].unique().tolist()}")
print(f"  Genotype values: {meta['genotype'].unique().tolist()}")

# Align
common = auc.index.intersection(meta.index)
auc  = auc.loc[common]
meta = meta.loc[common]
print(f"  Aligned: {len(common)} cells")

# Strip (+ ) suffix from regulon names -> TF names
def strip_suffix(col):
    return col.split("(")[0].strip()

auc.columns = [strip_suffix(c) for c in auc.columns]
# Remove duplicates (keep max)
auc = auc.T.groupby(level=0).max().T

# =============================================================================
# Group by sex x genotype
# =============================================================================
GROUPS = [
    ("Female", "WT"),
    ("Male",   "WT"),
    ("Female", "KO"),
    ("Male",   "KO"),
]
GROUP_LABELS = ["WT ♀", "WT ♂", "KO ♀", "KO ♂"]

group_cells = {}
for sex, geno in GROUPS:
    mask = (meta["sex"] == sex) & (meta["genotype"] == geno)
    group_cells[(sex, geno)] = meta.index[mask].tolist()
    print(f"  {sex} {geno}: {mask.sum()} cells")

# =============================================================================
# Statistical test: WT Female vs WT Male
# =============================================================================
print("\nTesting WT Female vs WT Male...")

wt_f = auc.loc[group_cells[("Female", "WT")]]
wt_m = auc.loc[group_cells[("Male",   "WT")]]

results = []
for tf in auc.columns:
    if tf in EXCLUDE:
        continue
    f_vals = wt_f[tf].values
    m_vals = wt_m[tf].values
    _, pval = stats.mannwhitneyu(f_vals, m_vals, alternative="two-sided")
    mean_f = f_vals.mean()
    mean_m = m_vals.mean()
    fc = (mean_f + 1e-6) / (mean_m + 1e-6)
    results.append({"tf": tf, "pval": pval, "mean_f": mean_f,
                    "mean_m": mean_m, "fc": fc})

res = pd.DataFrame(results)
res["log2fc"] = np.log2(res["fc"])

sig = res[res["pval"] < PVAL_THRESH].copy()
print(f"  Significant regulons: {len(sig)}")

female_top = (sig[sig["log2fc"] > 0]
              .nlargest(TOP_FEMALE, "log2fc")["tf"].tolist())
male_top   = (sig[sig["log2fc"] < 0]
              .nsmallest(TOP_MALE, "log2fc")["tf"].tolist())

selected = female_top + male_top
print(f"  Female-enriched: {female_top}")
print(f"  Male-enriched:   {male_top}")

if not selected:
    raise ValueError("No significant sex-differential regulons found. "
                     "Try lowering PVAL_THRESH.")

# =============================================================================
# Build heatmap matrix
# =============================================================================
group_means = {}
for (sex, geno), label in zip(GROUPS, GROUP_LABELS):
    cells = group_cells[(sex, geno)]
    group_means[label] = auc.loc[cells, selected].mean()

mat = pd.DataFrame(group_means).T   # groups x tfs

# Row-wise z-score across groups (per regulon)
row_mean = mat.mean(axis=0)
row_std  = mat.std(axis=0).replace(0, 1)
zmat = ((mat - row_mean) / row_std).clip(-ZSCORE_CAP, ZSCORE_CAP)

# =============================================================================
# Plot
# =============================================================================
n_rows = len(selected)
n_cols = len(GROUP_LABELS)
n_female = len(female_top)

cell_h    = 0.90
cell_w    = 0.80
right_pad = 1.5
top_pad   = 0.8
bot_pad   = 0.9
sidebar_w = 0.18   # sidebar bar width in inches
label_w   = 1.8    # space for y-axis label text
gap       = 0.10   # gap between sidebar right edge and heatmap

# Layout: [label_w][sidebar_w][gap][heatmap][right_pad]
left_pad = label_w + sidebar_w + gap

fig_h = n_rows * cell_h + top_pad + bot_pad
fig_w = left_pad + n_cols * cell_w + right_pad

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

ax_left   = left_pad / fig_w
ax_bottom = bot_pad / fig_h
ax_width  = (n_cols * cell_w) / fig_w
ax_height = (n_rows * cell_h) / fig_h
ax.set_position([ax_left, ax_bottom, ax_width, ax_height])

# Heatmap (tfs x groups)
plot_data = zmat.T.values   # tfs x groups

im = ax.imshow(plot_data, cmap=ZSCORE_CMAP, aspect="auto",
               vmin=-ZSCORE_CAP, vmax=ZSCORE_CAP)

# Axes labels
ax.set_xticks(range(n_cols))
ax.set_xticklabels(GROUP_LABELS, fontsize=11, fontweight="bold")
ax.set_yticks(range(n_rows))
ax.set_yticklabels(selected, fontsize=12, ha="right", fontweight="bold")
ax.yaxis.set_ticks_position("left")
ax.tick_params(axis="both", length=0)

# Vertical separator between WT and KO
ax.axvline(x=1.5, color="black", linewidth=2.0)

# White grid lines between rows
for i in range(1, n_rows):
    ax.axhline(y=i - 0.5, color="white", linewidth=0.8)
for j in range(1, n_cols):
    if j != 2:   # skip separator position (already drawn)
        ax.axvline(x=j - 0.5, color="white", linewidth=0.8)

# Color sidebar (female = purple, male = green)
FEMALE_COLOR = "#b084c8"
MALE_COLOR   = "#8fbc3b"

# Sidebar: sits just to the left of heatmap, after label text space
sidebar_ax_left = (label_w + 0.02) / fig_w
sidebar_ax = fig.add_axes([sidebar_ax_left, ax_bottom,
                            sidebar_w / fig_w, ax_height])
sidebar_colors = ([FEMALE_COLOR] * n_female +
                  [MALE_COLOR]   * len(male_top))
for i, color in enumerate(sidebar_colors):
    sidebar_ax.add_patch(plt.Rectangle((0, i - 0.5), 1, 1,
                                        color=color, linewidth=0))
sidebar_ax.set_xlim(0, 1)
sidebar_ax.set_ylim(-0.5, n_rows - 0.5)
sidebar_ax.invert_yaxis()   # match imshow top-to-bottom row order
sidebar_ax.axis("off")

# Title
ax.set_title("Eosinophil Sex Differential Regulons",
             fontsize=13, fontweight="bold", pad=12)

# Colorbar
cbar_ax = fig.add_axes([ax_left + ax_width + 0.02,
                         ax_bottom + ax_height * 0.2,
                         0.025,
                         ax_height * 0.6])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("Z-score", fontsize=9)
cbar.set_ticks([-2, -1, 0, 1, 2])
cbar.ax.tick_params(labelsize=8)

# Legend
import matplotlib.patches as mpatches
female_patch = mpatches.Patch(color=FEMALE_COLOR, label="Female-enriched in WT")
male_patch   = mpatches.Patch(color=MALE_COLOR,   label="Male-enriched in WT")
ax.legend(handles=[female_patch, male_patch], fontsize=8,
          loc="lower right", bbox_to_anchor=(1.0, -0.22),
          frameon=False)

out_pdf = os.path.join(FIG_DIR, "Fig5_sex_diff_regulons.pdf")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved: {out_pdf}")

res.to_csv(os.path.join(BASE, "scenic", "sex_diff_regulons_stats.csv"), index=False)
print("16e complete.")
