#!/usr/bin/env python3
"""
42b_regulon_compare_aucell.py
Cross-dataset regulon comparison using AUCell from BOTH datasets.
Apples-to-apples: both Qiu and Hu scored with pySCENIC AUCell.

Generates:
  - Boxplots of shared regulons across C0-C5 + eWAT/scWAT/BAT
  - Dendrogram based on AUCell regulon profiles
  - Combined heatmap

SLURM: --partition=batch --mem=32G --cpus-per-task=4 --time=00:30:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
HU_SCENIC  = os.path.join(BASE, "scenic", "hu_scenic")
HU_META    = os.path.join(BASE, "scenic", "hu_depot_metadata.csv")
FIG_DIR    = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]
TOP_PER_CLUSTER = 3
ZSCORE_CAP = 2.5

ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"])

EXCLUDE = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

CLUSTER_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4",
}
HU_COLORS = {"Hu_eWAT": "#E63946", "Hu_scWAT": "#457B9D", "Hu_BAT": "#2A9D8F"}
ALL_COLORS = {**CLUSTER_COLORS, **HU_COLORS}
GROUP_ORDER = ["C0", "C1", "C2", "C3", "C4", "C5", "Hu_eWAT", "Hu_scWAT", "Hu_BAT"]

print("=" * 60)
print("42b_regulon_compare_aucell.py")
print("  Both datasets scored with pySCENIC AUCell")
print("=" * 60)

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

# =============================================================================
# 1. Load our per-cluster AUCell
# =============================================================================
print("\nStep 1: Loading Qiu AUCell data...")

eos_auc_parts = []
eos_group_parts = []
for cl in EOS_CLUSTERS:
    auc_path = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(auc_path):
        df = pd.read_csv(auc_path, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T
        eos_auc_parts.append(df)
        eos_group_parts.append(pd.Series(cl, index=df.index))
        print(f"  {cl}: {df.shape[0]} cells x {df.shape[1]} regulons")
    else:
        raise FileNotFoundError(
            f"Missing auc_matrix.csv for cluster {cl}. Regenerate the "
            f"missing cluster's SCENIC output rather than silently "
            f"continuing with an incomplete analysis.")

# SCENIC is inferred independently per cluster, so each cluster's regulon set
# differs; an outer join + fillna(0) would record "not inferred in this
# cluster's independently fitted network" as a real AUCell score of zero,
# indistinguishable from "scored and found inactive" -- which systematically
# makes regulons look cluster-specific. Inner join instead: keep only
# regulons scored in every cluster.
qiu_auc = pd.concat(eos_auc_parts, join="inner")
qiu_group = pd.concat(eos_group_parts)
print(f"  Total Qiu: {len(qiu_auc)} cells, {qiu_auc.shape[1]} regulons "
      f"scored in all {len(EOS_CLUSTERS)} clusters")

# =============================================================================
# 2. Load Hu AUCell (combined SCENIC, split by depot)
# =============================================================================
print("\nStep 2: Loading Hu AUCell data...")

# Try combined first, then per-depot fallback
combined_auc_path = os.path.join(HU_SCENIC, "combined", "auc_matrix.csv")
if os.path.exists(combined_auc_path):
    print("  Loading combined Hu AUCell...")
    hu_auc = pd.read_csv(combined_auc_path, index_col=0)
    hu_auc.columns = [extract_tf(c) for c in hu_auc.columns]
    hu_auc = hu_auc.T.groupby(level=0).max().T
    print(f"  Hu: {hu_auc.shape[0]} cells x {hu_auc.shape[1]} regulons")

    # Assign depot labels from metadata
    hu_meta = pd.read_csv(HU_META).set_index("barcode")
    common_hu = hu_auc.index.intersection(hu_meta.index)
    hu_auc = hu_auc.loc[common_hu]
    hu_group = hu_meta.loc[common_hu, "depot"].map(lambda x: f"Hu_{x}")
else:
    # Fallback: load per-depot AUCell
    print("  Combined not found, loading per-depot...")
    hu_auc_parts = []
    hu_group_parts = []
    for depot in ["eWAT", "scWAT", "BAT"]:
        auc_path = os.path.join(HU_SCENIC, depot, "auc_matrix.csv")
        if not os.path.exists(auc_path):
            print(f"  WARNING: {auc_path} not found, skipping {depot}")
            continue
        df = pd.read_csv(auc_path, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T
        hu_auc_parts.append(df)
        hu_group_parts.append(pd.Series(f"Hu_{depot}", index=df.index))
        print(f"  Hu_{depot}: {df.shape[0]} cells x {df.shape[1]} regulons")

    if not hu_auc_parts:
        print("  ERROR: No Hu AUCell data found.")
        exit(1)
    # Same reasoning as the Qiu per-cluster merge above: inner join, not
    # outer + fillna(0), so a regulon not inferred for a given depot's
    # independently fitted network isn't recorded as a real zero score.
    hu_auc = pd.concat(hu_auc_parts, join="inner")
    hu_group = pd.concat(hu_group_parts)

print(f"  Total Hu: {len(hu_auc)} cells, {hu_auc.shape[1]} regulons")
print(f"  Depots: {hu_group.value_counts().to_dict()}")

# =============================================================================
# 3. Find shared regulons and select top per cluster
# =============================================================================
print("\nStep 3: Finding shared regulons...")

shared_all = sorted(set(qiu_auc.columns) & set(hu_auc.columns))
print(f"  Shared regulons between datasets: {len(shared_all)}")

# Compute Qiu cluster means on shared regulons
cluster_means = {}
for cl in EOS_CLUSTERS:
    mask = qiu_group == cl
    cluster_means[cl] = qiu_auc.loc[mask, shared_all].mean()
mean_df = pd.DataFrame(cluster_means)

# Select top per cluster
selected_tfs = {}
for cl in EOS_CLUSTERS:
    scores = []
    for tf in shared_all:
        if tf in EXCLUDE:
            continue
        max_val = mean_df.loc[tf, cl]
        if max_val < 0.02:
            continue
        if cl != mean_df.columns[np.argmax(mean_df.loc[tf].values)]:
            continue
        sorted_vals = sorted(mean_df.loc[tf].values, reverse=True)
        uniqueness = sorted_vals[0] - sorted_vals[1] if len(sorted_vals) >= 2 else sorted_vals[0]
        scores.append((tf, max_val, uniqueness))
    scores.sort(key=lambda x: x[2], reverse=True)
    for tf, _, _ in scores[:TOP_PER_CLUSTER]:
        selected_tfs[tf] = cl

print(f"\nSelected {len(selected_tfs)} shared regulons:")
for cl in EOS_CLUSTERS:
    tfs = [tf for tf, c in selected_tfs.items() if c == cl]
    if tfs:
        print(f"  {cl}: {', '.join(tfs)}")

sel_list = list(selected_tfs.keys())

# =============================================================================
# 4. Boxplots
# =============================================================================
print("\nStep 4: Generating boxplots...")

np.random.seed(42)
MAX_CELLS = 500

plot_rows = []
for cl in EOS_CLUSTERS:
    mask = qiu_group == cl
    cells = qiu_auc.index[mask]
    if len(cells) > MAX_CELLS:
        cells = np.random.choice(cells, MAX_CELLS, replace=False)
    for cell in cells:
        for tf in sel_list:
            plot_rows.append({"group": cl, "regulon": tf,
                              "AUCell": qiu_auc.loc[cell, tf]})

for depot_label in ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]:
    mask = hu_group == depot_label
    cells = hu_auc.index[mask]
    if len(cells) > MAX_CELLS:
        cells = np.random.choice(cells, MAX_CELLS, replace=False)
    for cell in cells:
        for tf in sel_list:
            plot_rows.append({"group": depot_label, "regulon": tf,
                              "AUCell": hu_auc.loc[cell, tf]})

plot_df = pd.DataFrame(plot_rows)
print(f"  Plot data: {len(plot_df)} rows")

# Order regulons by cluster
regulon_order = []
for cl in EOS_CLUSTERS:
    regulon_order.extend([tf for tf in sel_list if selected_tfs[tf] == cl])

n_reg = len(regulon_order)
n_cols = 3
n_rows_plot = int(np.ceil(n_reg / n_cols))

fig, axes = plt.subplots(n_rows_plot, n_cols, figsize=(4.5 * n_cols, 3.2 * n_rows_plot))
axes = axes.flatten()

for idx, tf in enumerate(regulon_order):
    ax = axes[idx]
    sub = plot_df[plot_df["regulon"] == tf]
    groups_present = [g for g in GROUP_ORDER if g in sub["group"].values]
    palette = [ALL_COLORS[g] for g in groups_present]

    sns.boxplot(data=sub, x="group", y="AUCell", order=groups_present,
                palette=palette, fliersize=0, linewidth=0.8,
                boxprops=dict(alpha=0.7), ax=ax)

    n_eos = len([g for g in groups_present if not g.startswith("Hu")])
    if n_eos < len(groups_present):
        ax.axvline(x=n_eos - 0.5, color="black", linewidth=1.5, linestyle="--")

    cl = selected_tfs[tf]
    ax.set_title(f"{tf} ({cl})", fontsize=10, fontweight="bold", fontstyle="italic")
    ax.set_xlabel("")
    ax.set_ylabel("AUCell" if idx % n_cols == 0 else "")
    ax.tick_params(axis="x", rotation=45, labelsize=7)

for idx in range(n_reg, len(axes)):
    axes[idx].axis("off")

fig.suptitle("Regulon Activity (AUCell): Qiu Clusters vs Hu Depots",
             fontsize=14, fontweight="bold", y=1.01)
fig.tight_layout()
out_box = os.path.join(FIG_DIR, "42b_aucell_boxplots_qiu_hu.pdf")
fig.savefig(out_box, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_box}")

# =============================================================================
# 5. Dendrogram
# =============================================================================
print("\nStep 5: Generating dendrogram...")

group_means = {}
for cl in EOS_CLUSTERS:
    mask = qiu_group == cl
    group_means[cl] = qiu_auc.loc[mask, sel_list].mean()

for depot_label in ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]:
    mask = hu_group == depot_label
    group_means[depot_label] = hu_auc.loc[mask, sel_list].mean()

mean_mat = pd.DataFrame(group_means).T
mean_z = (mean_mat - mean_mat.mean()) / mean_mat.std().replace(0, 1)

dist = pdist(mean_z.values, metric="correlation")
Z = linkage(dist, method="average")

fig, ax = plt.subplots(figsize=(10, 6))
dend = dendrogram(Z, labels=mean_z.index.tolist(), ax=ax,
                  leaf_font_size=12, leaf_rotation=0,
                  above_threshold_color="grey")

for lbl in ax.get_xticklabels():
    txt = lbl.get_text()
    lbl.set_color("#4DAF4A" if txt.startswith("Hu") else "#E63946")
    lbl.set_fontweight("bold")

ax.set_ylabel("1 - Correlation", fontsize=12)
ax.set_title("Regulon-Based Similarity (AUCell, both datasets)",
             fontsize=14, fontweight="bold")

from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#E63946',
           markersize=10, label='Qiu (HFD Adipose)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#4DAF4A',
           markersize=10, label='Hu (Healthy Adipose)')
], loc="upper right", fontsize=10, frameon=False)

out_dend = os.path.join(FIG_DIR, "42b_aucell_dendrogram_qiu_hu.pdf")
fig.savefig(out_dend, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_dend}")
print("  Leaf order:", [mean_z.index[i] for i in dend["leaves"]])

# =============================================================================
# 6. Combined heatmap
# =============================================================================
print("\nStep 6: Combined heatmap...")

combined = pd.DataFrame(group_means).T[regulon_order].T
row_mean = combined.mean(axis=1)
row_std = combined.std(axis=1).replace(0, 1)
zmat = ((combined.T - row_mean) / row_std).T.clip(-ZSCORE_CAP, ZSCORE_CAP)

col_order = [c for c in GROUP_ORDER if c in zmat.columns]
zmat = zmat[col_order]

n_r = len(regulon_order)
n_c = len(col_order)
n_eos = len([c for c in EOS_CLUSTERS if c in col_order])

fig, ax = plt.subplots(figsize=(n_c * 0.65 + 3.5, n_r * 0.4 + 2.5))
im = ax.imshow(zmat.values, cmap=ZSCORE_CMAP, aspect="auto",
               vmin=-ZSCORE_CAP, vmax=ZSCORE_CAP)

ax.set_xticks(range(n_c))
ax.set_xticklabels(col_order, fontsize=10, fontweight="bold", rotation=45, ha="right")
ax.set_yticks(range(n_r))
ax.set_yticklabels(regulon_order, fontsize=9, fontstyle="italic")
ax.tick_params(axis="both", length=0)

if n_eos < n_c:
    ax.axvline(x=n_eos - 0.5, color="black", linewidth=2.0)

# Cluster sidebar
y = 0
for cl in EOS_CLUSTERS:
    tfs = [tf for tf in regulon_order if selected_tfs.get(tf) == cl]
    for i, tf in enumerate(tfs):
        ax.add_patch(plt.Rectangle((-0.8, y + i - 0.5), 0.3, 1,
                                    color=CLUSTER_COLORS[cl], clip_on=False))
    if tfs:
        y += len(tfs)
        if y < n_r:
            ax.axhline(y=y - 0.5, color="white", linewidth=1.5)

ax.text(n_eos / 2 - 0.5, -1.5, "Qiu (HFD)",
        ha="center", fontsize=11, fontweight="bold", color="#E63946")
ax.text(n_eos + (n_c - n_eos) / 2 - 0.5, -1.5, "Hu (Healthy)",
        ha="center", fontsize=11, fontweight="bold", color="#4DAF4A")
ax.set_title("Cross-Dataset Regulon Comparison (AUCell)",
             fontsize=13, fontweight="bold", pad=25)

cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.04)
cbar.set_label("Z-score", fontsize=9)

out_heat = os.path.join(FIG_DIR, "42b_aucell_heatmap_qiu_hu.pdf")
fig.savefig(out_heat, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_heat}")

# =============================================================================
# 7. Summary
# =============================================================================
print("\n========== Summary ==========")
print(f"Qiu regulons: {qiu_auc.shape[1]}")
print(f"Hu regulons: {hu_auc.shape[1]}")
print(f"Shared: {len(shared_all)}")
print(f"Selected for comparison: {len(selected_tfs)}")

print("\n42b complete.")
