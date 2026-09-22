#!/usr/bin/env python3
"""
41b_hu_regulon_boxplot.py
Cross-dataset regulon comparison: boxplots + dendrogram
Qiu C0-C5 vs Hu eWAT/scWAT/BAT

SLURM: --partition=batch --mem=64G --cpus-per-task=4 --time=1:00:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
FIG_DIR    = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]
TOP_PER_CLUSTER = 3

EXCLUDE = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

# Colors
CLUSTER_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4",
}
HU_COLORS = {
    "Hu_eWAT": "#E63946", "Hu_scWAT": "#457B9D", "Hu_BAT": "#2A9D8F",
}
ALL_COLORS = {**CLUSTER_COLORS, **HU_COLORS}

GROUP_ORDER = ["C0", "C1", "C2", "C3", "C4", "C5", "Hu_eWAT", "Hu_scWAT", "Hu_BAT"]

print("=" * 60)
print("41b_hu_regulon_boxplot.py")
print("=" * 60)

# =============================================================================
# 1. Load our per-cluster AUCell data
# =============================================================================
print("\nStep 1: Loading our AUCell data...")

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

eos_data = {}
eos_meta = {}
for cl in EOS_CLUSTERS:
    auc_path = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(auc_path):
        df = pd.read_csv(auc_path, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T
        eos_data[cl] = df
        eos_meta[cl] = pd.Series(cl, index=df.index)
        print(f"  {cl}: {df.shape[0]} cells")

missing_clusters = [cl for cl in EOS_CLUSTERS if cl not in eos_data]
if missing_clusters:
    raise FileNotFoundError(
        f"Missing auc_matrix.csv for cluster(s) {missing_clusters}. "
        f"Regenerate the missing cluster's SCENIC output rather than "
        f"silently continuing with an incomplete analysis.")

# Merge all eos cells. SCENIC is inferred independently per cluster, so each
# cluster's regulon set differs; an outer join + fillna(0) would record "not
# inferred in this cluster's independently fitted network" as a real AUCell
# score of zero, indistinguishable from "scored and found inactive" -- which
# systematically makes regulons look cluster-specific. Inner join instead:
# keep only regulons scored in every cluster, so every remaining value is a
# genuine AUCell score.
all_eos_auc = pd.concat(eos_data.values(), join="inner")
all_eos_group = pd.concat(eos_meta.values())
print(f"  Total eos: {len(all_eos_auc)} cells, {all_eos_auc.shape[1]} regulons "
      f"scored in all {len(EOS_CLUSTERS)} clusters")

# Compute cluster means and select top regulons
all_tfs = sorted(all_eos_auc.columns)
cluster_means = {}
for cl in EOS_CLUSTERS:
    mask = all_eos_group == cl
    cluster_means[cl] = all_eos_auc.loc[mask].mean()
mean_df = pd.DataFrame(cluster_means)

selected_tfs = {}
for cl in EOS_CLUSTERS:
    scores = []
    for tf in all_tfs:
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

print(f"\nSelected {len(selected_tfs)} regulons:")
for cl in EOS_CLUSTERS:
    tfs = [tf for tf, c in selected_tfs.items() if c == cl]
    if tfs:
        print(f"  {cl}: {', '.join(tfs)}")

# =============================================================================
# 2. Load Hu regulon scores (computed by script 41)
# =============================================================================
print("\nStep 2: Loading Hu data...")

hu_expr_path = os.path.join(BASE, "scenic", "hu_expr_matrix.csv")
hu_meta_path = os.path.join(BASE, "scenic", "hu_depot_metadata.csv")

if not os.path.exists(hu_expr_path):
    print("  ERROR: Run script 41 first to export Hu expression matrix")
    exit(1)

hu_expr = pd.read_csv(hu_expr_path, index_col=0)
hu_meta = pd.read_csv(hu_meta_path).set_index("barcode")
print(f"  Hu: {hu_expr.shape[1]} cells x {hu_expr.shape[0]} genes")

# Load regulon gene sets
print("  Loading regulon gene sets...")
regulon_genes = {}
for cl in EOS_CLUSTERS:
    reg_path = os.path.join(SCENIC_DIR, cl, "regulons.csv")
    if not os.path.exists(reg_path):
        continue
    try:
        from pyscenic.utils import load_motifs
        from pyscenic.prune import df2regulons
        motifs = load_motifs(reg_path)
        regs = df2regulons(motifs)
        for r in regs:
            tf = r.name.split("(")[0].strip()
            if tf in selected_tfs:
                genes = list(r.genes)
                if tf not in regulon_genes or len(genes) > len(regulon_genes[tf]):
                    regulon_genes[tf] = genes
    except Exception as e:
        pass

# Score Hu cells
print("  Scoring Hu cells...")
hu_scores = {}
for tf, genes in regulon_genes.items():
    shared = [g for g in genes if g in hu_expr.index]
    if len(shared) >= 3:
        hu_scores[tf] = hu_expr.loc[shared].mean(axis=0)

hu_scores_df = pd.DataFrame(hu_scores)
hu_group = hu_meta["depot"].map(lambda x: f"Hu_{x}")
print(f"  Scored {len(hu_scores)} regulons")

# =============================================================================
# 3. Build combined per-cell data for boxplots
# =============================================================================
print("\nStep 3: Building combined data...")

# For our data: use AUCell scores directly
# For Hu data: use mean-expression scores (different scale, so we'll z-score)
shared_regulons = sorted(set(selected_tfs.keys()) & set(hu_scores_df.columns))
print(f"  Shared regulons for plotting: {len(shared_regulons)}")

# Subsample for plotting (too many cells makes boxplots slow)
np.random.seed(42)
MAX_CELLS = 500

plot_rows = []
for cl in EOS_CLUSTERS:
    mask = all_eos_group == cl
    cells = all_eos_auc.index[mask]
    if len(cells) > MAX_CELLS:
        cells = np.random.choice(cells, MAX_CELLS, replace=False)
    for cell in cells:
        for tf in shared_regulons:
            plot_rows.append({
                "group": cl,
                "regulon": tf,
                "score": all_eos_auc.loc[cell, tf],
                "dataset": "Qiu"
            })

for depot in ["eWAT", "scWAT", "BAT"]:
    mask = hu_meta["depot"] == depot
    cells = hu_meta.index[mask].intersection(hu_scores_df.index)
    if len(cells) > MAX_CELLS:
        cells = np.random.choice(cells, MAX_CELLS, replace=False)
    for cell in cells:
        for tf in shared_regulons:
            plot_rows.append({
                "group": f"Hu_{depot}",
                "regulon": tf,
                "score": hu_scores_df.loc[cell, tf],
                "dataset": "Hu"
            })

plot_df = pd.DataFrame(plot_rows)
print(f"  Plot data: {len(plot_df)} rows")

# =============================================================================
# 4. Boxplots — one per regulon, faceted
# =============================================================================
print("\nStep 4: Generating boxplots...")

# Order regulons by cluster
regulon_order = []
for cl in EOS_CLUSTERS:
    tfs = [tf for tf in shared_regulons if selected_tfs.get(tf) == cl]
    regulon_order.extend(sorted(tfs))

n_reg = len(regulon_order)
n_cols = 3
n_rows_plot = int(np.ceil(n_reg / n_cols))

fig, axes = plt.subplots(n_rows_plot, n_cols, figsize=(4.5 * n_cols, 3.2 * n_rows_plot))
axes = axes.flatten()

for idx, tf in enumerate(regulon_order):
    ax = axes[idx]
    sub = plot_df[plot_df["regulon"] == tf]

    # Only keep groups that exist
    groups_present = [g for g in GROUP_ORDER if g in sub["group"].values]
    palette = [ALL_COLORS[g] for g in groups_present]

    sns.boxplot(data=sub, x="group", y="score", order=groups_present,
                palette=palette, fliersize=0, linewidth=0.8,
                boxprops=dict(alpha=0.7), ax=ax)

    # Add separator line between datasets
    n_eos_groups = len([g for g in groups_present if not g.startswith("Hu")])
    if n_eos_groups < len(groups_present):
        ax.axvline(x=n_eos_groups - 0.5, color="black", linewidth=1.5, linestyle="--")

    cl = selected_tfs.get(tf, "")
    ax.set_title(f"{tf} ({cl})", fontsize=10, fontweight="bold", fontstyle="italic")
    ax.set_xlabel("")
    ax.set_ylabel("Score" if idx % n_cols == 0 else "")
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.tick_params(axis="y", labelsize=8)

# Hide empty axes
for idx in range(n_reg, len(axes)):
    axes[idx].axis("off")

fig.suptitle("Regulon Activity: Qiu Clusters vs Hu Adipose Depots",
             fontsize=14, fontweight="bold", y=1.01)
fig.tight_layout()

out_box = os.path.join(FIG_DIR, "41b_regulon_boxplots_qiu_hu.pdf")
fig.savefig(out_box, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_box}")

# =============================================================================
# 5. Dendrogram based on regulon activity profiles
# =============================================================================
print("\nStep 5: Generating dendrogram...")

# Compute group means for shared regulons
group_means = {}
for cl in EOS_CLUSTERS:
    mask = all_eos_group == cl
    vals = {}
    for tf in shared_regulons:
        # Safe without a fallback: shared_regulons is a subset of
        # selected_tfs, itself drawn only from all_eos_auc.columns (the
        # inner-joined, all-clusters-scored set above).
        vals[tf] = all_eos_auc.loc[mask, tf].mean()
    group_means[cl] = vals

for depot in ["eWAT", "scWAT", "BAT"]:
    mask = hu_meta["depot"] == depot
    cells = hu_meta.index[mask].intersection(hu_scores_df.index)
    vals = {}
    for tf in shared_regulons:
        vals[tf] = hu_scores_df.loc[cells, tf].mean()
    group_means[f"Hu_{depot}"] = vals

mean_mat = pd.DataFrame(group_means).T  # groups x regulons

# Z-score per regulon before computing distances
mean_z = (mean_mat - mean_mat.mean()) / mean_mat.std().replace(0, 1)

# Hierarchical clustering
dist = pdist(mean_z.values, metric="correlation")
Z = linkage(dist, method="average")

# Label colors
label_colors = {}
for g in mean_z.index:
    if g.startswith("Hu"):
        label_colors[g] = "#4DAF4A"
    else:
        label_colors[g] = "#E63946"

fig, ax = plt.subplots(figsize=(10, 6))
dend = dendrogram(Z, labels=mean_z.index.tolist(), ax=ax,
                  leaf_font_size=12, leaf_rotation=0,
                  above_threshold_color="grey")

# Color the labels
xlbls = ax.get_xticklabels()
for lbl in xlbls:
    txt = lbl.get_text()
    lbl.set_color(label_colors.get(txt, "black"))
    lbl.set_fontweight("bold")

ax.set_ylabel("1 - Correlation", fontsize=12)
ax.set_title("Regulon-Based Cross-Dataset Similarity", fontsize=14, fontweight="bold")

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#E63946',
           markersize=10, label='Qiu (HFD Adipose)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#4DAF4A',
           markersize=10, label='Hu (Healthy Adipose)')
]
ax.legend(handles=legend_elements, loc="upper right", fontsize=10, frameon=False)

out_dend = os.path.join(FIG_DIR, "41b_regulon_dendrogram_qiu_hu.pdf")
fig.savefig(out_dend, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_dend}")

# Print dendrogram groupings
print("\n  Dendrogram leaf order:", [mean_z.index[i] for i in dend["leaves"]])

print("\n41b complete.")
