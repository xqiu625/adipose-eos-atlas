#!/usr/bin/env python3
"""
41_hu_regulon_comparison.py
Cross-dataset regulon comparison: Qiu C0-C5 vs Hu eWAT/scWAT/BAT

Scores Hu eosinophils with our per-cluster SCENIC regulons (AUCell),
then compares regulon activity patterns across datasets.

Output:
  - Combined heatmap: regulons x (C0-C5 + eWAT/scWAT/BAT)
  - Spearman correlation: which Hu depot is most similar to which cluster
  - Per-depot regulon activity violin plots for top regulons

SLURM: --partition=batch --mem=64G --cpus-per-task=8 --time=2:00:00
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

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
HU_PATH    = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_adipose_eos.rds"
FIG_DIR    = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]
TOP_PER_CLUSTER = 5
ZSCORE_CAP = 2.5

ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"])

EXCLUDE = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

print("=" * 60)
print("41_hu_regulon_comparison.py")
print("=" * 60)

# =============================================================================
# 1. Load our per-cluster AUCell scores and select top regulons
# =============================================================================
print("\nStep 1: Loading our per-cluster AUCell data...")

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

eos_data = {}
for cl in EOS_CLUSTERS:
    auc_path = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(auc_path):
        df = pd.read_csv(auc_path, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T
        eos_data[cl] = df
        print(f"  {cl}: {df.shape[0]} cells x {df.shape[1]} regulons")

missing_clusters = [cl for cl in EOS_CLUSTERS if cl not in eos_data]
if missing_clusters:
    raise FileNotFoundError(
        f"Missing auc_matrix.csv for cluster(s) {missing_clusters}. "
        f"Regenerate the missing cluster's SCENIC output rather than "
        f"silently continuing with an incomplete analysis.")

# TF set = INTERSECTION (not union) of each cluster's independently inferred
# TFs. A TF absent from one cluster's network was "not inferred" there, not
# scored and found at AUCell 0 -- filling it as 0.0 would corrupt which TF
# looks "top" for that cluster and make regulons look artificially
# cluster-specific.
all_tfs = None
for df in eos_data.values():
    all_tfs = set(df.columns) if all_tfs is None else (all_tfs & set(df.columns))
all_tfs = sorted(all_tfs) if all_tfs else []
print(f"  TFs scored in all {len(eos_data)} loaded clusters: {len(all_tfs)}")

cluster_means = {}
for cl in EOS_CLUSTERS:
    if cl not in eos_data:
        continue
    df = eos_data[cl]
    cluster_means[cl] = {tf: df[tf].mean() for tf in all_tfs}

mean_df = pd.DataFrame(cluster_means)  # tfs x clusters

# Select top regulons per cluster by uniqueness
selected_tfs = {}
for cl in EOS_CLUSTERS:
    if cl not in mean_df.columns:
        continue
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
    for tf, auc_val, uniq in scores[:TOP_PER_CLUSTER]:
        selected_tfs[tf] = cl

print(f"\nSelected {len(selected_tfs)} regulons:")
for cl in EOS_CLUSTERS:
    tfs = [tf for tf, c in selected_tfs.items() if c == cl]
    if tfs:
        print(f"  {cl}: {', '.join(tfs)}")

# =============================================================================
# 2. Load regulon gene sets from per-cluster SCENIC
# =============================================================================
print("\nStep 2: Loading regulon gene sets...")

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
        print(f"  WARNING: Could not parse {reg_path}: {e}")

print(f"  Loaded gene sets for {len(regulon_genes)} regulons")
for tf in sorted(regulon_genes.keys()):
    print(f"    {tf}: {len(regulon_genes[tf])} target genes")

# =============================================================================
# 3. Load Hu eosinophil data and score with our regulons
# =============================================================================
print("\nStep 3: Loading Hu eosinophil data...")

# Load via R bridge (rds file)
import subprocess, tempfile

score_script = tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False)
score_script.write(f"""
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))
suppressPackageStartupMessages(library(Seurat))

hu <- readRDS("{HU_PATH}")
if (inherits(hu[["RNA"]], "Assay5")) hu[["RNA"]] <- JoinLayers(hu[["RNA"]])
hu <- NormalizeData(hu, verbose = FALSE)

# Find depot column
depot_col <- NULL
for (col in colnames(hu@meta.data)) {{
  vals <- as.character(unique(hu@meta.data[[col]]))
  if (any(c("BAT", "eWAT", "scWAT") %in% vals)) {{
    depot_col <- col
    break
  }}
}}

# Export normalized expression matrix
expr <- as.matrix(GetAssayData(hu, layer = "data"))
write.csv(expr, "{BASE}/scenic/hu_expr_matrix.csv")

# Export depot metadata
meta <- data.frame(
  barcode = colnames(hu),
  depot = if (!is.null(depot_col)) as.character(hu@meta.data[[depot_col]]) else "unknown"
)
write.csv(meta, "{BASE}/scenic/hu_depot_metadata.csv", row.names = FALSE)

cat(sprintf("Exported: %d cells x %d genes\\n", ncol(hu), nrow(hu)))
cat(sprintf("Depot column: %s\\n", depot_col))
cat("Depot distribution:\\n")
print(table(meta$depot))
""")
score_script.close()

# Check if Hu expression matrix already exported
hu_expr_path = os.path.join(BASE, "scenic", "hu_expr_matrix.csv")
hu_meta_path = os.path.join(BASE, "scenic", "hu_depot_metadata.csv")

if not os.path.exists(hu_expr_path):
    print("  Exporting Hu expression matrix via R...")
    result = subprocess.run(
        ["/bigdata/godziklab/shared/Xinru/cellrank_env/bin/Rscript", score_script.name],
        capture_output=True, text=True,
        env={**os.environ, "LD_LIBRARY_PATH": "/bigdata/godziklab/shared/Xinru/cellrank_env/lib:/opt/linux/rocky/8.x/x86_64/pkgs/gcc/12.2.0/lib64:" + os.environ.get("LD_LIBRARY_PATH", "")}
    )
    print(result.stdout)
    if result.returncode != 0:
        print(f"  R ERROR: {result.stderr}")
else:
    print("  Hu expression matrix already exported, loading...")

os.unlink(score_script.name)

# Load Hu expression
print("  Loading Hu expression matrix...")
hu_expr = pd.read_csv(hu_expr_path, index_col=0)
hu_meta = pd.read_csv(hu_meta_path)
hu_meta = hu_meta.set_index("barcode")
print(f"  Hu: {hu_expr.shape[1]} cells x {hu_expr.shape[0]} genes")
print(f"  Depots: {hu_meta['depot'].value_counts().to_dict()}")

# =============================================================================
# 4. Score Hu cells with our regulons using simple mean expression
# =============================================================================
print("\nStep 4: Scoring Hu cells with our regulons...")

# Use mean expression of regulon target genes as a simple AUC proxy
hu_regulon_scores = {}
for tf, genes in regulon_genes.items():
    # Find genes present in Hu data
    shared_genes = [g for g in genes if g in hu_expr.index]
    if len(shared_genes) < 3:
        print(f"  {tf}: only {len(shared_genes)} shared genes, skipping")
        continue
    # Mean expression across target genes per cell
    hu_regulon_scores[tf] = hu_expr.loc[shared_genes].mean(axis=0)
    print(f"  {tf}: {len(shared_genes)}/{len(genes)} genes found")

hu_scores_df = pd.DataFrame(hu_regulon_scores)
print(f"  Scored {hu_scores_df.shape[1]} regulons across {hu_scores_df.shape[0]} Hu cells")

# =============================================================================
# 5. Compute per-depot means and combine with our cluster means
# =============================================================================
print("\nStep 5: Computing per-depot means...")

# Hu depot means
hu_depot_means = {}
for depot in ["eWAT", "scWAT", "BAT"]:
    mask = hu_meta["depot"] == depot
    if mask.sum() == 0:
        continue
    cells = hu_meta.index[mask]
    common_cells = cells.intersection(hu_scores_df.index)
    if len(common_cells) == 0:
        continue
    hu_depot_means[f"Hu_{depot}"] = hu_scores_df.loc[common_cells].mean()
    print(f"  Hu_{depot}: {len(common_cells)} cells")

hu_mean_df = pd.DataFrame(hu_depot_means)

# Our cluster means (for the selected regulons only)
our_mean_df = mean_df.loc[mean_df.index.isin(hu_scores_df.columns),
                           [c for c in EOS_CLUSTERS if c in mean_df.columns]]

# Align regulons
shared_regulons = sorted(set(our_mean_df.index) & set(hu_mean_df.index))
print(f"\n  Shared regulons: {len(shared_regulons)}")

our_aligned = our_mean_df.loc[shared_regulons]
hu_aligned = hu_mean_df.loc[shared_regulons]

# Combine
combined = pd.concat([our_aligned, hu_aligned], axis=1)
print(f"  Combined matrix: {combined.shape[0]} regulons x {combined.shape[1]} groups")

# =============================================================================
# 6. Z-score and plot heatmap
# =============================================================================
print("\nStep 6: Plotting combined heatmap...")

# Z-score per regulon (row)
row_mean = combined.mean(axis=1)
row_std = combined.std(axis=1).replace(0, 1)
zmat = ((combined.T - row_mean) / row_std).T.clip(-ZSCORE_CAP, ZSCORE_CAP)

# Order regulons by cluster
regulon_order = []
for cl in EOS_CLUSTERS:
    tfs = [tf for tf in shared_regulons if selected_tfs.get(tf) == cl]
    regulon_order.extend(sorted(tfs))
# Add any remaining
for tf in shared_regulons:
    if tf not in regulon_order:
        regulon_order.append(tf)

zmat = zmat.loc[regulon_order]

# Column order
col_order = [c for c in EOS_CLUSTERS if c in zmat.columns] + \
            [c for c in ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"] if c in zmat.columns]
zmat = zmat[col_order]

# Plot
n_rows = len(regulon_order)
n_cols = len(col_order)
n_eos = len([c for c in EOS_CLUSTERS if c in col_order])

cell_h = 0.4
cell_w = 0.65
fig_h = n_rows * cell_h + 2.5
fig_w = n_cols * cell_w + 3.5

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

im = ax.imshow(zmat.values, cmap=ZSCORE_CMAP, aspect="auto",
               vmin=-ZSCORE_CAP, vmax=ZSCORE_CAP)

ax.set_xticks(range(n_cols))
ax.set_xticklabels(col_order, fontsize=10, fontweight="bold", rotation=45, ha="right")
ax.set_yticks(range(n_rows))
ax.set_yticklabels(regulon_order, fontsize=9, fontstyle="italic")
ax.tick_params(axis="both", length=0)

# Vertical separator between our clusters and Hu depots
if n_eos < n_cols:
    ax.axvline(x=n_eos - 0.5, color="black", linewidth=2.0)

# Horizontal separators between cluster groups
y_pos = 0
for cl in EOS_CLUSTERS:
    n_tfs = len([tf for tf in regulon_order if selected_tfs.get(tf) == cl])
    if n_tfs > 0:
        y_pos += n_tfs
        if y_pos < n_rows:
            ax.axhline(y=y_pos - 0.5, color="white", linewidth=1.5)

# Cluster sidebar
sidebar_colors = {"C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
                  "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4"}
for i, tf in enumerate(regulon_order):
    cl = selected_tfs.get(tf, "")
    if cl in sidebar_colors:
        ax.add_patch(plt.Rectangle((-0.8, i - 0.5), 0.3, 1,
                                    color=sidebar_colors[cl], clip_on=False))

# Group labels
ax.text(n_eos / 2 - 0.5, -1.5, "Qiu (HFD Adipose)",
        ha="center", fontsize=11, fontweight="bold", color="#E63946")
if n_eos < n_cols:
    ax.text(n_eos + (n_cols - n_eos) / 2 - 0.5, -1.5, "Hu (Healthy Adipose)",
            ha="center", fontsize=11, fontweight="bold", color="#4DAF4A")

ax.set_title("Cross-Dataset Regulon Activity Comparison",
             fontsize=13, fontweight="bold", pad=25)

# Colorbar
cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.04)
cbar.set_label("Z-score", fontsize=9)
cbar.set_ticks([-2, -1, 0, 1, 2])

out_pdf = os.path.join(FIG_DIR, "41_regulon_comparison_qiu_hu.pdf")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_pdf}")

# =============================================================================
# 7. Spearman correlation: which depot is closest to which cluster?
# =============================================================================
print("\nStep 7: Cross-dataset regulon correlation...")

cor_mat = combined.T.corr(method="spearman")

# Print depot-cluster correlations
print("\n  Regulon-based similarity (Spearman):")
for depot in ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]:
    if depot not in cor_mat.columns:
        continue
    for cl in EOS_CLUSTERS:
        if cl not in cor_mat.columns:
            continue
        r = cor_mat.loc[depot, cl]
        print(f"    {depot} vs {cl}: r = {r:.3f}")

# Save correlation matrix
cor_csv = os.path.join(BASE, "scenic", "41_regulon_correlation_qiu_hu.csv")
cor_mat.to_csv(cor_csv)
print(f"  Saved: {cor_csv}")

# =============================================================================
# 8. Correlation heatmap
# =============================================================================
print("\nStep 8: Correlation heatmap...")

fig, ax = plt.subplots(figsize=(7, 5.5))

# Subset to cross-dataset comparisons
depots = [d for d in ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"] if d in cor_mat.columns]
clusters = [c for c in EOS_CLUSTERS if c in cor_mat.columns]
cross_cor = cor_mat.loc[depots, clusters]

im = ax.imshow(cross_cor.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")

ax.set_xticks(range(len(clusters)))
ax.set_xticklabels(clusters, fontsize=11, fontweight="bold")
ax.set_yticks(range(len(depots)))
ax.set_yticklabels([d.replace("Hu_", "") for d in depots], fontsize=11, fontweight="bold")

for i in range(len(depots)):
    for j in range(len(clusters)):
        ax.text(j, i, f"{cross_cor.values[i, j]:.2f}",
                ha="center", va="center", fontsize=10,
                color="white" if abs(cross_cor.values[i, j]) > 0.5 else "black")

ax.set_title("Regulon Activity Similarity\nQiu Clusters vs Hu Adipose Depots",
             fontsize=12, fontweight="bold")
ax.set_xlabel("Qiu Eosinophil Cluster", fontsize=11)
ax.set_ylabel("Hu Adipose Depot", fontsize=11)

cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
cbar.set_label("Spearman r", fontsize=9)

out_pdf2 = os.path.join(FIG_DIR, "41_regulon_correlation_heatmap.pdf")
fig.savefig(out_pdf2, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_pdf2}")

print("\n41 complete.")
