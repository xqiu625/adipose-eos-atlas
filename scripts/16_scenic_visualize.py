#!/usr/bin/env python3
# ==============================================================================
# 16_scenic_visualize.py
# Visualize pySCENIC results: regulon activity heatmap, dot plot, network,
# and regulon specificity scores (RSS) per cluster
# ==============================================================================
#
# SLURM: --partition=batch --mem=64G --cpus-per-task=4 --time=2:00:00
# Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
#
# Inputs:
#   - eos_scenic_aucell.csv   (AUCell scores from pySCENIC)
#   - eos_for_scenic.h5ad     (original expression + metadata)
#
# Outputs:
#   - scenic/figures/*.pdf     (heatmap, dot plot, network)
#   - scenic/regulon_specificity_scores.csv
# ==============================================================================

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import jensenshannon
from collections import defaultdict
import json

warnings.filterwarnings("ignore")

# -- Publication style ---------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,  # TrueType fonts in PDF
    "ps.fonttype": 42,
})

# Eos cluster palette (matching R scripts)
EOS_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4", "C6": "#2F3C7E"
}
EOS_CMAP = sns.color_palette([EOS_COLORS[f"C{i}"] for i in range(7)])

# -- Configuration -------------------------------------------------------------
BASE = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
AUCELL_CSV = os.path.join(BASE, "scenic", "eos_scenic_aucell.csv")
H5AD_PATH = os.path.join(BASE, "scenic", "eos_for_scenic.h5ad")
FIG_DIR = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

print("=" * 60)
print("16_scenic_visualize.py")
print("Visualize pySCENIC results")
print("=" * 60)
print()

# -- Step 1: Load data ---------------------------------------------------------
print("Step 1: Loading data...")

# Load original h5ad for metadata and UMAP
print(f"  Loading h5ad: {H5AD_PATH}")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Load SCENIC AUCell CSV output
print(f"  Loading AUCell CSV: {AUCELL_CSV}")
auc_mtx = pd.read_csv(AUCELL_CSV, index_col=0)
# Clean column names: remove trailing "(+)" from regulon names if present
auc_mtx.columns = [c.strip() for c in auc_mtx.columns]

print(f"  AUCell matrix: {auc_mtx.shape[0]} cells x {auc_mtx.shape[1]} regulons")
print(f"  Regulons: {', '.join(auc_mtx.columns[:5])}... ({len(auc_mtx.columns)} total)")

# Align cells between h5ad and AUCell matrix
common_cells = adata.obs_names.intersection(auc_mtx.index)
print(f"  Common cells: {len(common_cells)}")

if len(common_cells) == 0:
    raise ValueError("No overlapping cell IDs between h5ad and loom files. "
                     "Check that cell IDs match.")

adata = adata[common_cells, :].copy()
auc_mtx = auc_mtx.loc[common_cells, :]

# Add cluster info to AUCell matrix
cluster_col = "cluster" if "cluster" in adata.obs.columns else "eos_cluster"
auc_mtx["cluster"] = adata.obs[cluster_col].values

print()

# -- Step 2: Regulon Specificity Score (RSS) -----------------------------------
print("Step 2: Computing Regulon Specificity Scores (RSS)...")


def compute_rss(auc_matrix, clusters):
    """
    Compute Regulon Specificity Score using Jensen-Shannon divergence.
    For each regulon, compare the distribution of AUC scores in each cluster
    versus the overall background.
    """
    regulon_cols = [c for c in auc_matrix.columns if c != "cluster"]
    unique_clusters = sorted(clusters.unique())
    rss_dict = {}

    for cluster in unique_clusters:
        rss_dict[cluster] = {}
        # Binary indicator: 1 if cell is in cluster, 0 otherwise
        cluster_indicator = (clusters == cluster).astype(float).values
        cluster_indicator = cluster_indicator / cluster_indicator.sum()

        for regulon in regulon_cols:
            # Regulon activity distribution (normalized to sum to 1)
            regulon_activity = auc_matrix[regulon].values.copy()
            regulon_activity = regulon_activity - regulon_activity.min()
            total = regulon_activity.sum()
            if total > 0:
                regulon_activity = regulon_activity / total
            else:
                regulon_activity = np.ones_like(regulon_activity) / len(regulon_activity)

            # RSS = 1 - JS divergence (higher = more specific)
            jsd = jensenshannon(cluster_indicator, regulon_activity)
            rss_dict[cluster][regulon] = 1.0 - jsd

    rss_df = pd.DataFrame(rss_dict)
    rss_df.index.name = "regulon"
    rss_df.columns.name = "cluster"
    return rss_df


clusters_series = auc_mtx["cluster"]
rss = compute_rss(auc_mtx, clusters_series)

# Save RSS table
rss_path = os.path.join(BASE, "scenic", "regulon_specificity_scores.csv")
rss.to_csv(rss_path)
print(f"  Saved RSS table: {rss_path}")
print(f"  RSS matrix: {rss.shape[0]} regulons x {rss.shape[1]} clusters")
print()

# -- Step 3: Top regulons per cluster ------------------------------------------
print("Step 3: Identifying top regulons per cluster...")

top_n = 5
top_regulons_per_cluster = {}
all_top_regulons = set()

for cluster in rss.columns:
    top = rss[cluster].nlargest(top_n)
    top_regulons_per_cluster[cluster] = list(top.index)
    all_top_regulons.update(top.index)
    print(f"  {cluster}: {', '.join(top.index)}")

all_top_regulons = sorted(all_top_regulons)
print(f"  Total unique top regulons: {len(all_top_regulons)}")
print()

# -- Step 4: Heatmap of top regulons (z-scored AUCell) -------------------------
print("Step 4: Plotting heatmap (z-scored AUCell, top regulons)...")

# Subset AUCell to top regulons
auc_top = auc_mtx[all_top_regulons].copy()

# Z-score across cells
auc_z = (auc_top - auc_top.mean()) / auc_top.std()

# Compute mean z-score per cluster
auc_z["cluster"] = clusters_series.values
cluster_means = auc_z.groupby("cluster").mean()

# Sort regulons by cluster with highest specificity
regulon_order = []
for cluster in sorted(cluster_means.index):
    if cluster in top_regulons_per_cluster:
        for reg in top_regulons_per_cluster[cluster]:
            if reg not in regulon_order:
                regulon_order.append(reg)
# Add any remaining
for reg in all_top_regulons:
    if reg not in regulon_order:
        regulon_order.append(reg)

cluster_means = cluster_means[regulon_order]

fig, ax = plt.subplots(figsize=(max(12, len(regulon_order) * 0.5),
                                max(5, len(cluster_means) * 0.8)))
sns.heatmap(
    cluster_means.T,
    cmap="RdBu_r",
    center=0,
    linewidths=0.5,
    linecolor="white",
    ax=ax,
    cbar_kws={"label": "Mean z-scored AUCell"},
    xticklabels=True,
    yticklabels=True
)
ax.set_title("Top Regulon Activity per Eosinophil Cluster", fontsize=14, pad=15)
ax.set_xlabel("Cluster", fontsize=12)
ax.set_ylabel("Regulon", fontsize=12)
ax.tick_params(axis="y", labelsize=9)
ax.tick_params(axis="x", labelsize=10)
plt.tight_layout()

heatmap_path = os.path.join(FIG_DIR, "16_regulon_heatmap.pdf")
fig.savefig(heatmap_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {heatmap_path}")

# -- Step 5: Dot plot of regulon activity by cluster ---------------------------
print("Step 5: Plotting dot plot (regulon activity by cluster)...")

# For each cluster-regulon pair: mean AUCell and fraction of cells > threshold
threshold = auc_mtx[all_top_regulons].quantile(0.75).mean()

dot_data = []
for cluster in sorted(clusters_series.unique()):
    mask = clusters_series == cluster
    for reg in regulon_order:
        vals = auc_mtx.loc[mask, reg]
        dot_data.append({
            "cluster": cluster,
            "regulon": reg,
            "mean_auc": vals.mean(),
            "frac_active": (vals > threshold).mean()
        })

dot_df = pd.DataFrame(dot_data)

fig, ax = plt.subplots(figsize=(max(12, len(regulon_order) * 0.5),
                                max(5, len(sorted(clusters_series.unique())) * 0.8)))

# Create grid
clusters_sorted = sorted(clusters_series.unique())
for i, cluster in enumerate(clusters_sorted):
    sub = dot_df[dot_df["cluster"] == cluster]
    for j, reg in enumerate(regulon_order):
        row = sub[sub["regulon"] == reg]
        if len(row) > 0:
            size = row["frac_active"].values[0] * 300
            color = row["mean_auc"].values[0]
            ax.scatter(j, i, s=max(size, 5), c=color, cmap="viridis",
                       vmin=dot_df["mean_auc"].min(),
                       vmax=dot_df["mean_auc"].max(),
                       edgecolors="grey", linewidth=0.5)

ax.set_xticks(range(len(regulon_order)))
ax.set_xticklabels(regulon_order, rotation=90, fontsize=8)
ax.set_yticks(range(len(clusters_sorted)))
ax.set_yticklabels(clusters_sorted, fontsize=10)
ax.set_xlabel("Regulon", fontsize=12)
ax.set_ylabel("Cluster", fontsize=12)
ax.set_title("Regulon Activity Dot Plot", fontsize=14, pad=15)

# Add colorbar for mean AUC
sm = plt.cm.ScalarMappable(
    cmap="viridis",
    norm=plt.Normalize(vmin=dot_df["mean_auc"].min(),
                       vmax=dot_df["mean_auc"].max())
)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.04)
cbar.set_label("Mean AUCell score", fontsize=10)

plt.tight_layout()

dotplot_path = os.path.join(FIG_DIR, "16_regulon_dotplot.pdf")
fig.savefig(dotplot_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {dotplot_path}")

# -- Step 6: Network visualization (top regulons + target genes) ---------------
print("Step 6: Plotting regulon-target gene networks per cluster...")

# Extract regulon-target gene mappings from regulons.csv
REGULONS_CSV = os.path.join(BASE, "scenic", "regulons.csv")
regulon_targets = {}
try:
    print(f"  Loading regulon targets from: {REGULONS_CSV}")
    reg_df = pd.read_csv(REGULONS_CSV)
    # pySCENIC ctx output: columns are TF, MotifID, AUC, ..., TargetGenes
    if "TF" in reg_df.columns:
        for tf_name in reg_df["TF"].unique():
            tf_rows = reg_df[reg_df["TF"] == tf_name]
            # TargetGenes column contains comma-separated gene lists or JSON
            targets = set()
            for _, row in tf_rows.iterrows():
                if "TargetGenes" in row and pd.notna(row["TargetGenes"]):
                    tg = str(row["TargetGenes"])
                    # Try parsing as list of (gene, score) tuples
                    try:
                        gene_list = eval(tg) if tg.startswith("[") else tg.split(",")
                        for item in gene_list:
                            if isinstance(item, (list, tuple)):
                                targets.add(str(item[0]).strip())
                            else:
                                targets.add(str(item).strip())
                    except Exception:
                        targets.update(t.strip() for t in tg.split(","))
            if targets:
                reg_key = f"{tf_name}(+)"
                regulon_targets[reg_key] = list(targets)[:50]  # cap at 50
    print(f"  Found target genes for {len(regulon_targets)} regulons")
except Exception as e:
    print(f"  Warning: Could not extract regulon targets: {e}")

if regulon_targets:
    print(f"  Found target genes for {len(regulon_targets)} regulons")

    # Plot network for each cluster's top regulons
    try:
        import networkx as nx

        for cluster in sorted(top_regulons_per_cluster.keys()):
            top_regs = top_regulons_per_cluster[cluster]
            fig, ax = plt.subplots(figsize=(10, 8))

            G = nx.Graph()

            for reg in top_regs:
                # Clean regulon name (remove (+) suffix if present)
                reg_clean = reg.rstrip("(+)").rstrip("(-)")
                G.add_node(reg, node_type="regulon")

                targets = regulon_targets.get(reg, [])
                # Show top 10 target genes per regulon
                for target in targets[:10]:
                    G.add_node(target, node_type="target")
                    G.add_edge(reg, target)

            if len(G.nodes()) == 0:
                plt.close(fig)
                continue

            # Layout
            pos = nx.spring_layout(G, k=2.0, iterations=50, seed=42)

            # Draw regulon nodes (larger, colored)
            reg_nodes = [n for n, d in G.nodes(data=True)
                         if d.get("node_type") == "regulon"]
            target_nodes = [n for n, d in G.nodes(data=True)
                            if d.get("node_type") == "target"]

            nx.draw_networkx_nodes(G, pos, nodelist=reg_nodes, node_size=800,
                                   node_color="#E84A5F", alpha=0.9, ax=ax)
            nx.draw_networkx_nodes(G, pos, nodelist=target_nodes, node_size=200,
                                   node_color="#6CA6CD", alpha=0.7, ax=ax)
            nx.draw_networkx_edges(G, pos, alpha=0.3, ax=ax)
            nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

            ax.set_title(f"Top Regulon Networks — {cluster}", fontsize=14)
            ax.axis("off")
            plt.tight_layout()

            net_path = os.path.join(FIG_DIR,
                                    f"16_regulon_network_{cluster}.pdf")
            fig.savefig(net_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {net_path}")

    except ImportError:
        print("  Warning: networkx not available; skipping network plots.")
else:
    print("  No regulon-target mappings found; skipping network plots.")

# -- Step 7: RSS panel plot ----------------------------------------------------
print("Step 7: Plotting RSS panel...")

n_clusters = len(rss.columns)
n_cols = min(4, n_clusters)
n_rows = int(np.ceil(n_clusters / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
if n_clusters == 1:
    axes = np.array([[axes]])
axes = np.atleast_2d(axes)

for idx, cluster in enumerate(sorted(rss.columns)):
    row_idx = idx // n_cols
    col_idx = idx % n_cols
    ax = axes[row_idx, col_idx]

    scores = rss[cluster].sort_values(ascending=False)
    top_labels = scores.head(top_n)

    ax.scatter(range(len(scores)), scores.values, s=10, c="grey", alpha=0.5)

    # Highlight top regulons
    for rank, (reg, val) in enumerate(top_labels.items()):
        overall_rank = list(scores.index).index(reg)
        ax.scatter(overall_rank, val, s=40, c="#E84A5F", zorder=5)
        ax.annotate(reg, (overall_rank, val), fontsize=6,
                    xytext=(5, 5), textcoords="offset points")

    ax.set_title(cluster, fontsize=12, fontweight="bold")
    ax.set_xlabel("Regulon rank", fontsize=9)
    ax.set_ylabel("RSS", fontsize=9)
    ax.tick_params(labelsize=8)

# Hide unused subplots
for idx in range(n_clusters, n_rows * n_cols):
    row_idx = idx // n_cols
    col_idx = idx % n_cols
    axes[row_idx, col_idx].set_visible(False)

plt.suptitle("Regulon Specificity Score (RSS) per Cluster",
             fontsize=16, y=1.02)
plt.tight_layout()

rss_panel_path = os.path.join(FIG_DIR, "16_rss_panel.pdf")
fig.savefig(rss_panel_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {rss_panel_path}")

# -- Summary -------------------------------------------------------------------
print()
print("=" * 60)
print("Summary")
print("=" * 60)
print(f"  Cells:            {auc_mtx.shape[0]}")
print(f"  Regulons:         {auc_mtx.shape[1] - 1}")  # minus cluster column
print(f"  Clusters:         {n_clusters}")
print(f"  Top regulons/cl:  {top_n}")
print()
print("Output files:")
print(f"  {rss_path}")
print(f"  {heatmap_path}")
print(f"  {dotplot_path}")
print(f"  {rss_panel_path}")
if regulon_targets:
    print(f"  {FIG_DIR}/16_regulon_network_*.pdf")
print()
print("Done! 16_scenic_visualize.py completed successfully.")
