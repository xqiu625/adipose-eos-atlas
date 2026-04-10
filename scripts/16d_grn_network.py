#!/usr/bin/env python3
"""
16d_grn_network.py
Eosinophil Gene Regulatory Network visualization (Figure 5B equivalent).

Nodes = TF regulons colored by dominant Eos subtype.
Node size = mean AUC activity.
Edges = Pearson co-activity correlation >= threshold across subsampled cells.

SLURM: --partition=batch --mem=32G --cpus-per-task=4 --time=00:30:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
REG_CSV    = os.path.join(BASE, "scenic", "eos_unique_regulons_per_cluster.csv")
FIG_DIR    = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]

TOP_PER_CLUSTER = 4       # TFs per cluster to include as nodes
CORR_THRESH     = 0.4     # Pearson r threshold for edges (higher = fewer edges)
SUBSAMPLE_N     = 500     # cells to subsample for correlation
RANDOM_SEED     = 42
MIN_AUC         = 0.02    # minimum mean AUC to include a TF

# Cluster colors — designed to be visually distinct
CLUSTER_COLORS = {
    "C0": "#7eb3d8",   # light blue
    "C1": "#e86850",   # coral/red
    "C2": "#f5a623",   # orange
    "C3": "#c8a2d0",   # lavender
    "C4": "#4fb08a",   # teal
    "C5": "#2c5a9e",   # dark blue
}

CLUSTER_LABELS = {
    "C0": "Eos C0", "C1": "Eos C1", "C2": "Eos C2",
    "C3": "Eos C3", "C4": "Eos C4", "C5": "Eos C5",
}

print("=" * 60)
print("16d_grn_network.py")
print("=" * 60)

np.random.seed(RANDOM_SEED)

# =============================================================================
# Step 1: Load per-cluster AUCell and select TFs
# =============================================================================
print("\nStep 1: Loading AUCell data and selecting TFs...")

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

# Load all per-cluster AUCell matrices
eos_data = {}
for cl in EOS_CLUSTERS:
    auc_path = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(auc_path):
        df = pd.read_csv(auc_path, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T   # deduplicate TFs
        eos_data[cl] = df
        print(f"  {cl}: {df.shape[0]} cells x {df.shape[1]} TFs")

# Build unified TF set and compute per-cluster means
all_tfs = set()
for df in eos_data.values():
    all_tfs.update(df.columns)
all_tfs = sorted(all_tfs)

cluster_means = {}
for cl in EOS_CLUSTERS:
    if cl not in eos_data:
        continue
    df = eos_data[cl]
    means = {}
    for tf in all_tfs:
        means[tf] = df[tf].mean() if tf in df.columns else 0.0
    cluster_means[cl] = means

mean_df = pd.DataFrame(cluster_means)  # tfs x clusters

# Select top TFs per cluster by uniqueness
selected_tfs = {}  # tf -> cluster assignment
for cl in EOS_CLUSTERS:
    if cl not in mean_df.columns:
        continue
    scores = []
    for tf in all_tfs:
        vals = [mean_df.loc[tf, c] for c in EOS_CLUSTERS if c in mean_df.columns]
        max_val = mean_df.loc[tf, cl]
        if max_val < MIN_AUC:
            continue
        # Only consider if this cluster is the max
        if cl != mean_df.columns[np.argmax(mean_df.loc[tf].values)]:
            continue
        sorted_vals = sorted(vals, reverse=True)
        uniqueness = sorted_vals[0] - sorted_vals[1] if len(sorted_vals) >= 2 else sorted_vals[0]
        scores.append((tf, max_val, uniqueness))

    scores.sort(key=lambda x: x[2], reverse=True)
    for tf, auc_val, uniq in scores[:TOP_PER_CLUSTER]:
        selected_tfs[tf] = cl

print(f"\nSelected {len(selected_tfs)} TFs:")
for cl in EOS_CLUSTERS:
    tfs = [tf for tf, c in selected_tfs.items() if c == cl]
    if tfs:
        print(f"  {cl}: {', '.join(tfs)}")

# =============================================================================
# Step 2: Merge AUCell and compute correlations
# =============================================================================
print("\nStep 2: Computing co-activity correlations...")

# Merge all cells, align to selected TFs
all_parts = []
for cl in EOS_CLUSTERS:
    if cl not in eos_data:
        continue
    df = eos_data[cl]
    sub = pd.DataFrame(index=df.index)
    for tf in selected_tfs:
        sub[tf] = df[tf] if tf in df.columns else 0.0
    all_parts.append(sub)

merged = pd.concat(all_parts)
print(f"  Total cells: {merged.shape[0]}")

# Subsample
if len(merged) > SUBSAMPLE_N:
    idx = np.random.choice(len(merged), SUBSAMPLE_N, replace=False)
    sampled = merged.iloc[idx]
else:
    sampled = merged
print(f"  Subsampled: {sampled.shape[0]} cells")

# Pairwise Pearson correlation
tf_list = list(selected_tfs.keys())
corr_matrix = sampled[tf_list].corr(method="pearson")

# =============================================================================
# Step 3: Build network
# =============================================================================
print("\nStep 3: Building network...")

G = nx.Graph()

# Add nodes with attributes
global_means = merged[tf_list].mean()
for tf in tf_list:
    G.add_node(tf, cluster=selected_tfs[tf], auc=global_means[tf])

# Add edges where correlation >= threshold
edges_added = 0
for i, tf1 in enumerate(tf_list):
    for j, tf2 in enumerate(tf_list):
        if j <= i:
            continue
        r = corr_matrix.loc[tf1, tf2]
        if abs(r) >= CORR_THRESH:
            G.add_edge(tf1, tf2, weight=r)
            edges_added += 1

print(f"  Nodes: {G.number_of_nodes()}")
print(f"  Edges: {edges_added}")

# =============================================================================
# Step 4: Layout and plot
# =============================================================================
print("\nStep 4: Plotting...")

fig, ax = plt.subplots(figsize=(12, 11))

# Spring layout — high k for spacing, many iterations to stabilize
pos = nx.spring_layout(G, k=5.0, iterations=200, seed=RANDOM_SEED, scale=2.0)

# Node colors and sizes
node_colors = [CLUSTER_COLORS[G.nodes[n]["cluster"]] for n in G.nodes]
auc_vals = np.array([G.nodes[n]["auc"] for n in G.nodes])
node_sizes = 150 + auc_vals * 1200   # smaller nodes to reduce overlap

# Draw edges — thin and transparent
for u, v in G.edges:
    r = abs(G[u][v]["weight"])
    alpha = 0.15 + r * 0.3
    lw = 0.5 + r * 1.5
    x = [pos[u][0], pos[v][0]]
    y = [pos[u][1], pos[v][1]]
    ax.plot(x, y, color="#aaaaaa", linewidth=lw, alpha=alpha, zorder=1)

# Draw nodes
ax.scatter([pos[n][0] for n in G.nodes],
           [pos[n][1] for n in G.nodes],
           s=node_sizes, c=node_colors,
           edgecolors="white", linewidths=1.0, zorder=2)

# Node labels — offset below each node to avoid overlap
for n in G.nodes:
    r = np.sqrt(node_sizes[list(G.nodes).index(n)]) / 180  # offset proportional to size
    ax.annotate(n, (pos[n][0], pos[n][1] - r * 0.12),
                fontsize=8, fontweight="bold",
                ha="center", va="top", zorder=3)

# Cluster labels — place at outermost node, pushed outward from graph center
graph_cx = np.mean([pos[n][0] for n in G.nodes])
graph_cy = np.mean([pos[n][1] for n in G.nodes])
for cl in EOS_CLUSTERS:
    cl_nodes = [n for n in G.nodes if G.nodes[n]["cluster"] == cl]
    if not cl_nodes:
        continue
    cx = np.mean([pos[n][0] for n in cl_nodes])
    cy = np.mean([pos[n][1] for n in cl_nodes])
    # Push label outward from graph center
    dx = cx - graph_cx
    dy = cy - graph_cy
    dist = max(np.sqrt(dx**2 + dy**2), 0.01)
    offset = 0.35
    lx = cx + dx / dist * offset
    ly = cy + dy / dist * offset
    ax.text(lx, ly, CLUSTER_LABELS[cl],
            fontsize=13, fontweight="bold",
            color=CLUSTER_COLORS[cl], ha="center", va="center",
            zorder=4)

# Node size legend — bottom right corner using axes fraction
ax_xlim = ax.get_xlim()
ax_ylim = ax.get_ylim()
leg_x = ax_xlim[1] - 0.3
leg_y = ax_ylim[0] + 0.2
ax.text(leg_x, leg_y + 0.35, "Node size", fontsize=9,
        fontweight="bold", va="bottom", ha="center")
for auc_val, label in [(0, "AUC =0"), (0.5, "AUC =0.5"), (1.0, "AUC =1")]:
    size = 150 + auc_val * 1200
    ax.scatter([leg_x], [leg_y], s=size, c="#aaaaaa",
               edgecolors="white", linewidths=0.8, zorder=2)
    ax.text(leg_x + 0.25, leg_y, label, fontsize=8,
            va="center", ha="left", zorder=3)
    leg_y -= 0.25

ax.set_title("B        Eosinophil Gene Regulatory Network",
             fontsize=14, fontweight="bold", pad=15, loc="left")
ax.axis("off")

# Expand limits for labels
xlim = ax.get_xlim()
ylim = ax.get_ylim()
margin = 0.5
ax.set_xlim(xlim[0] - margin, xlim[1] + margin)
ax.set_ylim(ylim[0] - margin, ylim[1] + margin)

out_pdf = os.path.join(FIG_DIR, "Fig5_grn_network.pdf")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved: {out_pdf}")

# Save edge list
edge_df = pd.DataFrame([(u, v, G[u][v]["weight"]) for u, v in G.edges],
                        columns=["TF1", "TF2", "pearson_r"])
edge_df.to_csv(os.path.join(BASE, "scenic", "grn_network_edges.csv"), index=False)
print("16d complete.")
