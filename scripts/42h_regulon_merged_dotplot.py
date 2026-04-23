#!/usr/bin/env python3
"""
42h_regulon_merged_dotplot.py
Merged regulon dot plot in Seurat DotPlot style (matching Fig 2B):
  - X axis = regulons (grouped by category)
  - Y axis = groups (C0-C5, Hu_eWAT, Hu_scWAT, Hu_BAT)
  - Dot size = % cells with regulon active (AUCell > threshold)
  - Dot color = Z-scored mean AUCell (blue-white-red diverging)

Combines:
  - Top cluster-specific shared regulons (from 42b)
  - HFD-only regulons (from 42g)
  - Shared HFD-enriched regulons (from 42g)

SLURM: --partition=batch --mem=32G --cpus-per-task=4 --time=00:30:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
HU_SCENIC  = os.path.join(BASE, "scenic", "hu_scenic")
FIG_DIR    = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]
HU_DEPOTS = ["eWAT", "scWAT", "BAT"]
ALL_GROUPS = EOS_CLUSTERS + ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]

AUC_THRESH = 0.01
ZSCORE_CAP = 2.0
TOP_PER_CLUSTER = 3

EXCLUDE = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

print("=" * 60)
print("42h_regulon_merged_dotplot.py")
print("=" * 60)

# =============================================================================
# 1. Load all AUCell data
# =============================================================================
print("\nStep 1: Loading data...")

qiu_parts, qiu_groups = [], []
for cl in EOS_CLUSTERS:
    p = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(p):
        df = pd.read_csv(p, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T
        qiu_parts.append(df)
        qiu_groups.append(pd.Series(cl, index=df.index))

qiu_auc = pd.concat(qiu_parts, join="outer").fillna(0)
qiu_group = pd.concat(qiu_groups)
qiu_regulons = set(qiu_auc.columns)

hu_parts, hu_groups = [], []
for depot in HU_DEPOTS:
    p = os.path.join(HU_SCENIC, depot, "auc_matrix.csv")
    if os.path.exists(p):
        df = pd.read_csv(p, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T
        hu_parts.append(df)
        hu_groups.append(pd.Series(f"Hu_{depot}", index=df.index))

hu_auc = pd.concat(hu_parts, join="outer").fillna(0)
hu_group = pd.concat(hu_groups)
hu_regulons = set(hu_auc.columns)

shared = sorted(qiu_regulons & hu_regulons)
qiu_only = sorted(qiu_regulons - hu_regulons)
hu_only = sorted(hu_regulons - qiu_regulons)

print(f"  Qiu: {len(qiu_auc)} cells, {len(qiu_regulons)} regulons")
print(f"  Hu: {len(hu_auc)} cells, {len(hu_regulons)} regulons")
print(f"  Shared: {len(shared)}, Qiu-only: {len(qiu_only)}, Hu-only: {len(hu_only)}")

# =============================================================================
# 2. Compute per-group mean AUCell and % active
# =============================================================================
print("\nStep 2: Computing per-group statistics...")

all_tfs = sorted(qiu_regulons | hu_regulons)
mean_scores = {}
pct_active = {}

for grp in EOS_CLUSTERS:
    mask = qiu_group == grp
    cells = qiu_auc.loc[mask]
    for tf in all_tfs:
        if tf in qiu_auc.columns:
            vals = cells[tf].values
            mean_scores[(tf, grp)] = np.mean(vals)
            pct_active[(tf, grp)] = np.mean(vals > AUC_THRESH) * 100
        else:
            mean_scores[(tf, grp)] = 0.0
            pct_active[(tf, grp)] = 0.0

for depot in HU_DEPOTS:
    grp = f"Hu_{depot}"
    mask = hu_group == grp
    cells = hu_auc.loc[mask]
    for tf in all_tfs:
        if tf in hu_auc.columns:
            vals = cells[tf].values
            mean_scores[(tf, grp)] = np.mean(vals)
            pct_active[(tf, grp)] = np.mean(vals > AUC_THRESH) * 100
        else:
            mean_scores[(tf, grp)] = 0.0
            pct_active[(tf, grp)] = 0.0

# =============================================================================
# 3. Select regulons: 3 categories
# =============================================================================
print("\nStep 3: Selecting regulons...")

# --- A) Cluster-specific SHARED regulons (from 42b logic) ---
shared_filtered = [t for t in shared if t not in EXCLUDE]
mean_df = pd.DataFrame({
    grp: {tf: mean_scores[(tf, grp)] for tf in shared_filtered}
    for grp in ALL_GROUPS
})

cluster_specific = {}
for cl in EOS_CLUSTERS:
    scores = []
    for tf in shared_filtered:
        max_val = mean_df.loc[tf, cl]
        if max_val < 0.02:
            continue
        if cl != mean_df.columns[np.argmax(mean_df.loc[tf].values)]:
            continue
        sorted_vals = sorted(mean_df.loc[tf].values, reverse=True)
        uniqueness = sorted_vals[0] - sorted_vals[1] if len(sorted_vals) >= 2 else sorted_vals[0]
        scores.append((tf, uniqueness))
    scores.sort(key=lambda x: x[1], reverse=True)
    for tf, _ in scores[:TOP_PER_CLUSTER]:
        if tf not in cluster_specific:
            cluster_specific[tf] = cl

cat_a = list(cluster_specific.keys())
print(f"  A) Cluster-specific shared: {len(cat_a)} regulons")
for cl in EOS_CLUSTERS:
    tfs = [t for t, c in cluster_specific.items() if c == cl]
    if tfs:
        print(f"     {cl}: {', '.join(tfs)}")

# --- B) HFD-only regulons (absent in Hu) ---
qiu_only_filtered = [t for t in qiu_only if t not in EXCLUDE]
qiu_only_max = {}
for tf in qiu_only_filtered:
    qiu_only_max[tf] = max(mean_scores.get((tf, cl), 0) for cl in EOS_CLUSTERS)
top_qiu_only = sorted(qiu_only_max, key=qiu_only_max.get, reverse=True)[:8]
cat_b = top_qiu_only
print(f"  B) HFD-only: {cat_b}")

# --- C) Shared but HFD-enriched (high fold-change Qiu/Hu) ---
enrichment = {}
for tf in shared_filtered:
    if tf in cat_a:
        continue  # already selected
    qiu_max_val = max(mean_scores.get((tf, cl), 0) for cl in EOS_CLUSTERS)
    hu_max_val = max(mean_scores.get((tf, f"Hu_{d}"), 0) for d in HU_DEPOTS)
    if qiu_max_val > 0.03:
        enrichment[tf] = qiu_max_val / (hu_max_val + 0.001)
top_enriched = sorted(enrichment, key=enrichment.get, reverse=True)[:6]
cat_c = top_enriched
print(f"  C) Shared HFD-enriched: {cat_c}")

# Combine all — ordered by category
all_selected = cat_a + cat_b + cat_c
# Remove duplicates
seen = set()
regulon_order = []
for tf in all_selected:
    if tf not in seen:
        regulon_order.append(tf)
        seen.add(tf)

print(f"\n  Total unique regulons: {len(regulon_order)}")

# Category boundaries for visual separators
n_cat_a = len([t for t in regulon_order if t in cat_a])
n_cat_b = len([t for t in regulon_order if t in cat_b and t not in cat_a])

# =============================================================================
# 4. Build matrices
# =============================================================================
print("\nStep 4: Building matrices...")

n_reg = len(regulon_order)
n_grp = len(ALL_GROUPS)

raw_mat = np.zeros((n_grp, n_reg))
pct_mat = np.zeros((n_grp, n_reg))

for i, grp in enumerate(ALL_GROUPS):
    for j, tf in enumerate(regulon_order):
        raw_mat[i, j] = mean_scores.get((tf, grp), 0)
        pct_mat[i, j] = pct_active.get((tf, grp), 0)

# Z-score per regulon (across groups) for color
z_mat = np.zeros_like(raw_mat)
for j in range(n_reg):
    col = raw_mat[:, j]
    mu = np.mean(col)
    sd = np.std(col)
    if sd > 0:
        z_mat[:, j] = (col - mu) / sd
    else:
        z_mat[:, j] = 0
z_mat = np.clip(z_mat, -ZSCORE_CAP, ZSCORE_CAP)

# =============================================================================
# 5. Plot: Seurat DotPlot style (X=regulons, Y=groups)
# =============================================================================
print("\nStep 5: Generating merged dot plot...")

# Color map: blue-white-red (matching Fig 2B)
cmap = LinearSegmentedColormap.from_list(
    "bwr", ["#2166ac", "#f7f7f7", "#b2182b"])
norm = TwoSlopeNorm(vmin=-ZSCORE_CAP, vcenter=0, vmax=ZSCORE_CAP)

# Dot sizing — scaled for 6.6cm x 12.67cm figure
max_dot = 80
min_dot = 3
pct_max = max(1, np.max(pct_mat))

# 6.6 cm height x 12.67 cm width
figwidth = 12.67 / 2.54   # 4.99 inches
figheight = 6.6 / 2.54    # 2.60 inches

fig, ax = plt.subplots(figsize=(figwidth, figheight))

for i in range(n_grp):
    for j in range(n_reg):
        z = z_mat[i, j]
        p = pct_mat[i, j]
        size = min_dot + (p / pct_max) * (max_dot - min_dot)
        color = cmap(norm(z))
        ax.scatter(j, n_grp - 1 - i, s=size, c=[color],
                   edgecolors="#404040", linewidths=0.2, zorder=3)

# Axes
ax.set_xticks(range(n_reg))
ax.set_xticklabels(regulon_order, fontsize=4, fontstyle="italic",
                   rotation=45, ha="right")
ax.set_yticks(range(n_grp))
ylabels = list(reversed(ALL_GROUPS))
ylabels_fmt = []
for g in ylabels:
    if g.startswith("Hu_"):
        ylabels_fmt.append(g.replace("_", " "))
    else:
        ylabels_fmt.append(g)
ax.set_yticklabels(ylabels_fmt, fontsize=5, fontweight="bold")

# Horizontal separator between Qiu and Hu
n_hu = 3
ax.axhline(y=n_hu - 0.5, color="black", linewidth=0.8, linestyle="-")

# Vertical separators between categories
sep1 = n_cat_a - 0.5
sep2 = n_cat_a + n_cat_b - 0.5
ax.axvline(x=sep1, color="#404040", linewidth=0.5, linestyle="--", alpha=0.6)
ax.axvline(x=sep2, color="#404040", linewidth=0.5, linestyle="--", alpha=0.6)

# Category labels at top
cat_a_center = (n_cat_a - 1) / 2
cat_b_center = n_cat_a + n_cat_b / 2 - 0.5
cat_c_center = n_cat_a + n_cat_b + len(cat_c) / 2 - 0.5

ax.text(cat_a_center, n_grp + 0.2, "Cluster-Specific\n(Shared)",
        ha="center", fontsize=4, fontweight="bold", color="#333333")
ax.text(cat_b_center, n_grp + 0.2, "HFD-Only",
        ha="center", fontsize=4, fontweight="bold", color="#E63946")
ax.text(cat_c_center, n_grp + 0.2, "Shared\n(HFD-Enriched)",
        ha="center", fontsize=4, fontweight="bold", color="#984ea3")

# Dataset labels on left
ax.text(-1.5, (n_grp - 1 + n_hu) / 2, "Qiu\n(HFD)",
        ha="center", va="center", fontsize=5, fontweight="bold",
        color="#E63946", rotation=90)
ax.text(-1.5, (n_hu - 1) / 2, "Hu\n(Healthy)",
        ha="center", va="center", fontsize=5, fontweight="bold",
        color="#2A9D8F", rotation=90)

# Grid
ax.set_xlim(-0.5, n_reg - 0.5)
ax.set_ylim(-0.5, n_grp + 0.6)
ax.set_axisbelow(True)
for i in range(n_grp):
    ax.axhline(y=i, color="#e0e0e0", linewidth=0.2, zorder=0)
for j in range(n_reg):
    ax.axvline(x=j, color="#e0e0e0", linewidth=0.2, zorder=0)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", length=1.5, width=0.3, pad=1)

# Colorbar for Z-score
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cax = fig.add_axes([0.92, 0.42, 0.012, 0.30])
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("Average\nExpression", fontsize=4, labelpad=4)
cbar.set_ticks([-2, -1, 0, 1, 2])
cbar.ax.tick_params(labelsize=3.5, length=1, width=0.3)

# Size legend for % active — use scatter so sizes match exactly
legend_pcts = [25, 50, 75, 100]
for p in legend_pcts:
    size = min_dot + (p / pct_max) * (max_dot - min_dot)
    ax.scatter([], [], s=size, c="#404040", edgecolors="#404040",
               linewidths=0.2, label=f'{p}')

leg = ax.legend(title="Percent\nExpressed",
                loc="upper left", bbox_to_anchor=(1.06, 0.40),
                frameon=False, fontsize=3.5, title_fontsize=4,
                handletextpad=0.5, labelspacing=2.0,
                scatterpoints=1)

fig.subplots_adjust(right=0.88, top=0.86, bottom=0.22, left=0.10)
fig.savefig(os.path.join(FIG_DIR, "42h_merged_regulon_dotplot.pdf"),
            dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved: 42h_merged_regulon_dotplot.pdf")

# =============================================================================
# 6. Summary table
# =============================================================================
print("\n" + "=" * 60)
print("Regulon categories:")
print("=" * 60)
print(f"\nCluster-specific shared ({len(cat_a)}):")
for tf in cat_a:
    cl = cluster_specific.get(tf, "?")
    qmax = max(mean_scores.get((tf, c), 0) for c in EOS_CLUSTERS)
    hmax = max(mean_scores.get((tf, f"Hu_{d}"), 0) for d in HU_DEPOTS)
    print(f"  {tf} ({cl}): Qiu max={qmax:.3f}, Hu max={hmax:.3f}")

print(f"\nHFD-only ({len(cat_b)}):")
for tf in cat_b:
    qmax = max(mean_scores.get((tf, c), 0) for c in EOS_CLUSTERS)
    print(f"  {tf}: Qiu max={qmax:.3f} (absent in Hu)")

print(f"\nShared HFD-enriched ({len(cat_c)}):")
for tf in cat_c:
    qmax = max(mean_scores.get((tf, c), 0) for c in EOS_CLUSTERS)
    hmax = max(mean_scores.get((tf, f"Hu_{d}"), 0) for d in HU_DEPOTS)
    fc = qmax / (hmax + 0.001)
    print(f"  {tf}: Qiu max={qmax:.3f}, Hu max={hmax:.3f}, FC={fc:.1f}x")

print("\n42h complete.")
