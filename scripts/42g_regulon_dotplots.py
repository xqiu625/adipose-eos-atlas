#!/usr/bin/env python3
"""
42g_regulon_dotplots.py
Dot plots for regulon comparison (Seurat DotPlot style):
  - Dot size = % cells with AUCell > threshold (activity)
  - Dot color = mean AUCell score (intensity)
  - Y axis = regulons
  - X axis = groups (C0-C5, Hu_eWAT, Hu_scWAT, Hu_BAT)

1. HFD-enriched regulons dot plot (replaces sidebyside heatmap)
2. HFD vs Hu eWAT dot plot (replaces barplot)

SLURM: --partition=batch --mem=32G --cpus-per-task=4 --time=00:30:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
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

EXCLUDE = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

# AUCell threshold for "active" (cell is considered to have regulon active)
AUC_THRESH = 0.01

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

print("=" * 60)
print("42g_regulon_dotplots.py")
print("=" * 60)

# =============================================================================
# 1. Load all AUCell data
# =============================================================================
print("\nStep 1: Loading data...")

# Qiu per-cluster
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

# Hu per-depot
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
# 2. Compute mean AUCell and % active per group for all regulons
# =============================================================================
print("\nStep 2: Computing per-group statistics...")

all_groups = EOS_CLUSTERS + ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]
all_regulons = sorted(qiu_regulons | hu_regulons)

mean_scores = {}   # {(regulon, group): mean_auc}
pct_active = {}    # {(regulon, group): fraction_active}

for grp in EOS_CLUSTERS:
    mask = qiu_group == grp
    cells = qiu_auc.loc[mask]
    n_cells = len(cells)
    for tf in all_regulons:
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
    n_cells = len(cells)
    for tf in all_regulons:
        if tf in hu_auc.columns:
            vals = cells[tf].values
            mean_scores[(tf, grp)] = np.mean(vals)
            pct_active[(tf, grp)] = np.mean(vals > AUC_THRESH) * 100
        else:
            mean_scores[(tf, grp)] = 0.0
            pct_active[(tf, grp)] = 0.0

print(f"  Computed stats for {len(all_regulons)} regulons x {len(all_groups)} groups")

# =============================================================================
# Helper: Seurat-style dot plot
# =============================================================================
def dotplot(regulons, groups, mean_scores, pct_active, title_str, outpath,
            cmap_name="Blues", vmax=None, figwidth=None):
    """
    regulons: list of TF names (y-axis, top to bottom)
    groups: list of group names (x-axis)
    """
    n_reg = len(regulons)
    n_grp = len(groups)

    # Build matrices
    mean_mat = np.zeros((n_reg, n_grp))
    pct_mat = np.zeros((n_reg, n_grp))
    for i, tf in enumerate(regulons):
        for j, grp in enumerate(groups):
            mean_mat[i, j] = mean_scores.get((tf, grp), 0)
            pct_mat[i, j] = pct_active.get((tf, grp), 0)

    if vmax is None:
        vmax = max(0.05, np.max(mean_mat))

    # Color map
    cmap = LinearSegmentedColormap.from_list("custom",
        ["#f7f7f7", "#d1e5f0", "#6baed6", "#2171b5", "#08306b"])
    norm = Normalize(vmin=0, vmax=vmax)

    # Dot sizing: scale % active to dot area
    max_dot_size = 250
    min_dot_size = 10
    # Normalize pct to [0, 1] for sizing
    pct_max = max(1, np.max(pct_mat))

    if figwidth is None:
        figwidth = max(6, n_grp * 0.9 + 3)
    figheight = max(4, n_reg * 0.38 + 2)

    fig, ax = plt.subplots(figsize=(figwidth, figheight))

    for i in range(n_reg):
        for j in range(n_grp):
            m = mean_mat[i, j]
            p = pct_mat[i, j]
            size = min_dot_size + (p / pct_max) * (max_dot_size - min_dot_size)
            color = cmap(norm(m))
            ax.scatter(j, n_reg - 1 - i, s=size, c=[color],
                       edgecolors="grey50", linewidths=0.5, zorder=3)

    # Axes
    ax.set_xticks(range(n_grp))
    # Format x labels: add line break for Hu labels
    xlabels = []
    for g in groups:
        if g.startswith("Hu_"):
            xlabels.append(g.replace("Hu_", "Hu\n"))
        else:
            xlabels.append(g)
    ax.set_xticklabels(xlabels, fontsize=10, fontweight="bold")
    ax.set_yticks(range(n_reg))
    ax.set_yticklabels(list(reversed(regulons)), fontsize=9, fontstyle="italic")

    # Separator line between Qiu and Hu
    n_qiu = len([g for g in groups if not g.startswith("Hu")])
    if n_qiu < n_grp:
        ax.axvline(x=n_qiu - 0.5, color="black", linewidth=1.5, linestyle="--")

    # Grid
    ax.set_xlim(-0.5, n_grp - 0.5)
    ax.set_ylim(-0.5, n_reg - 0.5)
    ax.set_axisbelow(True)
    ax.grid(True, axis="both", linewidth=0.3, color="grey85", zorder=0)

    ax.set_title(title_str, fontsize=13, fontweight="bold", pad=12)

    # Dataset labels at top
    ax.text((n_qiu - 1) / 2, n_reg + 0.3, "Qiu (HFD)", fontsize=11,
            fontweight="bold", ha="center", color="#E63946")
    hu_start = n_qiu
    hu_end = n_grp - 1
    ax.text((hu_start + hu_end) / 2, n_reg + 0.3, "Hu (Healthy)", fontsize=11,
            fontweight="bold", ha="center", color="#2A9D8F")

    # Color bar for mean AUCell
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.02, pad=0.08, aspect=20)
    cbar.set_label("Mean AUCell", fontsize=10)

    # Size legend
    legend_pcts = [10, 25, 50, 75, 100]
    legend_pcts = [p for p in legend_pcts if p <= pct_max + 5]
    if not legend_pcts:
        legend_pcts = [10, 25, 50]

    legend_elements = []
    for p in legend_pcts:
        size = min_dot_size + (p / pct_max) * (max_dot_size - min_dot_size)
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
                   markersize=np.sqrt(size) * 0.7, label=f'{p:.0f}%',
                   markeredgecolor='grey50', markeredgewidth=0.5)
        )

    leg = ax.legend(handles=legend_elements, title="% Active",
                    loc="upper left", bbox_to_anchor=(1.15, 0.5),
                    frameon=False, fontsize=8, title_fontsize=9,
                    handletextpad=0.5)

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {os.path.basename(outpath)}")

# =============================================================================
# 3. Dot plot 1: HFD-enriched regulons (replaces sidebyside heatmap)
#    Top Qiu-only + top shared HFD-enriched
# =============================================================================
print("\nStep 3: HFD-enriched dot plot...")

# A) Top Qiu-only regulons
qiu_only_filtered = [t for t in qiu_only if t not in EXCLUDE]
qiu_only_max = {}
for tf in qiu_only_filtered:
    vals = [mean_scores.get((tf, cl), 0) for cl in EOS_CLUSTERS]
    qiu_only_max[tf] = max(vals)
top_qiu_only = sorted(qiu_only_max, key=qiu_only_max.get, reverse=True)[:10]

# B) Top shared regulons enriched in Qiu vs Hu
shared_filtered = [t for t in shared if t not in EXCLUDE]
enrichment = {}
for tf in shared_filtered:
    qiu_max_val = max(mean_scores.get((tf, cl), 0) for cl in EOS_CLUSTERS)
    hu_max_val = max(mean_scores.get((tf, f"Hu_{d}"), 0) for d in HU_DEPOTS)
    if qiu_max_val > 0.03:
        enrichment[tf] = qiu_max_val / (hu_max_val + 0.001)
top_shared_enriched = sorted(enrichment, key=enrichment.get, reverse=True)[:10]

hfd_regulons = top_qiu_only + top_shared_enriched

print(f"  HFD-only regulons: {top_qiu_only}")
print(f"  Shared HFD-enriched: {top_shared_enriched}")

dotplot(
    regulons=hfd_regulons,
    groups=all_groups,
    mean_scores=mean_scores,
    pct_active=pct_active,
    title_str="HFD-Enriched Regulons: Qiu (HFD) vs Hu (Healthy)",
    outpath=os.path.join(FIG_DIR, "42g_hfd_enriched_dotplot.pdf"),
    figwidth=10
)

# =============================================================================
# 4. Dot plot 2: Hu-enriched regulons (healthy-specific)
# =============================================================================
print("\nStep 4: Healthy-enriched dot plot...")

hu_only_filtered = [t for t in hu_only if t not in EXCLUDE]
hu_only_max = {}
for tf in hu_only_filtered:
    vals = [mean_scores.get((tf, f"Hu_{d}"), 0) for d in HU_DEPOTS]
    hu_only_max[tf] = max(vals)
top_hu_only = sorted(hu_only_max, key=hu_only_max.get, reverse=True)[:10]

# Shared regulons enriched in Hu
hu_enrichment = {}
for tf in shared_filtered:
    hu_max_val = max(mean_scores.get((tf, f"Hu_{d}"), 0) for d in HU_DEPOTS)
    qiu_max_val = max(mean_scores.get((tf, cl), 0) for cl in EOS_CLUSTERS)
    if hu_max_val > 0.03:
        hu_enrichment[tf] = hu_max_val / (qiu_max_val + 0.001)
top_shared_hu = sorted(hu_enrichment, key=hu_enrichment.get, reverse=True)[:10]

healthy_regulons = top_hu_only + top_shared_hu

print(f"  Hu-only regulons: {top_hu_only}")
print(f"  Shared Hu-enriched: {top_shared_hu}")

dotplot(
    regulons=healthy_regulons,
    groups=all_groups,
    mean_scores=mean_scores,
    pct_active=pct_active,
    title_str="Healthy-Enriched Regulons: Hu (Healthy) vs Qiu (HFD)",
    outpath=os.path.join(FIG_DIR, "42g_healthy_enriched_dotplot.pdf"),
    figwidth=10
)

# =============================================================================
# 5. Dot plot 3: Combined — top HFD + top Healthy in one figure
# =============================================================================
print("\nStep 5: Combined dot plot (HFD + Healthy)...")

# Take top 8 from each category
combined_regulons = (
    top_qiu_only[:5] +
    top_shared_enriched[:5] +
    top_hu_only[:5] +
    top_shared_hu[:5]
)
# Remove any duplicates while preserving order
seen = set()
combined_unique = []
for tf in combined_regulons:
    if tf not in seen:
        combined_unique.append(tf)
        seen.add(tf)

# Add category labels
n_hfd_only = min(5, len(top_qiu_only))
n_hfd_shared = min(5, len(top_shared_enriched))
n_hu_unique = min(5, len(top_hu_only))

dotplot(
    regulons=combined_unique,
    groups=all_groups,
    mean_scores=mean_scores,
    pct_active=pct_active,
    title_str="Regulon Comparison: HFD vs Healthy Adipose Eosinophils",
    outpath=os.path.join(FIG_DIR, "42g_combined_dotplot.pdf"),
    figwidth=10
)

# =============================================================================
# 6. Dot plot 4: Focused Qiu vs Hu_eWAT only
# =============================================================================
print("\nStep 6: Qiu vs Hu_eWAT focused dot plot...")

# Top regulons by absolute difference in mean AUCell (Qiu max vs Hu_eWAT)
all_filtered = [t for t in sorted(qiu_regulons | hu_regulons) if t not in EXCLUDE]
diff_data = []
for tf in all_filtered:
    qiu_max_val = max(mean_scores.get((tf, cl), 0) for cl in EOS_CLUSTERS)
    hu_ewat_val = mean_scores.get((tf, "Hu_eWAT"), 0)
    diff = qiu_max_val - hu_ewat_val
    if qiu_max_val > 0.02 or hu_ewat_val > 0.02:
        diff_data.append({"tf": tf, "qiu_max": qiu_max_val,
                         "hu_ewat": hu_ewat_val, "diff": diff})

diff_df = pd.DataFrame(diff_data).sort_values("diff", ascending=False)

# Top 12 HFD-enriched + top 8 Hu_eWAT-enriched
top_hfd_ewat = diff_df.head(12)["tf"].tolist()
top_hu_ewat = diff_df.tail(8)["tf"].tolist()[::-1]

ewat_regulons = top_hfd_ewat + top_hu_ewat

# Just show C0-C5 and Hu_eWAT for focused comparison
ewat_groups = EOS_CLUSTERS + ["Hu_eWAT"]

dotplot(
    regulons=ewat_regulons,
    groups=ewat_groups,
    mean_scores=mean_scores,
    pct_active=pct_active,
    title_str="Regulon Activity: HFD Clusters vs Healthy eWAT",
    outpath=os.path.join(FIG_DIR, "42g_qiu_vs_hu_ewat_dotplot.pdf"),
    figwidth=8
)

# =============================================================================
# 7. Summary
# =============================================================================
print("\n" + "=" * 60)
print("Summary of outputs:")
print("=" * 60)
print("  42g_hfd_enriched_dotplot.pdf    — HFD-enriched regulons (Qiu-only + shared enriched)")
print("  42g_healthy_enriched_dotplot.pdf — Healthy-enriched regulons (Hu-only + shared enriched)")
print("  42g_combined_dotplot.pdf         — Combined: top HFD + top Healthy regulons")
print("  42g_qiu_vs_hu_ewat_dotplot.pdf   — Focused: C0-C5 vs Hu_eWAT")
print(f"\n  AUCell threshold for % active: {AUC_THRESH}")

print("\n42g complete.")
