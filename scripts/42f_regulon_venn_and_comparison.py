#!/usr/bin/env python3
"""
42f_regulon_venn_and_comparison.py
1. Venn diagram: Shared / Qiu-only / Hu-only regulon counts
2. Side-by-side heatmap: top HFD-enriched regulons showing BOTH
   Qiu C0-C5 AND Hu eWAT/scWAT/BAT (so the contrast is visible)
3. Focused comparison: HFD-enriched regulons vs Hu eWAT specifically

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
from matplotlib.patches import Circle
from matplotlib.colors import LinearSegmentedColormap
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

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

print("=" * 60)
print("42f_regulon_venn_and_comparison.py")
print("=" * 60)

# =============================================================================
# 1. Load all AUCell data (same as 42e)
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
    else:
        raise FileNotFoundError(
            f"Missing auc_matrix.csv for cluster {cl}. Regenerate the "
            f"missing cluster's SCENIC output rather than silently "
            f"continuing with an incomplete analysis.")

# Inner join, not outer + fillna(0): SCENIC is fit independently per
# cluster, so a regulon absent from one cluster's network is "not inferred,"
# not a real AUCell score of zero. Keep only regulons scored in every
# cluster so qiu_regulons genuinely reflects Qiu-wide activity.
qiu_auc = pd.concat(qiu_parts, join="inner")
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

# Same reasoning as qiu_auc above: inner join across independently fitted
# per-depot networks.
hu_auc = pd.concat(hu_parts, join="inner")
hu_group = pd.concat(hu_groups)
hu_regulons = set(hu_auc.columns)

shared = sorted(qiu_regulons & hu_regulons)
qiu_only = sorted(qiu_regulons - hu_regulons)
hu_only = sorted(hu_regulons - qiu_regulons)

n_shared = len(shared)
n_qiu_only = len(qiu_only)
n_hu_only = len(hu_only)

print(f"  Qiu: {len(qiu_regulons)} regulons")
print(f"  Hu: {len(hu_regulons)} regulons")
print(f"  Shared: {n_shared}, Qiu-only: {n_qiu_only}, Hu-only: {n_hu_only}")

# =============================================================================
# 2. Venn diagram
# =============================================================================
print("\nStep 2: Venn diagram...")

fig, ax = plt.subplots(figsize=(7, 5))

# Two-circle Venn
c1_x, c2_x = -0.8, 0.8
c_y = 0
r = 1.8

circle1 = Circle((c1_x, c_y), r, facecolor="#E63946", alpha=0.25,
                  edgecolor="#E63946", linewidth=2.5)
circle2 = Circle((c2_x, c_y), r, facecolor="#2A9D8F", alpha=0.25,
                  edgecolor="#2A9D8F", linewidth=2.5)
ax.add_patch(circle1)
ax.add_patch(circle2)

# Numbers
ax.text(-1.5, 0, str(n_qiu_only), fontsize=28, fontweight="bold",
        ha="center", va="center", color="#E63946")
ax.text(0, 0, str(n_shared), fontsize=28, fontweight="bold",
        ha="center", va="center", color="#333333")
ax.text(1.5, 0, str(n_hu_only), fontsize=28, fontweight="bold",
        ha="center", va="center", color="#2A9D8F")

# Labels
ax.text(-1.5, -1.2, "HFD-specific", fontsize=12, ha="center", va="center",
        color="#E63946", fontstyle="italic")
ax.text(0, -1.2, "Shared", fontsize=12, ha="center", va="center",
        color="#333333", fontstyle="italic")
ax.text(1.5, -1.2, "Healthy-specific", fontsize=12, ha="center", va="center",
        color="#2A9D8F", fontstyle="italic")

# Dataset labels
ax.text(-2.2, 1.8, f"Qiu et al.\n(HFD adipose)\n{len(qiu_regulons)} total",
        fontsize=11, ha="center", va="center", fontweight="bold", color="#E63946")
ax.text(2.2, 1.8, f"Hu et al.\n(Healthy adipose)\n{len(hu_regulons)} total",
        fontsize=11, ha="center", va="center", fontweight="bold", color="#2A9D8F")

ax.set_xlim(-3.5, 3.5)
ax.set_ylim(-2.5, 3.0)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("Regulon Programs: HFD vs Healthy Adipose Eosinophils",
             fontsize=14, fontweight="bold", pad=10)

fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "42f_regulon_venn.pdf"),
            dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved: 42f_regulon_venn.pdf")

# =============================================================================
# 3. Stacked bar: percentage breakdown
# =============================================================================
print("\nStep 3: Stacked bar chart...")

fig, ax = plt.subplots(figsize=(5, 3))

categories = ["Qiu (HFD)", "Hu (Healthy)"]
shared_counts = [n_shared, n_shared]
unique_counts = [n_qiu_only, n_hu_only]
totals = [len(qiu_regulons), len(hu_regulons)]

shared_pct = [n_shared / t * 100 for t in totals]
unique_pct = [u / t * 100 for u, t in zip(unique_counts, totals)]

bars1 = ax.barh(categories, shared_pct, color="#999999", label="Shared", height=0.5)
bars2 = ax.barh(categories, unique_pct, left=shared_pct,
                color=["#E63946", "#2A9D8F"], height=0.5)

# Add count labels
for i, (sp, up, sc, uc) in enumerate(zip(shared_pct, unique_pct, shared_counts, unique_counts)):
    ax.text(sp / 2, i, f"{sc}\n({sp:.0f}%)", ha="center", va="center",
            fontsize=9, fontweight="bold", color="white")
    ax.text(sp + up / 2, i, f"{uc}\n({up:.0f}%)", ha="center", va="center",
            fontsize=9, fontweight="bold", color="white")

ax.set_xlabel("Percentage of regulons", fontsize=11)
ax.set_xlim(0, 100)
ax.legend(["Shared", "Dataset-specific"], loc="lower right", fontsize=9)
ax.set_title("Regulon Program Overlap", fontsize=12, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "42f_regulon_stacked_bar.pdf"),
            dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved: 42f_regulon_stacked_bar.pdf")

# =============================================================================
# 4. Side-by-side heatmap: HFD-enriched regulons showing BOTH datasets
# =============================================================================
print("\nStep 4: Side-by-side heatmap (HFD-enriched, both datasets)...")

# Compute means for all groups
all_columns = EOS_CLUSTERS + ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]
group_means = {}

for cl in EOS_CLUSTERS:
    mask = qiu_group == cl
    group_means[cl] = qiu_auc.loc[mask].mean()

for depot in HU_DEPOTS:
    mask = hu_group == f"Hu_{depot}"
    group_means[f"Hu_{depot}"] = hu_auc.loc[mask].mean()

# For Qiu-only regulons, Hu values are 0 (absent)
# For shared regulons, find those enriched in Qiu vs Hu

# A) Qiu-only regulons (top 10, excluding housekeeping)
qiu_only_filtered = [t for t in qiu_only if t not in EXCLUDE]
qiu_only_means = {}
for cl in EOS_CLUSTERS:
    mask = qiu_group == cl
    qiu_only_means[cl] = qiu_auc.loc[mask, [t for t in qiu_only_filtered if t in qiu_auc.columns]].mean()
qiu_only_df = pd.DataFrame(qiu_only_means)
qiu_only_df["max_auc"] = qiu_only_df.max(axis=1)
top_qiu_only = qiu_only_df.sort_values("max_auc", ascending=False).head(10).index.tolist()

# B) Shared regulons enriched in Qiu (higher in Qiu than Hu)
shared_filtered = [t for t in shared if t not in EXCLUDE]
qiu_mean_shared = {}
hu_mean_shared = {}
for tf in shared_filtered:
    qiu_vals = []
    for cl in EOS_CLUSTERS:
        mask = qiu_group == cl
        if tf in qiu_auc.columns:
            qiu_vals.append(qiu_auc.loc[mask, tf].mean())
    qiu_mean_shared[tf] = max(qiu_vals) if qiu_vals else 0

    hu_vals = []
    for depot in HU_DEPOTS:
        mask = hu_group == f"Hu_{depot}"
        if tf in hu_auc.columns:
            hu_vals.append(hu_auc.loc[mask, tf].mean())
    hu_mean_shared[tf] = max(hu_vals) if hu_vals else 0

# Fold enrichment in Qiu
enrichment = {}
for tf in shared_filtered:
    if hu_mean_shared[tf] > 0.005:
        enrichment[tf] = qiu_mean_shared[tf] / hu_mean_shared[tf]
    elif qiu_mean_shared[tf] > 0.02:
        enrichment[tf] = qiu_mean_shared[tf] / 0.005

top_shared_enriched = sorted(enrichment, key=enrichment.get, reverse=True)[:10]

# Combine: top Qiu-only + top shared-enriched
hfd_regulons = top_qiu_only + top_shared_enriched

# Build heatmap matrix
heat_data = []
for tf in hfd_regulons:
    row = {}
    for cl in EOS_CLUSTERS:
        mask = qiu_group == cl
        row[cl] = qiu_auc.loc[mask, tf].mean() if tf in qiu_auc.columns else 0
    for depot in HU_DEPOTS:
        mask = hu_group == f"Hu_{depot}"
        row[f"Hu_{depot}"] = hu_auc.loc[mask, tf].mean() if tf in hu_auc.columns else 0
    heat_data.append(row)

heat_df = pd.DataFrame(heat_data, index=hfd_regulons)
heat_df = heat_df[all_columns]

# Categorize regulons
categories = ["HFD-only"] * len(top_qiu_only) + ["Shared\n(HFD-enriched)"] * len(top_shared_enriched)

# Plot
fig, ax = plt.subplots(figsize=(8, 8))

cmap = LinearSegmentedColormap.from_list("white_blue_dark",
    ["#f7f7f7", "#d1e5f0", "#6baed6", "#2171b5", "#08306b"])

im = ax.imshow(heat_df.values, cmap=cmap, aspect="auto", vmin=0,
               vmax=max(0.2, heat_df.values.max()))

# Grid lines
ax.set_xticks(range(len(all_columns)))
ax.set_xticklabels(["C0", "C1", "C2", "C3", "C4", "C5",
                     "Hu\neWAT", "Hu\nscWAT", "Hu\nBAT"],
                    fontsize=10, fontweight="bold")
ax.set_yticks(range(len(hfd_regulons)))
ax.set_yticklabels(hfd_regulons, fontsize=9, fontstyle="italic")

# Separator line between Qiu and Hu
ax.axvline(x=5.5, color="black", linewidth=2.5, linestyle="-")

# Separator line between Qiu-only and shared-enriched
ax.axhline(y=len(top_qiu_only) - 0.5, color="black", linewidth=1.5, linestyle="--")

# Category labels on right
ax.text(len(all_columns) + 0.3, len(top_qiu_only) / 2 - 0.5, "HFD-only",
        fontsize=10, fontweight="bold", color="#E63946", va="center",
        rotation=270)
ax.text(len(all_columns) + 0.3, len(top_qiu_only) + len(top_shared_enriched) / 2 - 0.5,
        "Shared\n(HFD-enriched)", fontsize=9, fontweight="bold", color="#555555",
        va="center", rotation=270)

# Values in cells
for i in range(len(hfd_regulons)):
    for j in range(len(all_columns)):
        val = heat_df.values[i, j]
        color = "white" if val > 0.06 else "black"
        ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                fontsize=7, color=color)

# Dataset labels at top
ax.text(2.5, -1.5, "Qiu (HFD)", fontsize=12, fontweight="bold",
        ha="center", color="#E63946")
ax.text(7, -1.5, "Hu (Healthy)", fontsize=12, fontweight="bold",
        ha="center", color="#2A9D8F")

cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.08, label="Mean AUCell")

ax.set_title("HFD-Enriched Regulon Programs\nQiu (HFD) vs Hu (Healthy Adipose)",
             fontsize=13, fontweight="bold", pad=25)

fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "42f_hfd_enriched_sidebyside.pdf"),
            dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved: 42f_hfd_enriched_sidebyside.pdf")

# =============================================================================
# 5. Focused: HFD-enriched vs Hu eWAT specifically (bar chart)
# =============================================================================
print("\nStep 5: HFD vs Hu_eWAT bar comparison...")

# Top 15 regulons by max Qiu AUC (excluding housekeeping)
all_tfs = [t for t in sorted(qiu_auc.columns) if t not in EXCLUDE]
qiu_max = {}
hu_ewat_val = {}
tfs_with_hu_ewat_data = []
for tf in all_tfs:
    # Max across Qiu clusters (tf is drawn from qiu_auc.columns, so always present)
    vals = [qiu_auc.loc[qiu_group == cl, tf].mean() for cl in EOS_CLUSTERS]
    qiu_max[tf] = max(vals)

    # Hu eWAT mean -- only if this regulon was actually scored in Hu's
    # SCENIC network at all. A regulon not in hu_auc.columns was "not
    # inferred" in Hu, not measured and found at AUCell 0; defaulting it to
    # 0 here would fabricate a huge log2FC below (qv / ~0) from missing
    # data rather than a real Qiu-vs-Hu difference, so such regulons are
    # excluded from this fold-change comparison entirely.
    if tf in hu_auc.columns:
        mask = hu_group == "Hu_eWAT"
        hu_ewat_val[tf] = hu_auc.loc[mask, tf].mean()
        tfs_with_hu_ewat_data.append(tf)

# Sort by fold-change (Qiu / Hu_eWAT), restricted to regulons actually
# scored in both datasets
fc_data = []
for tf in tfs_with_hu_ewat_data:
    qv = qiu_max[tf]
    hv = hu_ewat_val[tf]
    if qv > 0.03:  # only meaningful regulons
        fc = np.log2((qv + 0.001) / (hv + 0.001))
        fc_data.append({"TF": tf, "Qiu_max": qv, "Hu_eWAT": hv, "log2FC": fc})

fc_df = pd.DataFrame(fc_data).sort_values("log2FC", ascending=False)

# Top 15 HFD-enriched + top 10 Hu_eWAT-enriched
top_hfd = fc_df.head(15)
top_hu = fc_df.tail(10).sort_values("log2FC", ascending=True)
plot_df = pd.concat([top_hfd, top_hu])

fig, ax = plt.subplots(figsize=(7, 8))
colors = ["#E63946" if x > 0 else "#2A9D8F" for x in plot_df["log2FC"]]
bars = ax.barh(range(len(plot_df)), plot_df["log2FC"].values, color=colors,
               edgecolor="white", linewidth=0.5, height=0.7)
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels(plot_df["TF"].values, fontsize=9, fontstyle="italic")
ax.axvline(x=0, color="black", linewidth=1)
ax.set_xlabel("log₂(Qiu max / Hu eWAT)", fontsize=11)
ax.set_title("Regulon Enrichment: HFD vs Healthy eWAT",
             fontsize=13, fontweight="bold")
ax.invert_yaxis()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#E63946", label="HFD-enriched"),
    Patch(facecolor="#2A9D8F", label="Healthy eWAT-enriched")
]
ax.legend(handles=legend_elements, loc="lower right", fontsize=10, frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "42f_hfd_vs_hu_ewat_barplot.pdf"),
            dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved: 42f_hfd_vs_hu_ewat_barplot.pdf")

# =============================================================================
# 6. Summary
# =============================================================================
print("\n" + "=" * 60)
print("Summary of outputs:")
print("=" * 60)
print(f"  42f_regulon_venn.pdf         — Venn diagram: {n_shared} shared / {n_qiu_only} HFD-only / {n_hu_only} healthy-only")
print(f"  42f_regulon_stacked_bar.pdf  — Percentage breakdown by dataset")
print(f"  42f_hfd_enriched_sidebyside.pdf — Side-by-side heatmap: Qiu C0-C5 vs Hu depots")
print(f"  42f_hfd_vs_hu_ewat_barplot.pdf  — log2FC bar chart: HFD vs Hu eWAT")

print(f"\nTop HFD-enriched vs Hu eWAT:")
for _, row in top_hfd.head(10).iterrows():
    print(f"  {row['TF']}: Qiu={row['Qiu_max']:.3f}, Hu_eWAT={row['Hu_eWAT']:.3f}, log2FC={row['log2FC']:.1f}")

print(f"\nTop Hu eWAT-enriched vs HFD:")
for _, row in top_hu.head(5).iterrows():
    print(f"  {row['TF']}: Qiu={row['Qiu_max']:.3f}, Hu_eWAT={row['Hu_eWAT']:.3f}, log2FC={row['log2FC']:.1f}")

print("\n42f complete.")
