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

# This script's Venn/stacked-bar/heatmap panels are fundamentally a
# presence/absence question ("which regulons are dataset-specific"), so
# regulon SETS must come from the UNION of what each cluster's/depot's
# independently fit SCENIC network inferred, not the intersection -- taking
# the intersection (as this script previously did) silently redefines
# "Qiu-only" as "present in EVERY Qiu cluster AND absent from EVERY Hu
# depot," excluding any regulon inferred in only some clusters/depots. Raw
# per-cluster/per-depot DataFrames are kept in dicts (not concatenated) so a
# regulon's value can still be looked up wherever it was actually scored;
# `regulon_value()` below returns NaN, not a fabricated 0, wherever it
# wasn't.

# Qiu per-cluster
qiu_parts = {}
for cl in EOS_CLUSTERS:
    p = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(p):
        df = pd.read_csv(p, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T
        qiu_parts[cl] = df
    else:
        raise FileNotFoundError(
            f"Missing auc_matrix.csv for cluster {cl}. Regenerate the "
            f"missing cluster's SCENIC output rather than silently "
            f"continuing with an incomplete analysis.")

qiu_regulons = set()
for df in qiu_parts.values():
    qiu_regulons |= set(df.columns)

# Hu per-depot
hu_parts = {}
for depot in HU_DEPOTS:
    p = os.path.join(HU_SCENIC, depot, "auc_matrix.csv")
    if os.path.exists(p):
        df = pd.read_csv(p, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T
        hu_parts[depot] = df
    else:
        raise FileNotFoundError(
            f"Missing auc_matrix.csv for Hu depot {depot}. One missing "
            f"depot changes the regulon union/intersection and every "
            f"downstream comparison without any other error -- regenerate "
            f"it, or explicitly support a declared depot subset instead of "
            f"silently continuing.")

hu_regulons = set()
for df in hu_parts.values():
    hu_regulons |= set(df.columns)


def regulon_value(parts, key, tf):
    """Mean AUCell for `tf` in group `key` (a cluster or depot name), or NaN
    if that group's own SCENIC network never inferred `tf`."""
    df = parts[key]
    return df[tf].mean() if tf in df.columns else np.nan


shared = sorted(qiu_regulons & hu_regulons)
qiu_only = sorted(qiu_regulons - hu_regulons)
hu_only = sorted(hu_regulons - qiu_regulons)

n_shared = len(shared)
n_qiu_only = len(qiu_only)
n_hu_only = len(hu_only)

print(f"  Qiu: {len(qiu_regulons)} regulons (inferred in >=1 cluster)")
print(f"  Hu: {len(hu_regulons)} regulons (inferred in >=1 depot)")
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

all_columns = EOS_CLUSTERS + ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]

# A) Qiu-only regulons (top 10, excluding housekeeping). Values are NaN
# (not 0) in any cluster that didn't score a given qiu_only regulon --
# union-based selection means not every qiu_only regulon was scored in
# every cluster.
qiu_only_filtered = [t for t in qiu_only if t not in EXCLUDE]
qiu_only_means = {
    cl: {tf: regulon_value(qiu_parts, cl, tf) for tf in qiu_only_filtered}
    for cl in EOS_CLUSTERS
}
qiu_only_df = pd.DataFrame(qiu_only_means)
qiu_only_df["max_auc"] = qiu_only_df.max(axis=1, skipna=True)
top_qiu_only = qiu_only_df.sort_values("max_auc", ascending=False).head(10).index.tolist()

# B) Shared regulons enriched in Qiu (higher in Qiu than Hu). shared is
# union-based, so a given shared regulon is guaranteed a real value in at
# least one Qiu cluster and at least one Hu depot -- but not necessarily
# every one, so NaN entries are filtered out before taking the max rather
# than being read as a real 0.
shared_filtered = [t for t in shared if t not in EXCLUDE]
qiu_mean_shared = {}
hu_mean_shared = {}
for tf in shared_filtered:
    qiu_vals = [v for v in (regulon_value(qiu_parts, cl, tf) for cl in EOS_CLUSTERS)
               if not np.isnan(v)]
    qiu_mean_shared[tf] = max(qiu_vals) if qiu_vals else np.nan

    hu_vals = [v for v in (regulon_value(hu_parts, depot, tf) for depot in HU_DEPOTS)
              if not np.isnan(v)]
    hu_mean_shared[tf] = max(hu_vals) if hu_vals else np.nan

# Fold enrichment in Qiu
enrichment = {}
for tf in shared_filtered:
    qv, hv = qiu_mean_shared[tf], hu_mean_shared[tf]
    if hv > 0.005:
        enrichment[tf] = qv / hv
    elif qv > 0.02:
        enrichment[tf] = qv / 0.005

top_shared_enriched = sorted(enrichment, key=enrichment.get, reverse=True)[:10]

# Combine: top Qiu-only + top shared-enriched
hfd_regulons = top_qiu_only + top_shared_enriched

# Build heatmap matrix. NaN cells (regulon not scored in that cluster/depot)
# are rendered blank below, not as a fabricated zero-activity value.
heat_data = []
for tf in hfd_regulons:
    row = {}
    for cl in EOS_CLUSTERS:
        row[cl] = regulon_value(qiu_parts, cl, tf)
    for depot in HU_DEPOTS:
        row[f"Hu_{depot}"] = regulon_value(hu_parts, depot, tf)
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
               vmax=max(0.2, np.nanmax(heat_df.values)))

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

# Values in cells (NaN = regulon not inferred in that cluster/depot -- leave blank)
for i in range(len(hfd_regulons)):
    for j in range(len(all_columns)):
        val = heat_df.values[i, j]
        if np.isnan(val):
            continue
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

# Top 15 regulons by max Qiu AUC (excluding housekeeping). qiu_regulons is
# the union across clusters, so a regulon here is guaranteed a real value in
# at least one cluster, just not necessarily every one.
all_tfs = [t for t in sorted(qiu_regulons) if t not in EXCLUDE]
qiu_max = {}
hu_ewat_val = {}
tfs_with_hu_ewat_data = []
for tf in all_tfs:
    vals = [v for v in (regulon_value(qiu_parts, cl, tf) for cl in EOS_CLUSTERS)
           if not np.isnan(v)]
    qiu_max[tf] = max(vals)

    # Hu eWAT mean -- only if this regulon was actually scored in the eWAT
    # depot's own SCENIC network. A regulon not scored there was "not
    # inferred," not measured and found at AUCell 0; defaulting it to 0
    # would fabricate a huge log2FC below (qv / ~0) from missing data
    # rather than a real Qiu-vs-Hu difference, so such regulons are
    # excluded from this fold-change comparison entirely.
    val = regulon_value(hu_parts, "eWAT", tf) if "eWAT" in hu_parts else np.nan
    if not np.isnan(val):
        hu_ewat_val[tf] = val
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
