#!/usr/bin/env python3
"""
42e_regulon_unique_analysis.py
Identify regulons unique to Hu (healthy) vs Qiu (HFD) datasets.
Also show top 5 per cluster instead of top 3.

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
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
HU_SCENIC  = os.path.join(BASE, "scenic", "hu_scenic")
HU_META    = os.path.join(BASE, "scenic", "hu_depot_metadata.csv")
FIG_DIR    = os.path.join(BASE, "scenic", "figures")

EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]
HU_DEPOTS = ["eWAT", "scWAT", "BAT"]

EXCLUDE = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"])

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

print("=" * 60)
print("42e_regulon_unique_analysis.py")
print("=" * 60)

# =============================================================================
# 1. Load all AUCell data
# =============================================================================
print("\nStep 1: Loading data...")

# This script's purpose is presence/absence (which regulons exist only in
# one dataset), not a magnitude comparison -- so regulon SETS must be built
# from the UNION of what each cluster's/depot's independently fit SCENIC
# network inferred, not the intersection. Taking the intersection (as this
# script previously did) silently redefines "Qiu-only" as "present in EVERY
# Qiu cluster AND absent from EVERY Hu depot," which excludes any regulon
# inferred in only one or two Qiu clusters -- not the normal meaning of
# "dataset-specific." Per-cluster/per-depot DataFrames are kept (not
# concatenated) so each regulon's value can still be looked up wherever it
# actually was scored; a cluster/depot that didn't score a given regulon
# contributes NaN there, not a fabricated 0.

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

qiu_union_regulons = set()
for df in qiu_parts.values():
    qiu_union_regulons |= set(df.columns)
print(f"  Qiu: {sum(len(df) for df in qiu_parts.values())} cells, "
      f"{len(qiu_union_regulons)} regulons inferred in >=1 cluster")

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

hu_union_regulons = set()
for df in hu_parts.values():
    hu_union_regulons |= set(df.columns)
print(f"  Hu: {sum(len(df) for df in hu_parts.values())} cells, "
      f"{len(hu_union_regulons)} regulons inferred in >=1 depot")

shared = sorted(qiu_union_regulons & hu_union_regulons)
qiu_only = sorted(qiu_union_regulons - hu_union_regulons)
hu_only = sorted(hu_union_regulons - qiu_union_regulons)

print(f"\n  Shared: {len(shared)}")
print(f"  Qiu-only (HFD-specific): {len(qiu_only)}")
print(f"  Hu-only (healthy-specific): {len(hu_only)}")

# =============================================================================
# 2. Hu-only regulons: active in healthy but absent in HFD
# =============================================================================
print("\n" + "=" * 60)
print("Step 2: Hu-only regulons (healthy adipose, absent in HFD)")
print("=" * 60)

hu_only_filtered = [t for t in hu_only if t not in EXCLUDE]
hu_only_means = {}
for depot in HU_DEPOTS:
    df = hu_parts[depot]
    # NaN, not 0, for a regulon this depot's own SCENIC network never
    # inferred -- every hu_only regulon is guaranteed a real value in at
    # least one depot (that's why it's in hu_only), just not necessarily
    # this one.
    hu_only_means[f"Hu_{depot}"] = {
        tf: (df[tf].mean() if tf in df.columns else np.nan) for tf in hu_only_filtered
    }

hu_only_df = pd.DataFrame(hu_only_means)
hu_only_df["max_auc"] = hu_only_df.max(axis=1, skipna=True)
hu_only_df["best_depot"] = hu_only_df[["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]].idxmax(axis=1, skipna=True)
hu_only_df = hu_only_df.sort_values("max_auc", ascending=False)

print(f"\nTop 20 Hu-only regulons (active in healthy, absent in HFD):")
print(hu_only_df.head(20).to_string())

# Save full list
hu_only_df.to_csv(os.path.join(FIG_DIR, "42e_hu_only_regulons.csv"))

# =============================================================================
# 3. Qiu-only regulons: HFD-specific, absent in healthy
# =============================================================================
print("\n" + "=" * 60)
print("Step 3: Qiu-only regulons (HFD-specific, absent in healthy)")
print("=" * 60)

qiu_only_filtered = [t for t in qiu_only if t not in EXCLUDE]
qiu_cluster_means = {}
for cl in EOS_CLUSTERS:
    df = qiu_parts[cl]
    # NaN, not 0, for a regulon this cluster's own SCENIC network never
    # inferred -- every qiu_only regulon is guaranteed a real value in at
    # least one cluster (that's why it's in qiu_only), just not necessarily
    # this one.
    qiu_cluster_means[cl] = {
        tf: (df[tf].mean() if tf in df.columns else np.nan) for tf in qiu_only_filtered
    }

qiu_only_df = pd.DataFrame(qiu_cluster_means)
qiu_only_df["max_auc"] = qiu_only_df.max(axis=1, skipna=True)
qiu_only_df["best_cluster"] = qiu_only_df[EOS_CLUSTERS].idxmax(axis=1, skipna=True)
qiu_only_df = qiu_only_df.sort_values("max_auc", ascending=False)

print(f"\nTop 20 Qiu-only regulons (HFD-specific, absent in healthy):")
print(qiu_only_df.head(20).to_string())

qiu_only_df.to_csv(os.path.join(FIG_DIR, "42e_qiu_only_regulons.csv"))

# =============================================================================
# 4. Heatmap of top Hu-only regulons
# =============================================================================
print("\nStep 4: Plotting Hu-only regulon heatmap...")

top_hu_only = hu_only_df.head(15).index.tolist()
if top_hu_only:
    # A hu_only regulon may not have been scored in every depot (union-based
    # selection) -- those cells are NaN, not a fabricated 0; imshow renders
    # them blank and the text loop below skips them.
    plot_mat = hu_only_df.loc[top_hu_only, ["Hu_eWAT", "Hu_scWAT", "Hu_BAT"]]

    fig, ax = plt.subplots(figsize=(5, 6))
    im = ax.imshow(plot_mat.values, cmap="Reds", aspect="auto", vmin=0)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["eWAT", "scWAT", "BAT"], fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(top_hu_only)))
    ax.set_yticklabels(top_hu_only, fontsize=9, fontstyle="italic")
    ax.set_title("Healthy-Specific Regulons\n(active in Hu, absent in HFD)",
                 fontsize=12, fontweight="bold")
    for i in range(len(top_hu_only)):
        for j in range(3):
            val = plot_mat.values[i, j]
            if np.isnan(val):
                continue  # regulon not inferred in this depot -- leave blank
            ax.text(j, i, f"{val:.3f}",
                    ha="center", va="center", fontsize=7,
                    color="white" if val > 0.05 else "black")
    fig.colorbar(im, ax=ax, fraction=0.03, label="AUCell")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "42e_hu_only_regulon_heatmap.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: 42e_hu_only_regulon_heatmap.pdf")

# =============================================================================
# 5. Heatmap of top Qiu-only regulons
# =============================================================================
print("\nStep 5: Plotting Qiu-only regulon heatmap...")

top_qiu_only = qiu_only_df.head(15).index.tolist()
if top_qiu_only:
    # A qiu_only regulon may not have been scored in every cluster
    # (union-based selection) -- those cells are NaN, not a fabricated 0;
    # imshow renders them blank and the text loop below skips them.
    plot_mat = qiu_only_df.loc[top_qiu_only, EOS_CLUSTERS]

    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(plot_mat.values, cmap="Blues", aspect="auto", vmin=0)
    ax.set_xticks(range(6))
    ax.set_xticklabels(EOS_CLUSTERS, fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(top_qiu_only)))
    ax.set_yticklabels(top_qiu_only, fontsize=9, fontstyle="italic")
    ax.set_title("HFD-Specific Regulons\n(active in Qiu, absent in healthy Hu)",
                 fontsize=12, fontweight="bold")
    for i in range(len(top_qiu_only)):
        for j in range(6):
            val = plot_mat.values[i, j]
            if np.isnan(val):
                continue  # regulon not inferred in this cluster -- leave blank
            ax.text(j, i, f"{val:.3f}",
                    ha="center", va="center", fontsize=7,
                    color="white" if val > 0.05 else "black")
    fig.colorbar(im, ax=ax, fraction=0.03, label="AUCell")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "42e_qiu_only_regulon_heatmap.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: 42e_qiu_only_regulon_heatmap.pdf")

# =============================================================================
# 6. Summary
# =============================================================================
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Total Qiu regulons: {len(qiu_regulons)}")
print(f"Total Hu regulons: {len(hu_regulons)}")
print(f"Shared: {len(shared)}")
print(f"Qiu-only (HFD-driven programs): {len(qiu_only)}")
print(f"Hu-only (healthy programs lost in HFD): {len(hu_only)}")

print(f"\nTop 5 Hu-only (suppressed by HFD):")
for i, tf in enumerate(hu_only_df.head(5).index):
    depot = hu_only_df.loc[tf, "best_depot"]
    auc = hu_only_df.loc[tf, "max_auc"]
    print(f"  {tf}: AUC={auc:.4f} in {depot}")

print(f"\nTop 5 Qiu-only (induced by HFD):")
for i, tf in enumerate(qiu_only_df.head(5).index):
    cl = qiu_only_df.loc[tf, "best_cluster"]
    auc = qiu_only_df.loc[tf, "max_auc"]
    print(f"  {tf}: AUC={auc:.4f} in {cl}")

print("\n42e complete.")
