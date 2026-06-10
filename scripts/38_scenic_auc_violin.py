#!/usr/bin/env python3
"""
38_scenic_auc_violin.py
SCENIC regulon AUC violin plots by condition and replicate.

Condition-resolved SCENIC regulon AUC analysis:
  - "How eosinophil regulons differ between conditions"
  - "violin plots showing distributions of AUCs across cells in each condition/replicate"

Reads per-cluster AUCell matrices, merges with cell metadata,
plots violin + strip plots of top regulons per cluster across conditions.

SLURM: --partition=batch --mem=32G --cpus-per-task=4 --time=01:00:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
META_FILE  = os.path.join(BASE, "scenic", "eos_cell_metadata.csv")
FIG_DIR    = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]

COND_ORDER = ["WT_Female_HFD", "KO_Female_HFD", "WT_Male_HFD", "KO_Male_HFD",
              "WT_Female_HFDInfected", "KO_Female_HFDInfected"]
COND_COLORS = {
    "WT_Female_HFD": "#B2182B", "KO_Female_HFD": "#E08080",
    "WT_Male_HFD": "#2166AC",   "KO_Male_HFD": "#6BAED6",
    "WT_Female_HFDInfected": "#D6604D", "KO_Female_HFDInfected": "#FC8D59",
}
COND_SHORT = {
    "WT_Female_HFD": "WT\u2640", "KO_Female_HFD": "KO\u2640",
    "WT_Male_HFD": "WT\u2642",   "KO_Male_HFD": "KO\u2642",
    "WT_Female_HFDInfected": "WT\u2640\nInf", "KO_Female_HFDInfected": "KO\u2640\nInf",
}

TOP_PER_CLUSTER = 3   # top regulons per cluster to plot
MIN_AUC         = 0.02

EXCLUDE = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

print("=" * 60)
print("38_scenic_auc_violin.py")
print("=" * 60)

# =============================================================================
# 1. Load per-cluster AUCell matrices
# =============================================================================
print("\nStep 1: Loading per-cluster AUCell data...")

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

auc_parts = []
for cl in EOS_CLUSTERS:
    auc_path = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(auc_path):
        df = pd.read_csv(auc_path, index_col=0)
        df.columns = [extract_tf(c) for c in df.columns]
        df = df.T.groupby(level=0).max().T  # deduplicate TFs
        auc_parts.append(df)
        print(f"  {cl}: {df.shape[0]} cells x {df.shape[1]} regulons")

auc = pd.concat(auc_parts, join="outer").fillna(0.0)
print(f"  Combined: {auc.shape[0]} cells x {auc.shape[1]} regulons")

# =============================================================================
# 2. Load metadata
# =============================================================================
print("\nStep 2: Loading metadata...")
meta = pd.read_csv(META_FILE)
if "barcode" in meta.columns:
    meta = meta.set_index("barcode")

common = auc.index.intersection(meta.index)
auc = auc.loc[common]
meta = meta.loc[common]
print(f"  Aligned: {len(common)} cells")

# Derive condition (strip trailing _N replicate)
if "condition" not in meta.columns:
    meta["condition"] = meta["orig.ident"].str.replace(r"_\d+$", "", regex=True)

# Keep only conditions in our order
valid = meta["condition"].isin(COND_ORDER)
auc = auc.loc[valid]
meta = meta.loc[valid]
print(f"  After condition filter: {len(meta)} cells")
print(f"  Conditions: {meta['condition'].unique().tolist()}")

# Cluster assignment
if "eos_cluster" in meta.columns:
    meta["cluster"] = meta["eos_cluster"]
elif "seurat_clusters" in meta.columns:
    meta["cluster"] = meta["seurat_clusters"]

# =============================================================================
# 3. Select top regulons per cluster
# =============================================================================
print("\nStep 3: Selecting top regulons per cluster...")

cluster_means = {}
for cl in EOS_CLUSTERS:
    mask = meta["cluster"] == cl
    if mask.sum() == 0:
        continue
    cluster_means[cl] = auc.loc[mask].mean()

mean_df = pd.DataFrame(cluster_means)

selected = {}
for cl in EOS_CLUSTERS:
    if cl not in mean_df.columns:
        continue
    scores = []
    for tf in auc.columns:
        if tf in EXCLUDE:
            continue
        max_val = mean_df.loc[tf, cl]
        if max_val < MIN_AUC:
            continue
        if cl != mean_df.columns[np.argmax(mean_df.loc[tf].values)]:
            continue
        sorted_vals = sorted(mean_df.loc[tf].values, reverse=True)
        uniqueness = sorted_vals[0] - sorted_vals[1] if len(sorted_vals) >= 2 else sorted_vals[0]
        scores.append((tf, max_val, uniqueness))

    scores.sort(key=lambda x: x[2], reverse=True)
    for tf, auc_val, uniq in scores[:TOP_PER_CLUSTER]:
        selected[tf] = cl

print(f"  Selected {len(selected)} regulons:")
for cl in EOS_CLUSTERS:
    tfs = [tf for tf, c in selected.items() if c == cl]
    if tfs:
        print(f"    {cl}: {', '.join(tfs)}")

# =============================================================================
# 4. Plot violin plots — one figure per cluster, showing its top regulons
# =============================================================================
print("\nStep 4: Generating violin plots...")

for cl in EOS_CLUSTERS:
    tfs = [tf for tf, c in selected.items() if c == cl]
    if not tfs:
        continue

    n_tfs = len(tfs)
    fig, axes = plt.subplots(1, n_tfs, figsize=(4.5 * n_tfs, 4.5))
    if n_tfs == 1:
        axes = [axes]

    for ax, tf in zip(axes, tfs):
        plot_df = pd.DataFrame({
            "AUC": auc[tf].values,
            "condition": meta["condition"].values,
            "sample": meta["orig.ident"].values if "orig.ident" in meta.columns else "unknown",
        })
        plot_df["cond_short"] = plot_df["condition"].map(COND_SHORT)
        plot_df["cond_short"] = pd.Categorical(
            plot_df["cond_short"],
            categories=[COND_SHORT[c] for c in COND_ORDER if c in COND_SHORT],
            ordered=True
        )

        palette = [COND_COLORS[c] for c in COND_ORDER if c in plot_df["condition"].values]
        order = [COND_SHORT[c] for c in COND_ORDER if c in plot_df["condition"].values]

        sns.violinplot(data=plot_df, x="cond_short", y="AUC", order=order,
                       palette=palette, cut=0, inner=None, linewidth=0.8,
                       alpha=0.7, ax=ax)

        # Overlay per-sample means as strip dots
        sample_means = plot_df.groupby(["cond_short", "sample"])["AUC"].mean().reset_index()
        sns.stripplot(data=sample_means, x="cond_short", y="AUC", order=order,
                      color="black", size=4, alpha=0.8, jitter=0.15, ax=ax)

        ax.set_title(tf, fontsize=13, fontweight="bold", fontstyle="italic")
        ax.set_xlabel("")
        ax.set_ylabel("AUCell Score" if ax == axes[0] else "")
        ax.tick_params(axis="x", rotation=0)

    fig.suptitle(f"{cl} — Top Regulon Activity by Condition",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()

    out_pdf = os.path.join(FIG_DIR, f"38_regulon_violin_{cl}.pdf")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_pdf}")

# =============================================================================
# 5. Combined figure — all clusters, one row per cluster
# =============================================================================
print("\nStep 5: Combined multi-panel figure...")

all_tfs_ordered = []
for cl in EOS_CLUSTERS:
    tfs = [tf for tf, c in selected.items() if c == cl]
    all_tfs_ordered.extend(tfs)

n_total = len(all_tfs_ordered)
if n_total > 0:
    n_cols = min(TOP_PER_CLUSTER, max(len([t for t, c in selected.items() if c == cl])
                                       for cl in EOS_CLUSTERS))
    n_rows = len(EOS_CLUSTERS)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 3.5 * n_rows))

    for row_i, cl in enumerate(EOS_CLUSTERS):
        tfs = [tf for tf, c in selected.items() if c == cl]
        for col_j in range(n_cols):
            ax = axes[row_i, col_j] if n_rows > 1 else axes[col_j]
            if col_j >= len(tfs):
                ax.axis("off")
                continue

            tf = tfs[col_j]
            plot_df = pd.DataFrame({
                "AUC": auc[tf].values,
                "condition": meta["condition"].values,
                "sample": meta["orig.ident"].values if "orig.ident" in meta.columns else "unknown",
            })
            plot_df["cond_short"] = plot_df["condition"].map(COND_SHORT)

            order = [COND_SHORT[c] for c in COND_ORDER if c in plot_df["condition"].values]
            palette = [COND_COLORS[c] for c in COND_ORDER if c in plot_df["condition"].values]

            sns.violinplot(data=plot_df, x="cond_short", y="AUC", order=order,
                           palette=palette, cut=0, inner=None, linewidth=0.6,
                           alpha=0.7, ax=ax)

            sample_means = plot_df.groupby(["cond_short", "sample"])["AUC"].mean().reset_index()
            sns.stripplot(data=sample_means, x="cond_short", y="AUC", order=order,
                          color="black", size=3, alpha=0.7, jitter=0.12, ax=ax)

            ax.set_title(f"{tf} ({cl})", fontsize=10, fontweight="bold", fontstyle="italic")
            ax.set_xlabel("")
            ax.set_ylabel("AUCell" if col_j == 0 else "")
            ax.tick_params(axis="x", labelsize=8, rotation=0)
            ax.tick_params(axis="y", labelsize=8)

    fig.suptitle("Regulon Activity by Condition (per-sample means as dots)",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    out_pdf = os.path.join(FIG_DIR, "38_regulon_violin_all_clusters.pdf")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_pdf}")

# =============================================================================
# 6. Statistical tests — Kruskal-Wallis across conditions per regulon
# =============================================================================
print("\nStep 6: Statistical tests...")

stat_rows = []
for tf in all_tfs_ordered:
    groups = []
    for cond in COND_ORDER:
        mask = meta["condition"] == cond
        if mask.sum() > 0:
            groups.append(auc.loc[mask, tf].values)

    if len(groups) >= 2:
        h_stat, p_val = stats.kruskal(*groups)
    else:
        h_stat, p_val = np.nan, np.nan

    cl = selected.get(tf, "")
    stat_rows.append({
        "regulon": tf, "cluster": cl,
        "kruskal_H": h_stat, "kruskal_p": p_val,
    })

    # Pairwise: WT vs KO within each sex
    for sex_label, wt_cond, ko_cond in [
        ("Female", "WT_Female_HFD", "KO_Female_HFD"),
        ("Male",   "WT_Male_HFD",   "KO_Male_HFD"),
    ]:
        wt_mask = meta["condition"] == wt_cond
        ko_mask = meta["condition"] == ko_cond
        if wt_mask.sum() > 0 and ko_mask.sum() > 0:
            u_stat, u_p = stats.mannwhitneyu(
                auc.loc[wt_mask, tf].values,
                auc.loc[ko_mask, tf].values,
                alternative="two-sided"
            )
            stat_rows[-1][f"WTvsKO_{sex_label}_U"] = u_stat
            stat_rows[-1][f"WTvsKO_{sex_label}_p"] = u_p

stat_df = pd.DataFrame(stat_rows)
stat_csv = os.path.join(BASE, "scenic", "38_regulon_condition_stats.csv")
stat_df.to_csv(stat_csv, index=False)
print(f"  Saved: {stat_csv}")
print(stat_df.to_string(index=False))

print("\n38 complete.")
