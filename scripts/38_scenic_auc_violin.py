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

import math
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")


def mannwhitney_floor_p(n1, n2):
    """Minimum possible exact two-sided Mann-Whitney p-value achievable at
    sample sizes (n1, n2), regardless of effect size: 2 / C(n1+n2, n1). E.g.
    0.1 at 3-vs-3, 0.333 at 2-vs-2. Matches 22_scenic_by_condition.py."""
    return 2.0 / math.comb(n1 + n2, min(n1, n2))


MIN_GROUP_N_FOR_INFERENCE = 3  # Kruskal-Wallis heuristic, matches 22_scenic_by_condition.py

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

# No hand-picked TF exclusion list here (this script has no regulon
# target-gene-set sizes available to apply an objective size-based filter,
# same as 16e_sex_diff_regulons.py) -- selection below is MIN_AUC + uniqueness
# only.

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
    else:
        raise FileNotFoundError(
            f"Missing auc_matrix.csv for cluster {cl}. Regenerate the "
            f"missing cluster's SCENIC output rather than silently "
            f"continuing with an incomplete analysis.")

# Inner join, not outer + fillna(0.0): SCENIC is fit independently per
# cluster, so a regulon absent from one cluster's network is "not inferred,"
# not a real AUCell score of zero -- an outer join would make it look like
# every OTHER cluster's cells were scored and found inactive for that
# regulon, systematically inflating apparent cluster-specificity.
auc = pd.concat(auc_parts, join="inner")
if auc.shape[1] == 0:
    raise ValueError(
        f"Zero regulons survive the intersection across all {len(EOS_CLUSTERS)} "
        f"clusters -- each cluster's SCENIC network shares no regulon in "
        f"common. This analysis needs a common regulon set to test against; "
        f"consider rescoring AUCell for a shared/union gene-set instead of "
        f"relying on per-cluster independent inference.")
print(f"  Combined: {auc.shape[0]} cells x {auc.shape[1]} regulons "
      f"(scored in all {len(EOS_CLUSTERS)} clusters)")

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

# Resolve the per-sample identity column: 14_scenic_export.R renames
# orig.ident -> "sample" in this CSV; fall back to "orig.ident" for older
# exports. Previously this script referenced "orig.ident" directly with a
# silent fallback to the literal string "unknown" wherever it was missing,
# which would collapse the per-sample dot overlay (and any per-sample
# statistic) into a single group per condition with no warning.
if "sample" in meta.columns:
    SAMPLE_COL = "sample"
elif "orig.ident" in meta.columns:
    SAMPLE_COL = "orig.ident"
else:
    raise ValueError("eos_cell_metadata.csv has neither 'sample' nor "
                     "'orig.ident' -- cannot identify biological replicates "
                     "for per-sample statistics/plotting.")
print(f"  Using '{SAMPLE_COL}' as the sample/replicate identity column")

# Derive condition (strip trailing _N replicate)
if "condition" not in meta.columns:
    meta["condition"] = meta[SAMPLE_COL].str.replace(r"_\d+$", "", regex=True)

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
            "sample": meta[SAMPLE_COL].values,
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
                "sample": meta[SAMPLE_COL].values,
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

# Per-sample mean AUCell ("pseudobulk" summary), not per-cell values: testing
# across individual cells treats them as independent biological replicates
# (pseudoreplication) when there are only a handful of samples per condition,
# and produces artificially small p-values. The unit of replication is the
# sample.
#
# Critically, each regulon's per-sample mean is computed ONLY over cells in
# that regulon's home cluster (selected[tf]), not pooled across all
# eosinophils. `selected` holds each TF's top-scoring cluster; testing a
# cluster-selected regulon over all clusters combined would confound the
# condition effect with each sample's cluster composition -- a sample with
# proportionally fewer cells from that cluster would show a lower "mean
# AUCell" largely because of composition, not the condition being tested.
def cluster_sample_means(tf, cl, cond):
    """Per-sample mean AUCell for `tf`, restricted to cells in cluster `cl`
    and condition `cond`."""
    mask = (meta["cluster"] == cl) & (meta["condition"] == cond)
    if mask.sum() == 0:
        return np.array([])
    sub_meta = meta.loc[mask]
    sub_auc = auc.loc[mask, tf]
    return sub_auc.groupby(sub_meta[SAMPLE_COL].values).mean().values


stat_rows = []
for tf in all_tfs_ordered:
    cl = selected.get(tf, "")
    if not cl:
        continue

    # Require >=2 samples per group -- a cluster present in only one sample
    # of a condition would otherwise enter Kruskal-Wallis/Mann-Whitney as an
    # n=1 group, which is not a valid two-sample comparison.
    groups = []
    group_ns = []
    for cond in COND_ORDER:
        vals = cluster_sample_means(tf, cl, cond)
        if len(vals) >= 2:
            groups.append(vals)
            group_ns.append(len(vals))

    if len(groups) >= 2:
        h_stat, p_val = stats.kruskal(*groups)
    else:
        h_stat, p_val = np.nan, np.nan

    stat_rows.append({
        "regulon": tf, "cluster": cl,
        "kruskal_H": h_stat, "kruskal_p": p_val,
        "n_samples_per_condition": ",".join(str(n) for n in group_ns),
        # Heuristic, not an exact combinatorial floor (unlike Mann-Whitney
        # below): any group below this size cannot support a meaningful
        # significance claim, so treat as a descriptive effect-size
        # comparison rather than a confirmed significant result.
        "underpowered": (len(groups) < 2
                         or min(group_ns) < MIN_GROUP_N_FOR_INFERENCE),
    })

    # Pairwise: WT vs KO within each sex, on per-sample means within cl.
    # Require >=2 samples per side for the same reason as above.
    for sex_label, wt_cond, ko_cond in [
        ("Female", "WT_Female_HFD", "KO_Female_HFD"),
        ("Male",   "WT_Male_HFD",   "KO_Male_HFD"),
    ]:
        wt_vals = cluster_sample_means(tf, cl, wt_cond)
        ko_vals = cluster_sample_means(tf, cl, ko_cond)
        if len(wt_vals) >= 2 and len(ko_vals) >= 2:
            u_stat, u_p = stats.mannwhitneyu(
                wt_vals, ko_vals, alternative="two-sided"
            )
            floor_p = mannwhitney_floor_p(len(wt_vals), len(ko_vals))
            stat_rows[-1][f"WTvsKO_{sex_label}_U"] = u_stat
            stat_rows[-1][f"WTvsKO_{sex_label}_p"] = u_p
            stat_rows[-1][f"WTvsKO_{sex_label}_floor_p"] = floor_p
            stat_rows[-1][f"WTvsKO_{sex_label}_underpowered"] = floor_p > 0.05
            stat_rows[-1][f"WTvsKO_{sex_label}_n"] = f"{len(wt_vals)}v{len(ko_vals)}"

stat_df = pd.DataFrame(stat_rows)

# BH-FDR correction, each test family (omnibus Kruskal-Wallis, and each
# pairwise WT-vs-KO comparison) corrected separately across all regulons
# tested in that family -- raw p-values alone don't account for the
# multiple-testing burden of testing every selected regulon.
def _add_fdr(df, pcol, fdrcol):
    valid = df[pcol].notna()
    df[fdrcol] = np.nan
    if valid.any():
        df.loc[valid, fdrcol] = multipletests(df.loc[valid, pcol], method="fdr_bh")[1]

_add_fdr(stat_df, "kruskal_p", "kruskal_fdr")
stat_df["kruskal_statistically_testable"] = ~stat_df["underpowered"]
for sex_label in ["Female", "Male"]:
    pcol = f"WTvsKO_{sex_label}_p"
    if pcol in stat_df.columns:
        _add_fdr(stat_df, pcol, f"WTvsKO_{sex_label}_fdr")
        upcol = f"WTvsKO_{sex_label}_underpowered"
        # NaN (not True/False) where the pairwise test never ran at all
        # (fewer than 2 samples on one side) -- there is no test to call
        # testable or not.
        stat_df[f"WTvsKO_{sex_label}_statistically_testable"] = stat_df[upcol].map(
            lambda x: (not x) if pd.notna(x) else np.nan)

n_kw_underpowered = int(stat_df["underpowered"].sum())
print(f"  Kruskal-Wallis: {len(stat_df) - n_kw_underpowered}/{len(stat_df)} "
      f"statistically testable (smallest group >= {MIN_GROUP_N_FOR_INFERENCE} samples)")

stat_csv = os.path.join(BASE, "scenic", "38_regulon_condition_stats.csv")
stat_df.to_csv(stat_csv, index=False)
print(f"  Saved: {stat_csv}")
print(stat_df.to_string(index=False))

print("\n38 complete.")
