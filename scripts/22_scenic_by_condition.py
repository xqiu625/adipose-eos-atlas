#!/usr/bin/env python3
# ==============================================================================
# 22_scenic_by_condition.py
# Compare regulon activity across conditions, not just clusters
# (addresses Reviewer 1 points 20 and 22)
# ==============================================================================
#
# SLURM: --partition=batch --mem=64G --cpus-per-task=8 --time=2:00:00
# Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
#
# Inputs:
#   - scenic/eos_scenic_output.loom  (AUCell scores from pySCENIC)
#   - scenic/eos_cell_metadata.csv   (cell metadata with condition info)
#
# Outputs:
#   - scenic/figures/22_*.pdf         (condition-specific regulon figures)
#   - scenic/regulon_condition_*.csv  (statistical results)
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
from scipy import stats
from itertools import combinations
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
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Shared palettes (matching R scripts)
EOS_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4"
}
COND_COLORS = {
    "WT_Female_HFD": "#B2182B", "WT_Male_HFD": "#2166AC",
    "KO_Female_HFD": "#E08080", "KO_Male_HFD": "#6BAED6",
    "WT_Female_HFDInfected": "#D6604D", "KO_Female_HFDInfected": "#FC8D59"
}

# -- Configuration -------------------------------------------------------------
BASE = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
LOOM_PATH = os.path.join(BASE, "scenic", "eos_scenic_output.loom")
META_PATH = os.path.join(BASE, "scenic", "eos_cell_metadata.csv")
H5AD_PATH = os.path.join(BASE, "scenic", "eos_for_scenic.h5ad")
FIG_DIR = os.path.join(BASE, "scenic", "figures")
OUT_DIR = os.path.join(BASE, "scenic")
os.makedirs(FIG_DIR, exist_ok=True)

print("=" * 60)
print("22_scenic_by_condition.py")
print("Compare regulon activity across conditions")
print("=" * 60)
print()

# ==============================================================================
# Step 1: Load data
# ==============================================================================
print("Step 1: Loading data...")

# Load cell metadata
print(f"  Loading metadata: {META_PATH}")
if os.path.exists(META_PATH):
    cell_meta = pd.read_csv(META_PATH)
    # The CSV from 14_scenic_export.R has: cluster,condition,sex,genotype,sample,barcode,UMAP_1,UMAP_2
    # Use barcode as index to align with AUCell matrix
    if "barcode" in cell_meta.columns:
        cell_meta.index = cell_meta["barcode"]
    print(f"  Metadata: {cell_meta.shape[0]} cells, columns: {list(cell_meta.columns)}")
else:
    # Fallback: load from h5ad
    print(f"  Metadata CSV not found. Loading from h5ad: {H5AD_PATH}")
    adata_tmp = sc.read_h5ad(H5AD_PATH)
    cell_meta = adata_tmp.obs.copy()
    del adata_tmp
    print(f"  Metadata: {cell_meta.shape[0]} cells, columns: {list(cell_meta.columns)}")

# Identify relevant columns (cluster, condition, sex, genotype, infection, sample)
cluster_col = "cluster" if "cluster" in cell_meta.columns else "eos_cluster"
condition_col = "condition" if "condition" in cell_meta.columns else None
sex_col = "sex" if "sex" in cell_meta.columns else None
genotype_col = "genotype" if "genotype" in cell_meta.columns else None
infection_col = "infection" if "infection" in cell_meta.columns else None
sample_col = "sample" if "sample" in cell_meta.columns else "orig.ident"

print(f"  Cluster column: {cluster_col}")
print(f"  Condition column: {condition_col}")
print(f"  Sample column: {sample_col}")

# If condition not available but genotype/sex/infection are, construct it
if condition_col is None and genotype_col and sex_col and infection_col:
    cell_meta["condition"] = (cell_meta[genotype_col] + "_" +
                              cell_meta[sex_col] + "_" +
                              cell_meta[infection_col])
    condition_col = "condition"
    print(f"  Constructed condition column from genotype+sex+infection")

# Load SCENIC AUCell scores (try CSV first, fall back to loom)
AUCELL_CSV = os.path.join(BASE, "scenic", "eos_scenic_aucell.csv")

if os.path.exists(AUCELL_CSV):
    print(f"  Loading AUCell CSV: {AUCELL_CSV}")
    auc_mtx = pd.read_csv(AUCELL_CSV, index_col=0)
    auc_mtx.columns = [c.strip() for c in auc_mtx.columns]
    regulon_names = list(auc_mtx.columns)
elif os.path.exists(LOOM_PATH):
    print(f"  Loading loom: {LOOM_PATH}")
    import loompy
    lf = loompy.connect(LOOM_PATH, mode="r")
    if hasattr(lf.ca["RegulonsAUC"], "dtype") and lf.ca["RegulonsAUC"].dtype.names:
        regulon_names = list(lf.ca["RegulonsAUC"].dtype.names)
        cell_ids = lf.ca["CellID"]
        auc_data = np.column_stack([lf.ca["RegulonsAUC"][r] for r in regulon_names])
        auc_mtx = pd.DataFrame(auc_data, index=cell_ids, columns=regulon_names)
    else:
        cell_ids = lf.ca["CellID"]
        auc_mtx = pd.DataFrame(lf.ca["RegulonsAUC"], index=cell_ids)
        regulon_names = list(auc_mtx.columns)
    lf.close()
else:
    raise FileNotFoundError(
        f"Neither {AUCELL_CSV} nor {LOOM_PATH} found. Run pySCENIC first."
    )

# Load regulon-target gene mappings from regulons.csv
REGULONS_CSV = os.path.join(BASE, "scenic", "regulons.csv")
regulon_targets = {}
try:
    if os.path.exists(REGULONS_CSV):
        print(f"  Loading regulon targets from: {REGULONS_CSV}")
        reg_df = pd.read_csv(REGULONS_CSV)
        if "TF" in reg_df.columns:
            for tf_name in reg_df["TF"].unique():
                tf_rows = reg_df[reg_df["TF"] == tf_name]
                targets = set()
                for _, row in tf_rows.iterrows():
                    if "TargetGenes" in row and pd.notna(row["TargetGenes"]):
                        tg = str(row["TargetGenes"])
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
                    regulon_targets[f"{tf_name}(+)"] = list(targets)[:50]
        print(f"  Found target genes for {len(regulon_targets)} regulons")
except Exception as e:
    print(f"  Warning: Could not extract regulon targets: {e}")

print(f"  AUCell matrix: {auc_mtx.shape[0]} cells x {auc_mtx.shape[1]} regulons")

# Align cells
common_cells = cell_meta.index.intersection(auc_mtx.index)
print(f"  Common cells: {len(common_cells)}")

if len(common_cells) == 0:
    raise ValueError("No overlapping cell IDs between metadata and AUCell matrix. "
                     "Check that cell IDs match.")

cell_meta = cell_meta.loc[common_cells, :]
auc_mtx = auc_mtx.loc[common_cells, :]

clusters = sorted(cell_meta[cluster_col].unique())
conditions = sorted(cell_meta[condition_col].unique()) if condition_col else []
print(f"  Clusters: {clusters}")
print(f"  Conditions: {conditions}")
print()

# ==============================================================================
# Helper functions
# ==============================================================================


def cohens_d(group1, group2):
    """Compute Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    if n1 < 2 or n2 < 2:
        return 0.0
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(group1) - np.mean(group2)) / pooled_std


def multipletests_bh(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.array(pvals)
    n = len(pvals)
    if n == 0:
        return np.array([])
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    adjusted = np.zeros(n)
    cum_min = 1.0
    for i in range(n - 1, -1, -1):
        adjusted_p = sorted_pvals[i] * n / (i + 1)
        cum_min = min(cum_min, adjusted_p)
        adjusted[sorted_idx[i]] = min(cum_min, 1.0)
    return adjusted


# ==============================================================================
# Analysis 1: Regulon activity by condition (within each cluster)
# ==============================================================================
print("Analysis 1: Regulon activity by condition within each cluster...")

all_kw_results = []

for cl in clusters:
    cl_mask = cell_meta[cluster_col] == cl
    cl_meta = cell_meta.loc[cl_mask]
    cl_auc = auc_mtx.loc[cl_mask]

    cl_conditions = sorted(cl_meta[condition_col].unique())
    if len(cl_conditions) < 2:
        print(f"  {cl}: Only {len(cl_conditions)} condition(s), skipping")
        continue

    print(f"  {cl}: {cl_mask.sum()} cells, {len(cl_conditions)} conditions")

    for reg in regulon_names:
        groups = [cl_auc.loc[cl_meta[condition_col] == cond, reg].values
                  for cond in cl_conditions]
        # Filter out groups with < 3 cells
        valid_groups = [g for g in groups if len(g) >= 3]
        if len(valid_groups) < 2:
            continue

        try:
            stat, pval = stats.kruskal(*valid_groups)
        except ValueError:
            continue

        all_kw_results.append({
            "cluster": cl,
            "regulon": reg,
            "kruskal_stat": stat,
            "p_value": pval,
            "n_conditions": len(valid_groups),
            "n_cells": sum(len(g) for g in valid_groups)
        })

kw_df = pd.DataFrame(all_kw_results)

if len(kw_df) > 0:
    kw_df["fdr"] = multipletests_bh(kw_df["p_value"].values)

    # Identify condition-responsive regulons per cluster
    sig_kw = kw_df[kw_df["fdr"] < 0.05]
    print(f"  Total tests: {len(kw_df)}")
    print(f"  Significant (FDR < 0.05): {len(sig_kw)}")

    for cl in clusters:
        cl_sig = sig_kw[sig_kw["cluster"] == cl]
        if len(cl_sig) > 0:
            top_regs = cl_sig.nsmallest(5, "fdr")["regulon"].tolist()
            print(f"    {cl}: {len(cl_sig)} sig regulons, top: {', '.join(top_regs)}")

    # Save
    kw_path = os.path.join(OUT_DIR, "regulon_condition_kruskal.csv")
    kw_df.to_csv(kw_path, index=False)
    print(f"  Saved: {kw_path}")

    # Pairwise Mann-Whitney for significant regulons
    print("  Running pairwise Mann-Whitney for significant regulons...")
    pairwise_results = []

    for _, row in sig_kw.iterrows():
        cl = row["cluster"]
        reg = row["regulon"]
        cl_mask = cell_meta[cluster_col] == cl
        cl_meta_sub = cell_meta.loc[cl_mask]
        cl_auc_sub = auc_mtx.loc[cl_mask]

        cl_conditions = sorted(cl_meta_sub[condition_col].unique())
        for cond1, cond2 in combinations(cl_conditions, 2):
            g1 = cl_auc_sub.loc[cl_meta_sub[condition_col] == cond1, reg].values
            g2 = cl_auc_sub.loc[cl_meta_sub[condition_col] == cond2, reg].values
            if len(g1) < 3 or len(g2) < 3:
                continue
            try:
                stat, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")
                d = cohens_d(g1, g2)
            except ValueError:
                continue
            pairwise_results.append({
                "cluster": cl,
                "regulon": reg,
                "condition1": cond1,
                "condition2": cond2,
                "U_stat": stat,
                "p_value": pval,
                "cohens_d": d,
                "mean1": np.mean(g1),
                "mean2": np.mean(g2),
                "n1": len(g1),
                "n2": len(g2)
            })

    pw_df = pd.DataFrame(pairwise_results)
    if len(pw_df) > 0:
        pw_df["fdr"] = multipletests_bh(pw_df["p_value"].values)
        pw_path = os.path.join(OUT_DIR, "regulon_condition_pairwise.csv")
        pw_df.to_csv(pw_path, index=False)
        print(f"  Saved: {pw_path} ({len(pw_df)} tests)")
else:
    print("  No Kruskal-Wallis results to report")
    kw_df = pd.DataFrame()
    pw_df = pd.DataFrame()

print()

# ==============================================================================
# Analysis 2: Condition-specific regulon heatmap
# ==============================================================================
print("Analysis 2: Condition x cluster regulon heatmap...")

# Compute mean AUCell per condition x cluster
mean_auc = {}
for cl in clusters:
    for cond in conditions:
        mask = (cell_meta[cluster_col] == cl) & (cell_meta[condition_col] == cond)
        n_cells = mask.sum()
        if n_cells >= 3:
            mean_auc[(cl, cond)] = auc_mtx.loc[mask].mean()

mean_auc_df = pd.DataFrame(mean_auc).T
mean_auc_df.index = pd.MultiIndex.from_tuples(mean_auc_df.index,
                                               names=["cluster", "condition"])

# Z-score across condition x cluster combinations within each regulon
z_scored = (mean_auc_df - mean_auc_df.mean()) / mean_auc_df.std()

# Select top regulons: most variable across conditions
regulon_var = z_scored.var(axis=0).sort_values(ascending=False)
top_n_regs = min(30, len(regulon_var))
top_regulons = regulon_var.head(top_n_regs).index.tolist()

# Create column labels
z_plot = z_scored[top_regulons].copy()
z_plot.index = [f"{cl}_{cond}" for cl, cond in z_plot.index]

# Figure 1: Heatmap
print(f"  Plotting heatmap with top {top_n_regs} regulons...")

fig, ax = plt.subplots(figsize=(max(14, len(z_plot) * 0.6),
                                max(8, top_n_regs * 0.35)))
sns.heatmap(
    z_plot.T,
    cmap="RdBu_r",
    center=0,
    linewidths=0.3,
    linecolor="white",
    ax=ax,
    cbar_kws={"label": "Z-scored mean AUCell"},
    xticklabels=True,
    yticklabels=True
)
ax.set_title("Regulon Activity: Condition x Cluster", fontsize=14, pad=15)
ax.set_xlabel("Cluster_Condition", fontsize=11)
ax.set_ylabel("Regulon", fontsize=11)
ax.tick_params(axis="x", labelsize=7, rotation=90)
ax.tick_params(axis="y", labelsize=8)
plt.tight_layout()

heatmap_path = os.path.join(FIG_DIR, "22_regulon_condition_heatmap.pdf")
fig.savefig(heatmap_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {heatmap_path}")

# Save mean AUCell table
mean_auc_path = os.path.join(OUT_DIR, "regulon_mean_auc_by_condition_cluster.csv")
mean_auc_df.to_csv(mean_auc_path)
print(f"  Saved: {mean_auc_path}")

print()

# ==============================================================================
# Analysis 3: Violin plots for top condition-responsive regulons
# ==============================================================================
print("Analysis 3: Violin plots for top condition-responsive regulons...")

# Identify top 10 most variable regulons across conditions (globally)
top10_regulons = regulon_var.head(10).index.tolist()
print(f"  Top 10 variable regulons: {', '.join(top10_regulons)}")

# Build long-format data
violin_data = []
for reg in top10_regulons:
    for cl in clusters:
        for cond in conditions:
            mask = (cell_meta[cluster_col] == cl) & (cell_meta[condition_col] == cond)
            vals = auc_mtx.loc[mask, reg].values
            for v in vals:
                violin_data.append({
                    "regulon": reg,
                    "cluster": cl,
                    "condition": cond,
                    "auc_score": v
                })

violin_df = pd.DataFrame(violin_data)

if len(violin_df) > 0:
    # Plot: one row per regulon, columns = clusters, violins = conditions
    n_regs = len(top10_regulons)
    n_cols = len(clusters)

    fig, axes = plt.subplots(n_regs, n_cols,
                              figsize=(3 * n_cols, 2.5 * n_regs),
                              squeeze=False)

    # Determine unique conditions for consistent coloring
    cond_palette = {c: COND_COLORS.get(c, "#999999") for c in conditions}

    for i, reg in enumerate(top10_regulons):
        for j, cl in enumerate(clusters):
            ax = axes[i, j]
            sub = violin_df[(violin_df["regulon"] == reg) &
                            (violin_df["cluster"] == cl)]

            if len(sub) > 0:
                conds_present = sorted(sub["condition"].unique())
                positions = list(range(len(conds_present)))
                parts_data = [sub.loc[sub["condition"] == c, "auc_score"].values
                              for c in conds_present]

                # Filter out empty groups
                valid = [(p, d, c) for p, d, c in
                         zip(positions, parts_data, conds_present) if len(d) > 0]
                if valid:
                    vp = ax.violinplot([d for _, d, _ in valid],
                                       positions=[p for p, _, _ in valid],
                                       showmeans=True, showextrema=False)
                    for idx, body in enumerate(vp["bodies"]):
                        c = valid[idx][2]
                        body.set_facecolor(cond_palette.get(c, "grey"))
                        body.set_alpha(0.7)

                    ax.set_xticks([p for p, _, _ in valid])
                    ax.set_xticklabels([c for _, _, c in valid],
                                       rotation=90, fontsize=5)

            if i == 0:
                ax.set_title(cl, fontsize=10, fontweight="bold")
            if j == 0:
                ax.set_ylabel(reg, fontsize=8)
            else:
                ax.set_ylabel("")
            ax.tick_params(labelsize=6)

    plt.suptitle("Top Condition-Responsive Regulons by Cluster",
                 fontsize=14, y=1.01)
    plt.tight_layout()

    violin_path = os.path.join(FIG_DIR, "22_regulon_condition_violins.pdf")
    fig.savefig(violin_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {violin_path}")
else:
    print("  No data for violin plots")

print()

# ==============================================================================
# Analysis 4: RELMa-dependent regulon changes (WT vs KO)
# ==============================================================================
print("Analysis 4: RELMa-dependent regulon changes (WT vs KO)...")

# Parse genotype from condition or use genotype column
if genotype_col and genotype_col in cell_meta.columns:
    geno_series = cell_meta[genotype_col]
elif condition_col:
    geno_series = cell_meta[condition_col].str.extract(r"^(WT|KO)", expand=False)
else:
    geno_series = None

if sex_col and sex_col in cell_meta.columns:
    sex_series = cell_meta[sex_col]
elif condition_col:
    sex_series = cell_meta[condition_col].str.extract(r"_(Male|Female)_", expand=False)
else:
    sex_series = None

if infection_col and infection_col in cell_meta.columns:
    inf_series = cell_meta[infection_col]
elif condition_col:
    inf_series = cell_meta[condition_col].apply(
        lambda x: "HFDInfected" if "Infected" in str(x) else "HFD"
    )
else:
    inf_series = None

relma_results = []

# Compare WT vs KO within each sex (HFD only, non-infected)
comparisons = []
if geno_series is not None and sex_series is not None and inf_series is not None:
    for sex_val in ["Female", "Male"]:
        wt_mask = ((geno_series == "WT") & (sex_series == sex_val) &
                   (inf_series == "HFD"))
        ko_mask = ((geno_series == "KO") & (sex_series == sex_val) &
                   (inf_series == "HFD"))

        if wt_mask.sum() >= 3 and ko_mask.sum() >= 3:
            comparisons.append({
                "name": f"WT_vs_KO_{sex_val}_HFD",
                "wt_mask": wt_mask,
                "ko_mask": ko_mask,
                "sex": sex_val
            })
            print(f"  Comparison: WT vs KO {sex_val} HFD "
                  f"({wt_mask.sum()} WT, {ko_mask.sum()} KO cells)")

for comp in comparisons:
    comp_name = comp["name"]

    for reg in regulon_names:
        wt_vals = auc_mtx.loc[comp["wt_mask"], reg].values
        ko_vals = auc_mtx.loc[comp["ko_mask"], reg].values

        if len(wt_vals) < 3 or len(ko_vals) < 3:
            continue

        try:
            stat, pval = stats.mannwhitneyu(wt_vals, ko_vals,
                                            alternative="two-sided")
        except ValueError:
            continue

        d = cohens_d(wt_vals, ko_vals)

        # Log2 fold change of mean AUCell
        mean_wt = np.mean(wt_vals)
        mean_ko = np.mean(ko_vals)
        if mean_ko > 0 and mean_wt > 0:
            log2fc = np.log2(mean_wt / mean_ko)
        elif mean_wt > 0:
            log2fc = np.inf
        elif mean_ko > 0:
            log2fc = -np.inf
        else:
            log2fc = 0.0

        relma_results.append({
            "comparison": comp_name,
            "sex": comp["sex"],
            "regulon": reg,
            "mean_WT": mean_wt,
            "mean_KO": mean_ko,
            "log2fc_WT_vs_KO": log2fc,
            "cohens_d": d,
            "U_stat": stat,
            "p_value": pval,
            "n_WT": len(wt_vals),
            "n_KO": len(ko_vals)
        })

relma_df = pd.DataFrame(relma_results)

if len(relma_df) > 0:
    relma_df["fdr"] = multipletests_bh(relma_df["p_value"].values)

    sig_relma = relma_df[relma_df["fdr"] < 0.05]
    print(f"  Total tests: {len(relma_df)}")
    print(f"  Significant (FDR < 0.05): {len(sig_relma)}")

    for comp_name in relma_df["comparison"].unique():
        comp_sig = sig_relma[sig_relma["comparison"] == comp_name]
        if len(comp_sig) > 0:
            top_regs = comp_sig.nsmallest(5, "fdr")["regulon"].tolist()
            print(f"    {comp_name}: {len(comp_sig)} sig, "
                  f"top: {', '.join(top_regs)}")

    # Save
    relma_path = os.path.join(OUT_DIR, "regulon_relma_dependent.csv")
    relma_df.to_csv(relma_path, index=False)
    print(f"  Saved: {relma_path}")

    # Figure 3: Volcano-style plot (effect size vs -log10 p)
    print("  Plotting volcano plots...")

    unique_comps = relma_df["comparison"].unique()
    n_comps = len(unique_comps)

    fig, axes = plt.subplots(1, n_comps, figsize=(7 * n_comps, 6), squeeze=False)

    for idx, comp_name in enumerate(unique_comps):
        ax = axes[0, idx]
        sub = relma_df[relma_df["comparison"] == comp_name].copy()

        # Replace inf log2fc with max finite value + margin
        finite_mask = np.isfinite(sub["log2fc_WT_vs_KO"])
        if finite_mask.any():
            max_abs = sub.loc[finite_mask, "log2fc_WT_vs_KO"].abs().max()
            sub.loc[~finite_mask & (sub["log2fc_WT_vs_KO"] > 0),
                    "log2fc_WT_vs_KO"] = max_abs * 1.2
            sub.loc[~finite_mask & (sub["log2fc_WT_vs_KO"] < 0),
                    "log2fc_WT_vs_KO"] = -max_abs * 1.2

        sub["neg_log10_fdr"] = -np.log10(sub["fdr"].clip(lower=1e-300))

        # Color by significance
        sub["sig_group"] = "NS"
        sub.loc[(sub["fdr"] < 0.05) & (sub["log2fc_WT_vs_KO"] > 0),
                "sig_group"] = "Up in WT"
        sub.loc[(sub["fdr"] < 0.05) & (sub["log2fc_WT_vs_KO"] < 0),
                "sig_group"] = "Up in KO"

        color_map = {"NS": "#B3B3B3", "Up in WT": "#E84A5F", "Up in KO": "#2166AC"}

        for group, color in color_map.items():
            mask = sub["sig_group"] == group
            ax.scatter(sub.loc[mask, "log2fc_WT_vs_KO"],
                       sub.loc[mask, "neg_log10_fdr"],
                       c=color, s=20, alpha=0.7, label=group, edgecolors="none")

        # Label top significant regulons
        top_sig = sub[sub["fdr"] < 0.05].nsmallest(10, "fdr")
        for _, row in top_sig.iterrows():
            ax.annotate(row["regulon"],
                        (row["log2fc_WT_vs_KO"], row["neg_log10_fdr"]),
                        fontsize=7, ha="center", va="bottom",
                        xytext=(0, 5), textcoords="offset points")

        ax.axhline(y=-np.log10(0.05), color="grey", linewidth=0.5,
                   linestyle="--")
        ax.axvline(x=0, color="grey", linewidth=0.5, linestyle="--")
        ax.set_xlabel("Log2 FC (WT / KO)", fontsize=11)
        ax.set_ylabel("-log10(FDR)", fontsize=11)
        ax.set_title(comp_name, fontsize=12, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right")

    plt.suptitle("RELMa-Dependent Regulon Changes (WT vs KO)",
                 fontsize=14, y=1.02)
    plt.tight_layout()

    volcano_path = os.path.join(FIG_DIR, "22_regulon_relma_volcano.pdf")
    fig.savefig(volcano_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {volcano_path}")

    # Figure 4: Network of RELMa-dependent regulons and target genes
    print("  Plotting RELMa-dependent regulon networks...")

    try:
        import networkx as nx

        for comp_name in unique_comps:
            comp_sig = sig_relma[sig_relma["comparison"] == comp_name]
            if len(comp_sig) == 0:
                continue

            # Take top 10 by significance
            top_regs = comp_sig.nsmallest(10, "fdr")

            fig, ax = plt.subplots(figsize=(12, 10))
            G = nx.Graph()

            for _, row in top_regs.iterrows():
                reg = row["regulon"]
                G.add_node(reg, node_type="regulon",
                           direction="up_WT" if row["log2fc_WT_vs_KO"] > 0
                           else "up_KO")

                targets = regulon_targets.get(reg, [])
                for target in targets[:8]:  # Top 8 targets per regulon
                    G.add_node(target, node_type="target")
                    G.add_edge(reg, target)

            if len(G.nodes()) == 0:
                plt.close(fig)
                continue

            pos = nx.spring_layout(G, k=2.5, iterations=50, seed=42)

            # Regulon nodes colored by direction
            up_wt_nodes = [n for n, d in G.nodes(data=True)
                           if d.get("direction") == "up_WT"]
            up_ko_nodes = [n for n, d in G.nodes(data=True)
                           if d.get("direction") == "up_KO"]
            target_nodes = [n for n, d in G.nodes(data=True)
                            if d.get("node_type") == "target"]

            nx.draw_networkx_nodes(G, pos, nodelist=up_wt_nodes, node_size=600,
                                   node_color="#E84A5F", alpha=0.9, ax=ax,
                                   label="Up in WT")
            nx.draw_networkx_nodes(G, pos, nodelist=up_ko_nodes, node_size=600,
                                   node_color="#6CA6CD", alpha=0.9, ax=ax,
                                   label="Up in KO")
            nx.draw_networkx_nodes(G, pos, nodelist=target_nodes, node_size=150,
                                   node_color="#CBAACB", alpha=0.7, ax=ax,
                                   label="Target genes")
            nx.draw_networkx_edges(G, pos, alpha=0.2, ax=ax)
            nx.draw_networkx_labels(G, pos, font_size=6, ax=ax)

            ax.set_title(f"RELMa-Dependent Regulons: {comp_name}",
                         fontsize=13, fontweight="bold")
            ax.legend(fontsize=9, loc="upper left")
            ax.axis("off")
            plt.tight_layout()

            net_path = os.path.join(
                FIG_DIR,
                f"22_regulon_relma_network_{comp_name}.pdf"
            )
            fig.savefig(net_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {net_path}")

    except ImportError:
        print("  Warning: networkx not available; skipping network plots.")

else:
    print("  No RELMa comparison data available.")
    relma_df = pd.DataFrame()

print()

# ==============================================================================
# Summary
# ==============================================================================
print("=" * 60)
print("Summary")
print("=" * 60)
print(f"  Cells:            {len(common_cells)}")
print(f"  Regulons:         {len(regulon_names)}")
print(f"  Clusters:         {len(clusters)}")
print(f"  Conditions:       {len(conditions)}")

if len(kw_df) > 0:
    print(f"  KW tests:         {len(kw_df)} (sig: {(kw_df['fdr'] < 0.05).sum()})")
if len(relma_df) > 0:
    print(f"  RELMa tests:      {len(relma_df)} "
          f"(sig: {(relma_df['fdr'] < 0.05).sum()})")

print()
print("Output files:")
print(f"  {os.path.join(OUT_DIR, 'regulon_condition_kruskal.csv')}")
print(f"  {os.path.join(OUT_DIR, 'regulon_condition_pairwise.csv')}")
print(f"  {os.path.join(OUT_DIR, 'regulon_mean_auc_by_condition_cluster.csv')}")
print(f"  {os.path.join(OUT_DIR, 'regulon_relma_dependent.csv')}")
print(f"  {os.path.join(FIG_DIR, '22_regulon_condition_heatmap.pdf')}")
print(f"  {os.path.join(FIG_DIR, '22_regulon_condition_violins.pdf')}")
print(f"  {os.path.join(FIG_DIR, '22_regulon_relma_volcano.pdf')}")
print(f"  {os.path.join(FIG_DIR, '22_regulon_relma_network_*.pdf')}")
print()
print("Done! 22_scenic_by_condition.py completed successfully.")
