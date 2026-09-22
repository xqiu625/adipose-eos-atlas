#!/usr/bin/env python3
"""
16e_sex_diff_regulons.py
Eosinophil Sex Differential Regulons heatmap (Figure 5C equivalent).

Identifies regulons enriched in WT Female vs WT Male eosinophils,
then shows z-scored AUCell means across WT♀, WT♂, KO♀, KO♂.

SLURM: --partition=epyc --mem=16G --cpus-per-task=2 --time=00:30:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"])

BASE      = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
META_FILE = os.path.join(BASE, "scenic", "eos_cell_metadata.csv")
FIG_DIR   = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)
EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]

TOP_FEMALE  = 2      # top female-enriched regulons
TOP_MALE    = 3      # top male-enriched regulons
PVAL_THRESH = 0.05
ZSCORE_CAP  = 2.5

# Validated TFs: consistent direction in both SCENIC regulon activity AND
# SCENIC target gene expression (module score, naive-only, Wilcoxon p<0.05)
# Female: Tbx21 (T-bet, type-1 immunity), Stat5a (IL-5 signaling)
# Male: Gata2 (progenitor identity), Xbp1 (UPR/secretory), Stat3 (IL-6/JAK-STAT)
FIXED_FEMALE = ["Tbx21", "Stat5a"]
FIXED_MALE   = ["Gata2", "Xbp1", "Stat3"]

# No hand-picked TF exclusion list: this script does not have regulon
# target-gene-set sizes available (unlike 16b_scenic_eos_mac_heatmap.py) to
# apply an objective size-based filter, so every regulon is tested and the
# FDR-corrected q-value threshold below is the only selection criterion.

print("=" * 60)
print("16e_sex_diff_regulons.py")
print("=" * 60)

# =============================================================================
# Load data — merge per-cluster AUCell matrices
# =============================================================================
print("\nLoading per-cluster AUCell scores...")
auc_parts = []
union_regulons = set()
missing_clusters = []
for cl in EOS_CLUSTERS:
    auc_path = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(auc_path):
        df = pd.read_csv(auc_path, index_col=0)
        print(f"  {cl}: {df.shape[0]} cells x {df.shape[1]} regulons")
        union_regulons |= set(df.columns)
        auc_parts.append(df)
    else:
        print(f"  WARNING: {auc_path} not found")
        missing_clusters.append(cl)

if missing_clusters:
    raise FileNotFoundError(
        f"Missing auc_matrix.csv for cluster(s) {missing_clusters} -- the "
        f"downstream intersection is described as 'regulons scored in all "
        f"{len(EOS_CLUSTERS)} clusters,' which would be false if computed "
        f"over only {len(EOS_CLUSTERS) - len(missing_clusters)} of them. "
        f"Regenerate the missing cluster's SCENIC output rather than "
        f"silently continuing with an incomplete analysis.")

# SCENIC is inferred independently per cluster, so each cluster's regulon set
# differs. An outer join + fillna(0.0) would record "not inferred in this
# cluster's independently fitted network" as if it were a real AUCell score
# of zero -- indistinguishable from "scored and found inactive" -- which
# systematically makes regulons look cluster-specific. Use an inner join
# instead: keep only regulons that were actually scored in every cluster, so
# every remaining value is a genuine AUCell score, not an inference artifact.
auc = pd.concat(auc_parts, join="inner")
print(f"  Combined: {auc.shape[0]} cells x {auc.shape[1]} regulons "
      f"(inner join: regulons scored in ALL {len(EOS_CLUSTERS)} clusters; "
      f"{len(union_regulons)} total distinct regulons existed across clusters, "
      f"{auc.shape[1]} survive the intersection)")

print("Loading metadata...")
meta = pd.read_csv(META_FILE)
# set barcode as index
if "barcode" in meta.columns:
    meta = meta.set_index("barcode")
print(f"  {meta.shape[0]} cells")
print(f"  Sex values: {meta['sex'].unique().tolist()}")
print(f"  Genotype values: {meta['genotype'].unique().tolist()}")

# Align
common = auc.index.intersection(meta.index)
auc  = auc.loc[common]
meta = meta.loc[common]
print(f"  Aligned: {len(common)} cells")

# Strip (+ ) suffix from regulon names -> TF names
def strip_suffix(col):
    return col.split("(")[0].strip()

auc.columns = [strip_suffix(c) for c in auc.columns]
# Remove duplicates (keep max)
auc = auc.T.groupby(level=0).max().T

# =============================================================================
# Group by sex x genotype — EXCLUDE infected samples for clean sex comparison
# =============================================================================
GROUPS = [
    ("Female", "WT"),
    ("Male",   "WT"),
    ("Female", "KO"),
    ("Male",   "KO"),
]
GROUP_LABELS = ["WT ♀", "WT ♂", "KO ♀", "KO ♂"]

# Determine infection column name
if "infection" in meta.columns:
    inf_col = "infection"
elif "condition" in meta.columns:
    inf_col = "condition"
else:
    inf_col = None

if inf_col == "infection":
    naive_mask = meta[inf_col] == "HFD"
    print(f"  Filtering to naive (HFD-only) samples using '{inf_col}' column")
    print(f"  Infection values: {meta[inf_col].unique().tolist()}")
    print(f"  Naive cells: {naive_mask.sum()} / {len(meta)}")
elif inf_col == "condition":
    # Exclude conditions containing "Infected"
    naive_mask = ~meta[inf_col].str.contains("Infected", case=False, na=False)
    print(f"  Filtering to naive samples using '{inf_col}' column")
    print(f"  Condition values: {meta[inf_col].unique().tolist()}")
    print(f"  Naive cells: {naive_mask.sum()} / {len(meta)}")
else:
    naive_mask = pd.Series(True, index=meta.index)
    print("  WARNING: No infection/condition column found, using all cells")

group_cells = {}
for sex, geno in GROUPS:
    mask = (meta["sex"] == sex) & (meta["genotype"] == geno) & naive_mask
    group_cells[(sex, geno)] = meta.index[mask].tolist()
    print(f"  {sex} {geno} (naive only): {mask.sum()} cells")

# =============================================================================
# Statistical test: WT Female vs WT Male
# =============================================================================
print("\nTesting WT Female vs WT Male...")

if "sample" not in meta.columns:
    raise ValueError("Expected a 'sample' column in eos_cell_metadata.csv -- "
                     "testing must be done on per-sample summaries, not "
                     "individual cells (see below).")

wt_f_cells = group_cells[("Female", "WT")]
wt_m_cells = group_cells[("Male",   "WT")]
wt_f_samples = sorted(meta.loc[wt_f_cells, "sample"].unique())
wt_m_samples = sorted(meta.loc[wt_m_cells, "sample"].unique())
print(f"  WT Female samples ({len(wt_f_samples)}): {wt_f_samples}")
print(f"  WT Male samples   ({len(wt_m_samples)}): {wt_m_samples}")

# Per-sample mean AUCell ("pseudobulk" summary), not per-cell values: with
# only ~3 biological samples per sex, testing across the thousands of
# individual cells instead treats cells as independent replicates
# (pseudoreplication) and produces artificially small p-values. The unit of
# replication here is the sample.
auc_sample = auc.copy()
auc_sample["sample"] = meta.loc[auc.index, "sample"].values
wt_f_sample_means = auc_sample.loc[wt_f_cells].groupby("sample").mean(numeric_only=True)
wt_m_sample_means = auc_sample.loc[wt_m_cells].groupby("sample").mean(numeric_only=True)

results = []
for tf in auc.columns:
    f_vals = wt_f_sample_means[tf].values
    m_vals = wt_m_sample_means[tf].values
    _, pval = stats.mannwhitneyu(f_vals, m_vals, alternative="two-sided")
    mean_f = f_vals.mean()
    mean_m = m_vals.mean()
    fc = (mean_f + 1e-6) / (mean_m + 1e-6)
    results.append({"tf": tf, "pval": pval, "mean_f": mean_f,
                    "mean_m": mean_m, "fc": fc,
                    "n_samples_f": len(f_vals), "n_samples_m": len(m_vals)})

res = pd.DataFrame(results)
res["log2fc"] = np.log2(res["fc"])

# BH-FDR correction across all tested regulons -- raw p<0.05 alone does not
# account for the multiple-testing burden of testing every regulon.
res["qval"] = np.nan
valid = res["pval"].notna()
if valid.any():
    res.loc[valid, "qval"] = multipletests(res.loc[valid, "pval"], method="fdr_bh")[1]

sig = res[res["qval"] < PVAL_THRESH].copy()
print(f"  Significant regulons (FDR q < {PVAL_THRESH}): {len(sig)}")
n_f, n_m = len(wt_f_samples), len(wt_m_samples)
if min(n_f, n_m) <= 3:
    print(f"  NOTE: Mann-Whitney on n={n_f} vs n={n_m} samples has a small-n "
          f"exact two-tailed p-value floor (e.g. 0.1 for a balanced 3-vs-3 "
          f"comparison) that may sit above {PVAL_THRESH} regardless of effect "
          f"size -- see per-TF pval/qval below rather than expecting 'sig' "
          f"to be populated.")

res_by_tf = res.set_index("tf")


def _validate_fixed_tfs(tf_list, expect_female_higher):
    """Predeclared selection criteria for the fixed/"validated" TF list: the
    direction of effect in the CURRENT run's sample-level data must match the
    claimed direction. This is the hard gate (always computable, unlike
    nominal significance at n=3 vs n=3 -- see NOTE above); failures are
    dropped and reported explicitly rather than displayed anyway."""
    passed, failed = [], []
    for tf in tf_list:
        if tf not in res_by_tf.index:
            failed.append((tf, "not found in tested regulons"))
            continue
        row = res_by_tf.loc[tf]
        direction_ok = (row["mean_f"] > row["mean_m"]) == expect_female_higher
        if direction_ok:
            passed.append(tf)
        else:
            failed.append((tf, f"direction mismatch (mean_f={row['mean_f']:.4f}, "
                               f"mean_m={row['mean_m']:.4f}, qval={row['qval']:.4f})"))
    return passed, failed


if FIXED_FEMALE is not None and FIXED_MALE is not None:
    # Use manually specified TFs from reference figure, intersected with the
    # predeclared direction criterion computed from this run's data.
    available = set(auc.columns)
    fixed_female_avail = [tf for tf in FIXED_FEMALE if tf in available]
    fixed_male_avail   = [tf for tf in FIXED_MALE   if tf in available]
    missing = [tf for tf in FIXED_FEMALE + FIXED_MALE if tf not in available]
    if missing:
        print(f"  WARNING: TFs not found in AUCell data: {missing}")

    female_top, female_failed = _validate_fixed_tfs(fixed_female_avail, expect_female_higher=True)
    male_top,   male_failed   = _validate_fixed_tfs(fixed_male_avail,   expect_female_higher=False)

    for tf, reason in female_failed + male_failed:
        print(f"  EXCLUDED fixed TF {tf}: {reason}")

    print(f"  Using fixed TFs (direction-validated against this run) — "
          f"Female: {female_top}, Male: {male_top}")
    for tf in female_top + male_top:
        row = res_by_tf.loc[tf]
        print(f"    {tf}: mean_f={row['mean_f']:.4f}, mean_m={row['mean_m']:.4f}, "
              f"pval={row['pval']:.4f}, qval={row['qval']:.4f}")
else:
    female_top = (sig[sig["log2fc"] > 0]
                  .nlargest(TOP_FEMALE, "log2fc")["tf"].tolist())
    male_top   = (sig[sig["log2fc"] < 0]
                  .nsmallest(TOP_MALE, "log2fc")["tf"].tolist())
    print(f"  Female-enriched: {female_top}")
    print(f"  Male-enriched:   {male_top}")

selected = female_top + male_top

if not selected:
    raise ValueError("No significant sex-differential regulons found. "
                     "Try lowering PVAL_THRESH.")

# =============================================================================
# Build heatmap matrix
# =============================================================================
group_means = {}
for (sex, geno), label in zip(GROUPS, GROUP_LABELS):
    cells = group_cells[(sex, geno)]
    group_means[label] = auc.loc[cells, selected].mean()

mat = pd.DataFrame(group_means).T   # groups x tfs

# Row-wise z-score across groups (per regulon)
row_mean = mat.mean(axis=0)
row_std  = mat.std(axis=0).replace(0, 1)
zmat = ((mat - row_mean) / row_std).clip(-ZSCORE_CAP, ZSCORE_CAP)

# =============================================================================
# Plot
# =============================================================================
n_rows = len(selected)
n_cols = len(GROUP_LABELS)
n_female = len(female_top)

cell_h    = 0.90
cell_w    = 0.80
right_pad = 1.5
top_pad   = 0.8
bot_pad   = 0.9
sidebar_w = 0.18   # sidebar bar width in inches
label_w   = 1.8    # space for y-axis label text
gap       = 0.10   # gap between sidebar right edge and heatmap

# Layout: [label_w][sidebar_w][gap][heatmap][right_pad]
left_pad = label_w + sidebar_w + gap

fig_h = n_rows * cell_h + top_pad + bot_pad
fig_w = left_pad + n_cols * cell_w + right_pad

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

ax_left   = left_pad / fig_w
ax_bottom = bot_pad / fig_h
ax_width  = (n_cols * cell_w) / fig_w
ax_height = (n_rows * cell_h) / fig_h
ax.set_position([ax_left, ax_bottom, ax_width, ax_height])

# Heatmap (tfs x groups)
plot_data = zmat.T.values   # tfs x groups

im = ax.imshow(plot_data, cmap=ZSCORE_CMAP, aspect="auto",
               vmin=-ZSCORE_CAP, vmax=ZSCORE_CAP)

# Axes labels
ax.set_xticks(range(n_cols))
ax.set_xticklabels(GROUP_LABELS, fontsize=11, fontweight="bold")
ax.set_yticks(range(n_rows))
ax.set_yticklabels(selected, fontsize=12, ha="right", fontweight="bold")
ax.yaxis.set_ticks_position("left")
ax.tick_params(axis="both", length=0)

# Vertical separator between WT and KO
ax.axvline(x=1.5, color="black", linewidth=2.0)

# White grid lines between rows
for i in range(1, n_rows):
    ax.axhline(y=i - 0.5, color="white", linewidth=0.8)
for j in range(1, n_cols):
    if j != 2:   # skip separator position (already drawn)
        ax.axvline(x=j - 0.5, color="white", linewidth=0.8)

# Color sidebar: positioned just left of heatmap, after label space
FEMALE_COLOR = "#b084c8"
MALE_COLOR   = "#8fbc3b"

sidebar_colors = ([FEMALE_COLOR] * n_female + [MALE_COLOR] * len(male_top))
sb_left = (label_w + 0.05) / fig_w   # after label space + small margin
sb_w_fig = sidebar_w / fig_w
sb_ax = fig.add_axes([sb_left, ax_bottom, sb_w_fig, ax_height], zorder=0)
for i, color in enumerate(sidebar_colors):
    sb_ax.axhspan(n_rows - i - 1, n_rows - i, color=color, linewidth=0)
sb_ax.set_xlim(0, 1)
sb_ax.set_ylim(0, n_rows)
sb_ax.axis("off")
ax.set_zorder(1)
ax.patch.set_alpha(0)   # transparent background so sidebar shows through

# Title. With FIXED_FEMALE/FIXED_MALE set, this panel is a preselected,
# direction-filtered visualization, NOT a statistically supported
# differential-expression result: at n=3 vs 3 samples the Mann-Whitney test
# above cannot reach even raw p<0.05 (exact two-tailed floor of 0.1), so
# every q-value here is nonsignificant regardless of effect size. Label the
# figure accordingly rather than calling it "differential" or "validated" --
# see res_by_tf/qval for the actual (nonsignificant) statistics.
using_fixed_tfs = FIXED_FEMALE is not None and FIXED_MALE is not None
if using_fixed_tfs:
    title = "Eosinophil Sex-Associated Regulons (preselected, descriptive)"
    caption = ("Preselected TFs shown by direction of effect only "
              f"(n={n_f} vs {n_m} samples; q-values nonsignificant at this "
              "sample size -- not a statistical differential-expression claim)")
else:
    title = "Eosinophil Sex Differential Regulons (data-driven, FDR q<0.05)"
    caption = f"Data-driven selection, FDR q < {PVAL_THRESH}"

ax.set_title(title, fontsize=12, fontweight="bold", pad=20)
ax.text(0.5, 1.06, caption, transform=ax.transAxes,
        ha="center", va="bottom", fontsize=7.5, style="italic", color="dimgrey")

# Colorbar
cbar_ax = fig.add_axes([ax_left + ax_width + 0.02,
                         ax_bottom + ax_height * 0.2,
                         0.025,
                         ax_height * 0.6])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("Z-score", fontsize=9)
cbar.set_ticks([-2, -1, 0, 1, 2])
cbar.ax.tick_params(labelsize=8)

# Legend
import matplotlib.patches as mpatches
female_label = ("Higher in WT♀ (direction only)" if using_fixed_tfs
               else "Female-enriched in WT (FDR q<0.05)")
male_label = ("Higher in WT♂ (direction only)" if using_fixed_tfs
             else "Male-enriched in WT (FDR q<0.05)")
female_patch = mpatches.Patch(color=FEMALE_COLOR, label=female_label)
male_patch   = mpatches.Patch(color=MALE_COLOR,   label=male_label)
ax.legend(handles=[female_patch, male_patch], fontsize=8,
          loc="lower right", bbox_to_anchor=(1.0, -0.22),
          frameon=False)

suffix  = "_preselected_descriptive" if using_fixed_tfs else "_datadriven"
out_pdf = os.path.join(FIG_DIR, f"Fig5_sex_diff_regulons{suffix}.pdf")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved: {out_pdf}")

res.to_csv(os.path.join(BASE, "scenic", "sex_diff_regulons_stats.csv"), index=False)
print("16e complete.")
