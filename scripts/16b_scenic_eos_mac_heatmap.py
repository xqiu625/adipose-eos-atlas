#!/usr/bin/env python3
"""
16b_scenic_eos_mac_heatmap.py
Combined Eos + Mac regulon heatmap (Fig 5 manuscript style)

Shows Eos-subcluster-unique regulons across:
  columns: C0, C1, C2, C3, C4, C5 | M1, M2 Funct., M2 Surf.

Logic adapted from 1208_Combined_Eos_Mac_Heatmap.py

SLURM: --partition=batch --mem=64G --cpus-per-task=8 --time=2:00:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import ast
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

warnings.filterwarnings("ignore")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
})

# Blue-white-red colormap matching manuscript style
ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"]
)

BASE      = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC    = os.path.join(BASE, "scenic")
FIG_DIR   = os.path.join(SCENIC, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

EOS_CLUSTERS  = ["C0", "C1", "C2", "C3", "C4", "C5"]
MAC_SUBTYPES  = ["M1", "M2_Functional", "M2_Surface"]
MAC_RENAME    = {"Mac_M1": "M1", "Mac_M2_Func": "M2_Functional", "Mac_M2_Surf": "M2_Surface"}
MAC_DISPLAY   = {"M1": "M1", "M2_Functional": "M2\nFunct.", "M2_Surface": "M2\nSurf."}

CELLS_PER_GROUP = 200
RANDOM_SEED     = 42
ZSCORE_CAP      = 2.5

# Eos-unique selection thresholds (matching relax level)
MIN_EOS_AUC        = 0.05   # lower threshold since these are global AUCell values
MIN_FOLD_CHANGE    = 1.2
MIN_UNIQUENESS     = 0.003
TOP_PER_CLUSTER    = 5      # top N per Eos cluster → shown in figure

# Regulons are excluded from selection only via an objective, reproducible
# rule computed from the data: a minimum target-gene-set size. Regulons with
# very few inferred targets are noisier and less interpretable as a "module."
# (Previously a hand-picked list of ~20 TF names -- including biologically
# relevant factors such as Stat6, Atf4, Atf3, Cebpb -- was excluded with no
# documented rationale; that list has been removed.)
MIN_REGULON_SIZE = 10

print("=" * 60)
print("16b_scenic_eos_mac_heatmap.py")
print("=" * 60)

# =============================================================================
# Step 1: Load Eos AUCell
# =============================================================================
print("\nStep 1: Loading Eos AUCell scores...")
auc_path = os.path.join(SCENIC, "eos_scenic_aucell.csv")
auc_mtx  = pd.read_csv(auc_path, index_col=0)
auc_mtx.columns = [c.strip() for c in auc_mtx.columns]
print(f"  Eos AUCell: {auc_mtx.shape[0]} cells x {auc_mtx.shape[1]} regulons")

# Load Eos metadata (cluster assignments)
meta_path = os.path.join(SCENIC, "eos_cell_metadata.csv")
clusters = None

if os.path.exists(meta_path):
    # Read without index_col so we can inspect all columns
    eos_meta = pd.read_csv(meta_path)
    print(f"  Metadata columns: {list(eos_meta.columns)}")

    # Find barcode and cluster columns by name
    bc_col  = next((c for c in eos_meta.columns if "barcode" in c.lower()), None)
    cl_col  = next((c for c in eos_meta.columns if c.lower() == "cluster"), None)

    # If not found as named columns, the CSV may have been saved with index as cluster
    # Try reading with index_col=0 and reset
    if bc_col is None or cl_col is None:
        eos_meta2 = pd.read_csv(meta_path, index_col=0).reset_index()
        bc_col = next((c for c in eos_meta2.columns if "barcode" in c.lower()), None)
        cl_col = next((c for c in eos_meta2.columns
                       if c.lower() in ["cluster", "index"]), None)
        if bc_col and cl_col:
            eos_meta = eos_meta2

    if bc_col and cl_col:
        eos_meta = eos_meta.set_index(bc_col)
        common   = auc_mtx.index.intersection(eos_meta.index)
        print(f"  Common cells: {len(common)}")
        if len(common) > 0:
            auc_mtx  = auc_mtx.loc[common]
            clusters = eos_meta.loc[common, cl_col].astype(str)
            clusters = clusters.str.replace("Eos_", "", regex=False)
    else:
        print(f"  WARNING: could not identify barcode/cluster columns")

if clusters is None or len(clusters) == 0:
    print("  Falling back to h5ad for cluster info...")
    h5ad_path = os.path.join(SCENIC, "eos_for_scenic.h5ad")
    adata_eos = sc.read_h5ad(h5ad_path)
    c_col = "cluster" if "cluster" in adata_eos.obs.columns else "eos_cluster"
    common    = auc_mtx.index.intersection(adata_eos.obs_names)
    auc_mtx   = auc_mtx.loc[common]
    clusters  = adata_eos.obs.loc[common, c_col].astype(str)
    clusters  = clusters.str.replace("Eos_", "", regex=False)
    del adata_eos

auc_mtx["cluster"] = clusters.values
print(f"  Clusters: {sorted(clusters.unique())}")

# =============================================================================
# Step 2: Parse regulon gene sets from regulons.csv
# =============================================================================
print("\nStep 2: Parsing regulon gene sets...")
regulon_gene_sets = {}

# Prefer adj.csv (GRN output: TF, target, importance) — clean format
adj_csv     = os.path.join(SCENIC, "adj.csv")
regulon_csv = os.path.join(SCENIC, "regulons.csv")

if os.path.exists(adj_csv):
    print(f"  Parsing adj.csv (GRN output)...")
    adj_df = pd.read_csv(adj_csv)
    print(f"  adj.csv columns: {list(adj_df.columns)}")
    tf_col  = next((c for c in adj_df.columns if c.lower() == "tf"), adj_df.columns[0])
    tgt_col = next((c for c in adj_df.columns if c.lower() in ["target","gene"]), adj_df.columns[1])
    for tf, grp in adj_df.groupby(tf_col):
        regulon_gene_sets[str(tf)] = list(grp[tgt_col].dropna().astype(str).unique())
    print(f"  Parsed {len(regulon_gene_sets)} TF gene sets from adj.csv")

elif os.path.exists(regulon_csv):
    print(f"  Falling back to regulons.csv...")
    # Try reading with second row as header (pySCENIC ctx multi-index format)
    try:
        reg_df = pd.read_csv(regulon_csv, header=[0, 1])
        # Look for TF and TargetGenes in multi-index columns
        tf_col  = [c for c in reg_df.columns if "TF" in str(c) or "tf" in str(c)]
        tgt_col = [c for c in reg_df.columns if "Target" in str(c) or "target" in str(c)]
        if tf_col and tgt_col:
            for _, row in reg_df.iterrows():
                tf = str(row[tf_col[0]]).split("(")[0].strip()
                tg_str = str(row[tgt_col[0]])
                try:
                    items = ast.literal_eval(tg_str)
                    genes = [i[0] if isinstance(i, (list, tuple)) else str(i) for i in items]
                except Exception:
                    genes = []
                if genes:
                    regulon_gene_sets.setdefault(tf, []).extend(genes)
    except Exception as e:
        print(f"  regulons.csv parse failed: {e}")
    print(f"  Parsed {len(regulon_gene_sets)} regulon gene sets")
else:
    print("  WARNING: neither adj.csv nor regulons.csv found")


def regulon_size(reg):
    """Number of inferred target genes for a regulon column name (handles
    the '(+)'/'(-)' suffix). Returns None if sizes are unknown (gene sets
    could not be parsed), in which case the size filter should be skipped
    rather than treated as size 0."""
    if not regulon_gene_sets:
        return None
    tf = reg.replace("(+)", "").replace("(-)", "").strip()
    return len(regulon_gene_sets.get(tf, []))


# =============================================================================
# Step 3: Compute Eos cluster means (subsampled)
# =============================================================================
print("\nStep 3: Computing Eos cluster mean AUCell (subsampled)...")
regulon_cols = [c for c in auc_mtx.columns if c != "cluster"]

np.random.seed(RANDOM_SEED)
eos_means = {}
for cl in EOS_CLUSTERS:
    mask = auc_mtx["cluster"] == cl
    sub  = auc_mtx.loc[mask, regulon_cols]
    n    = len(sub)
    if n == 0:
        print(f"  WARNING: no cells for cluster {cl}")
        continue
    idx = np.random.choice(n, min(n, CELLS_PER_GROUP),
                           replace=(n < CELLS_PER_GROUP))
    eos_means[cl] = sub.iloc[idx].mean()
    print(f"  {cl}: {n} cells")

eos_mean_df = pd.DataFrame(eos_means)  # regulons x clusters

# =============================================================================
# Step 4: Score Mac cells with regulon gene sets
# =============================================================================
print("\nStep 4: Scoring macrophage cells...")
mac_mean_df = None

mac_scored_csv = os.path.join(SCENIC, "mac_aucell_scored.csv")
mac_h5ad_path  = os.path.join(SCENIC, "mac_for_scenic_scoring.h5ad")
mac_meta_path  = os.path.join(BASE, "mac_subtyped_metadata.csv")

if os.path.exists(mac_scored_csv):
    print(f"  Loading pre-computed Mac scores: {mac_scored_csv}")
    mac_score_df = pd.read_csv(mac_scored_csv, index_col=0)
    mac_meta     = pd.read_csv(mac_meta_path, index_col=0)
    common_mac   = mac_meta.index.intersection(mac_score_df.index)
    mac_score_df = mac_score_df.loc[common_mac]
    subtype_col  = "mac_subtype" if "mac_subtype" in mac_meta.columns else mac_meta.columns[0]
    mac_score_df["subtype"] = mac_meta.loc[common_mac, subtype_col].map(
        lambda x: MAC_RENAME.get(x, x)).values
    mac_mean_df  = mac_score_df.groupby("subtype")[regulon_cols].mean()

elif os.path.exists(mac_h5ad_path) and len(regulon_gene_sets) > 0:
    print(f"  Loading Mac h5ad: {mac_h5ad_path}")
    mac_adata = sc.read_h5ad(mac_h5ad_path)
    print(f"  Mac: {mac_adata.n_obs} cells x {mac_adata.n_vars} genes")

    # Score each regulon using gene sets
    mac_score_df = pd.DataFrame(index=mac_adata.obs_names)
    scored = 0
    for reg in regulon_cols:
        tf = reg.replace("(+)", "").replace("(-)", "").strip()
        genes = regulon_gene_sets.get(tf, [])
        genes_ok = [g for g in genes if g in mac_adata.var_names]
        if len(genes_ok) < 3:
            mac_score_df[reg] = 0.0
            continue
        sc.tl.score_genes(mac_adata, gene_list=genes_ok,
                          score_name=f"_sc_{tf}", use_raw=False)
        mac_score_df[reg] = mac_adata.obs[f"_sc_{tf}"].values
        scored += 1
    print(f"  Scored {scored}/{len(regulon_cols)} regulons")

    mac_score_df.to_csv(mac_scored_csv)
    print(f"  Saved: {mac_scored_csv}")

    # Attach subtype labels
    if os.path.exists(mac_meta_path):
        mac_meta = pd.read_csv(mac_meta_path, index_col=0)
        subtype_col = "mac_subtype" if "mac_subtype" in mac_meta.columns else mac_meta.columns[0]
        common_mac  = mac_meta.index.intersection(mac_score_df.index)
        mac_score_df = mac_score_df.loc[common_mac]
        mac_score_df["subtype"] = mac_meta.loc[common_mac, subtype_col].map(
            lambda x: MAC_RENAME.get(x, x)).values
    elif "mac_subtype" in mac_adata.obs.columns:
        mac_score_df["subtype"] = mac_adata.obs["mac_subtype"].map(
            lambda x: MAC_RENAME.get(x, x)).values
    else:
        print("  WARNING: no mac subtype info found")
        mac_score_df["subtype"] = "Mac"

    mac_mean_df = mac_score_df.groupby("subtype")[regulon_cols].mean()
    del mac_adata

else:
    if not os.path.exists(mac_h5ad_path):
        print("  Mac h5ad not found — run 14b_export_mac_h5ad.R first")
    else:
        print("  No regulon gene sets — skipping Mac scoring")

if mac_mean_df is not None:
    # Reindex to standard subtype names
    present = [s for s in MAC_SUBTYPES if s in mac_mean_df.index]
    mac_mean_df = mac_mean_df.loc[present]
    print(f"  Mac subtypes available: {list(mac_mean_df.index)}")

# =============================================================================
# Step 5: Select top regulons per cluster using RSS (Regulon Specificity Score)
# =============================================================================
print("\nStep 5: Selecting top regulons per cluster via RSS...")

rss_path = os.path.join(BASE, "scenic", "regulon_specificity_scores.csv")
top_df = None

if os.path.exists(rss_path):
    rss = pd.read_csv(rss_path, index_col=0)
    print(f"  RSS loaded: {rss.shape[0]} regulons x {rss.shape[1]} clusters")

    # Normalize cluster column names (strip 'C' prefix if needed)
    rss.columns = [c.replace("Eos_", "") for c in rss.columns]

    top_rows = []
    for cl in EOS_CLUSTERS:
        if cl not in rss.columns:
            print(f"  WARNING: cluster {cl} not in RSS columns {list(rss.columns)}")
            continue
        # Exclude regulons with too few target genes to be interpretable
        # (skip the filter entirely if sizes are unknown, rather than
        # dropping everything)
        sizes = {r: regulon_size(r) for r in rss.index}
        if any(s is not None for s in sizes.values()):
            keep_idx = [r for r in rss.index if (sizes[r] or 0) >= MIN_REGULON_SIZE]
            valid = rss.loc[keep_idx]
        else:
            valid = rss
        top = valid[cl].nlargest(TOP_PER_CLUSTER)
        for reg, score in top.items():
            tf = reg.replace("(+)", "").replace("(-)", "").strip()
            top_rows.append({"regulon": reg, "tf": tf,
                             "max_cluster": cl, "rss": score})

    top_df = pd.DataFrame(top_rows).drop_duplicates("regulon")
    print(f"  Selected {len(top_df)} regulons ({TOP_PER_CLUSTER} per cluster)")
    for cl in EOS_CLUSTERS:
        sub = top_df[top_df["max_cluster"] == cl]
        if len(sub):
            print(f"    {cl}: {', '.join(sub['tf'].tolist())}")

else:
    # Fallback: uniqueness score from AUCell means
    print("  RSS not found — falling back to AUCell uniqueness")
    results = []
    for reg in regulon_cols:
        tf = reg.replace("(+)", "").replace("(-)", "").strip()
        size = regulon_size(reg)
        if (size is not None and size < MIN_REGULON_SIZE) or reg not in eos_mean_df.index:
            continue
        eos_vals = eos_mean_df.loc[reg, [c for c in EOS_CLUSTERS if c in eos_mean_df.columns]]
        if eos_vals.empty:
            continue
        sorted_eos = sorted(eos_vals.values, reverse=True)
        uniqueness = sorted_eos[0] - sorted_eos[1] if len(sorted_eos) >= 2 else sorted_eos[0]
        results.append({"regulon": reg, "tf": tf,
                        "max_cluster": eos_vals.idxmax(),
                        "max_auc": eos_vals.max(), "uniqueness": uniqueness})
    if not results:
        raise RuntimeError("No regulons found — check cluster labels")
    res_df = pd.DataFrame(results)
    top_rows = []
    for cl in EOS_CLUSTERS:
        sub = res_df[res_df["max_cluster"] == cl].nlargest(TOP_PER_CLUSTER, "uniqueness")
        top_rows.append(sub)
    top_df = pd.concat(top_rows, ignore_index=True).drop_duplicates("regulon")

# =============================================================================
# Step 6: Build heatmap matrix and plot
# =============================================================================
print("\nStep 6: Building combined heatmap...")

selected_regs = top_df["regulon"].tolist()

# Eos cluster means for selected regulons
eos_heat = eos_mean_df.loc[selected_regs,
                            [c for c in EOS_CLUSTERS if c in eos_mean_df.columns]].T

# Mac means for selected regulons (align to eos z-score scale)
if mac_mean_df is not None:
    mac_cols_avail = [r for r in selected_regs if r in mac_mean_df.columns]
    mac_heat = mac_mean_df.loc[
        [s for s in MAC_SUBTYPES if s in mac_mean_df.index],
        mac_cols_avail
    ]
    mac_heat.index = [MAC_DISPLAY.get(s, s) for s in mac_heat.index]
    # Fill missing regulons with 0
    for r in selected_regs:
        if r not in mac_heat.columns:
            mac_heat[r] = 0.0
    mac_heat = mac_heat[selected_regs]
else:
    mac_heat = None

# Z-score row-wise across ALL columns (Eos + Mac combined) for comparability
eos_z = eos_heat.copy()
eos_z.index = [str(c) for c in eos_z.index]

if mac_heat is not None:
    mac_z = mac_heat.copy()
    mac_z.index = [MAC_DISPLAY.get(s, s) for s in mac_z.index]
    combined_raw = pd.concat([eos_z, mac_z])  # cell_types x regulons
else:
    combined_raw = eos_z

# Row-wise z-score per regulon (transpose: regulons x cell_types, z-score per row)
mat = combined_raw.T  # regulons x cell_types
row_mean = mat.mean(axis=1)
row_std  = mat.std(axis=1).replace(0, 1)
combined_z = ((mat.T - row_mean) / row_std).T.clip(-ZSCORE_CAP, ZSCORE_CAP)
combined = combined_z.T  # back to cell_types x regulons

# Regulon labels with target gene count from adj.csv
def reg_label(reg):
    tf = reg.replace("(+)", "").replace("(-)", "").strip()
    n  = len(regulon_gene_sets.get(tf, []))
    return f"{tf}({n})" if n > 0 else tf

col_labels = [reg_label(r) for r in selected_regs]
plot_data  = combined.T  # regulons x cell_types
plot_data.index = col_labels

n_rows = len(col_labels)
n_cols = len(combined.index)
n_eos  = len(eos_z)

# Cell size matching reference image (~0.38" per row, 0.55" per col)
cell_h    = 0.38
cell_w    = 0.55
left_pad  = 1.8   # for y-axis labels
right_pad = 1.2   # for colorbar
top_pad   = 0.7   # for title
bot_pad   = 1.0   # for x-axis labels + group labels

fig_h = n_rows * cell_h + top_pad + bot_pad
fig_w = n_cols * cell_w + left_pad + right_pad

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

# Position axes precisely
ax_left   = left_pad / fig_w
ax_bottom = bot_pad  / fig_h
ax_width  = (n_cols * cell_w) / fig_w
ax_height = (n_rows * cell_h) / fig_h
ax.set_position([ax_left, ax_bottom, ax_width, ax_height])

im = ax.imshow(plot_data.values, cmap=ZSCORE_CMAP, aspect="auto",
               vmin=-ZSCORE_CAP, vmax=ZSCORE_CAP)

# X-axis: C0–C5, M1, M2, M2
x_labels = list(combined.index)
# Replace Mac display names with short versions for top tick
x_tick_labels = []
for lbl in x_labels:
    if lbl == "M1":
        x_tick_labels.append("M1")
    elif "Funct" in lbl or lbl == "M2_Functional":
        x_tick_labels.append("M2")
    elif "Surf" in lbl or lbl == "M2_Surface":
        x_tick_labels.append("M2")
    else:
        x_tick_labels.append(lbl)

ax.set_xticks(range(n_cols))
ax.set_xticklabels(x_tick_labels, fontsize=10, fontweight="bold")
ax.xaxis.set_ticks_position("bottom")

# Y-axis on LEFT, italic
ax.set_yticks(range(n_rows))
ax.set_yticklabels(col_labels, fontsize=9, ha="right")
ax.yaxis.set_ticks_position("left")
ax.tick_params(axis="both", length=0)

# Thick vertical separator between Eos and Mac
if mac_heat is not None and n_eos < n_cols:
    ax.axvline(x=n_eos - 0.5, color="black", linewidth=2.0)

# --- Group labels and underlines below x-axis (in figure coords) ---
# Transform: x in data coords → figure coords
def data_x_to_fig(x_data):
    return ax_left + (x_data + 0.5) / n_cols * ax_width

eos_x_center = data_x_to_fig((n_eos - 1) / 2)
mac_x_center = data_x_to_fig(n_eos + (n_cols - n_eos - 1) / 2)

line_y  = ax_bottom - 0.35 / fig_h   # just below x-tick labels
label_y = ax_bottom - 0.65 / fig_h

eos_x0 = data_x_to_fig(-0.5)
eos_x1 = data_x_to_fig(n_eos - 0.5)
mac_x0 = data_x_to_fig(n_eos - 0.5)
mac_x1 = data_x_to_fig(n_cols - 0.5)

# Eosinophils underline
fig.add_artist(plt.Line2D([eos_x0, eos_x1], [line_y, line_y],
                           transform=fig.transFigure,
                           color="black", linewidth=1.2, clip_on=False))
fig.text(eos_x_center, label_y, "Eosinophils",
         ha="center", va="top", fontsize=10, transform=fig.transFigure)

if mac_heat is not None and n_eos < n_cols:
    # Mac underline
    fig.add_artist(plt.Line2D([mac_x0, mac_x1], [line_y, line_y],
                               transform=fig.transFigure,
                               color="black", linewidth=1.2, clip_on=False))
    fig.text(mac_x_center, label_y, "Mac",
             ha="center", va="top", fontsize=10, transform=fig.transFigure)

    # "Funct." and "Surf." sub-labels
    sublabel_y = ax_bottom - 0.22 / fig_h
    if n_cols >= n_eos + 3:
        fig.text(data_x_to_fig(n_eos + 1), sublabel_y, "Funct.",
                 ha="center", va="top", fontsize=8, transform=fig.transFigure)
        fig.text(data_x_to_fig(n_eos + 2), sublabel_y, "Surf.",
                 ha="center", va="top", fontsize=8, transform=fig.transFigure)

# Title
ax.set_title("Eos Subcluster Unique Regulons", fontsize=13,
             fontweight="bold", pad=10)

# Colorbar
cbar_ax = fig.add_axes([ax_left + ax_width + 0.05 / fig_w,
                         ax_bottom + ax_height * 0.2,
                         0.03 / fig_w * fig_h,   # width proportional
                         ax_height * 0.6])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("Z-score", fontsize=9)
cbar.set_ticks([-2, -1, 0, 1, 2])
cbar.ax.tick_params(labelsize=8)

out_pdf = os.path.join(FIG_DIR, "Fig5_regulon_heatmap_eos_mac.pdf")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_pdf}")

# Save selection table
top_df.to_csv(os.path.join(SCENIC, "eos_unique_regulons_selected.csv"), index=False)
print(f"  Saved: {os.path.join(SCENIC, 'eos_unique_regulons_selected.csv')}")

print("\n16b complete.")
