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

EXCLUDE_REGULONS = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Tef", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

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
if os.path.exists(meta_path):
    eos_meta = pd.read_csv(meta_path, index_col=0)
    cluster_col = "cluster" if "cluster" in eos_meta.columns else eos_meta.columns[0]
    common = auc_mtx.index.intersection(eos_meta.index)
    auc_mtx  = auc_mtx.loc[common]
    clusters = eos_meta.loc[common, cluster_col].astype(str)
    # Normalize cluster names to C0–C5
    clusters = clusters.str.replace("Eos_", "", regex=False)
else:
    # Fall back: load h5ad
    h5ad_path = os.path.join(SCENIC, "eos_for_scenic.h5ad")
    adata_eos = sc.read_h5ad(h5ad_path)
    c_col = "cluster" if "cluster" in adata_eos.obs.columns else "eos_cluster"
    common = auc_mtx.index.intersection(adata_eos.obs_names)
    auc_mtx  = auc_mtx.loc[common]
    clusters = adata_eos.obs.loc[common, c_col].astype(str)
    clusters = clusters.str.replace("Eos_", "", regex=False)

auc_mtx["cluster"] = clusters.values
print(f"  Clusters: {sorted(clusters.unique())}")

# =============================================================================
# Step 2: Parse regulon gene sets from regulons.csv
# =============================================================================
print("\nStep 2: Parsing regulon gene sets...")
regulon_gene_sets = {}
regulon_csv = os.path.join(SCENIC, "regulons.csv")
if os.path.exists(regulon_csv):
    reg_df = pd.read_csv(regulon_csv)
    print(f"  regulons.csv columns: {list(reg_df.columns)}")

    tf_col = target_col = None
    for c in reg_df.columns:
        cl = c.lower()
        if cl in ["tf", "regulon", "transcription_factor"] and tf_col is None:
            tf_col = c
        if cl in ["target", "gene", "target_gene", "targetgenes"] and target_col is None:
            target_col = c

    if tf_col and target_col:
        if target_col.lower() in ["targetgenes", "target_genes"]:
            # pySCENIC ctx format: TargetGenes = "[(gene, score), ...]"
            for _, row in reg_df.iterrows():
                tf = str(row[tf_col]).split("(")[0].strip()
                tg_str = str(row[target_col])
                try:
                    items = ast.literal_eval(tg_str)
                    genes = [i[0] if isinstance(i, (list, tuple)) else str(i)
                             for i in items]
                except Exception:
                    genes = [g.strip() for g in tg_str.replace("[", "").replace("]", "")
                             .split(",") if g.strip()]
                if tf not in regulon_gene_sets:
                    regulon_gene_sets[tf] = []
                regulon_gene_sets[tf].extend(genes)
        else:
            # Long format: one row per TF-target pair
            for tf, grp in reg_df.groupby(tf_col):
                tf_name = str(tf).split("(")[0].strip()
                regulon_gene_sets[tf_name] = list(grp[target_col].dropna().unique())

    print(f"  Parsed {len(regulon_gene_sets)} regulon gene sets")
else:
    print("  WARNING: regulons.csv not found — Mac scoring will be skipped")

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
# Step 5: Find Eos-unique regulons
# =============================================================================
print("\nStep 5: Selecting Eos-unique regulons...")

results = []
for reg in regulon_cols:
    if reg not in eos_mean_df.index:
        continue
    tf = reg.replace("(+)", "").replace("(-)", "").strip()
    if tf in EXCLUDE_REGULONS:
        continue

    eos_vals = eos_mean_df.loc[reg, [c for c in EOS_CLUSTERS if c in eos_mean_df.columns]]
    if eos_vals.empty:
        continue
    max_eos_cluster = eos_vals.idxmax()
    max_eos_auc     = eos_vals.max()

    if mac_mean_df is not None and reg in mac_mean_df.columns:
        mac_vals = mac_mean_df[reg]
        mac_mean = mac_vals.mean()
    else:
        mac_mean = 0.0

    fold_change = max_eos_auc / mac_mean if mac_mean > 0 else 0
    sorted_eos  = sorted(eos_vals.values, reverse=True)
    uniqueness  = sorted_eos[0] - sorted_eos[1] if len(sorted_eos) >= 2 else sorted_eos[0]

    results.append({
        "regulon": reg, "tf": tf,
        "max_cluster": max_eos_cluster,
        "max_auc": max_eos_auc,
        "fold_change": fold_change,
        "uniqueness": uniqueness,
    })

res_df = pd.DataFrame(results)
if mac_mean_df is not None:
    eos_unique = res_df[
        (res_df["max_auc"] >= MIN_EOS_AUC) &
        (res_df["fold_change"] >= MIN_FOLD_CHANGE) &
        (res_df["uniqueness"] >= MIN_UNIQUENESS)
    ]
else:
    # No Mac data: select top N per cluster by uniqueness
    eos_unique = res_df[res_df["max_auc"] >= MIN_EOS_AUC]

# Top N per Eos cluster
top_rows = []
for cl in EOS_CLUSTERS:
    sub = eos_unique[eos_unique["max_cluster"] == cl].nlargest(TOP_PER_CLUSTER, "uniqueness")
    top_rows.append(sub)
top_df = pd.concat(top_rows, ignore_index=True).drop_duplicates("regulon")

print(f"  Selected {len(top_df)} Eos-unique regulons "
      f"({TOP_PER_CLUSTER} per cluster)")
for cl in EOS_CLUSTERS:
    sub = top_df[top_df["max_cluster"] == cl]
    if len(sub):
        print(f"    {cl}: {', '.join(sub['tf'].tolist())}")

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

# Z-score: compute per-regulon mean/std from ALL Eos cells (not just sampled)
reg_data = auc_mtx[selected_regs]
z_mean = reg_data.mean()
z_std  = reg_data.std().replace(0, 1)

eos_z = ((eos_heat - z_mean) / z_std).clip(-ZSCORE_CAP, ZSCORE_CAP)
eos_z.index = [f"C{c}" if not str(c).startswith("C") else str(c) for c in eos_z.index]

if mac_heat is not None:
    mac_z = ((mac_heat - z_mean) / z_std).clip(-ZSCORE_CAP, ZSCORE_CAP)
    combined = pd.concat([eos_z, mac_z])
else:
    combined = eos_z

# Regulon labels with target gene count
def reg_label(reg):
    tf = reg.replace("(+)", "").replace("(-)", "").strip()
    n  = len(regulon_gene_sets.get(tf, []))
    return f"{tf}({n})" if n > 0 else tf

col_labels = [reg_label(r) for r in selected_regs]
plot_data  = combined.T  # regulons x cell_types
plot_data.index = col_labels

n_rows = len(col_labels)
n_cols = len(combined.index)
fig_h  = max(5, n_rows * 0.30 + 2)
fig_w  = max(4, n_cols * 0.55 + 2.5)

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

im = ax.imshow(plot_data.values, cmap=ZSCORE_CMAP, aspect="auto",
               vmin=-ZSCORE_CAP, vmax=ZSCORE_CAP)

# X-axis (cell types)
ax.set_xticks(range(n_cols))
ax.set_xticklabels(list(combined.index), fontsize=8, rotation=0)

# Y-axis (regulons)
ax.set_yticks(range(n_rows))
ax.set_yticklabels(col_labels, fontsize=7, fontstyle="italic")

# Separator line between Eos and Mac
n_eos = len(eos_z)
if mac_heat is not None and n_eos < n_cols:
    ax.axvline(x=n_eos - 0.5, color="black", linewidth=1.5)
    # Group labels
    ax.text((n_eos - 1) / 2, -1.2, "Eosinophils",
            ha="center", va="top", fontsize=9, fontweight="bold",
            transform=ax.get_xaxis_transform())
    ax.text(n_eos + (n_cols - n_eos - 1) / 2, -1.2, "Mac",
            ha="center", va="top", fontsize=9, fontweight="bold",
            transform=ax.get_xaxis_transform())

# Horizontal separators between Eos clusters
row_idx = 0
for cl in EOS_CLUSTERS:
    n_cl = len(top_df[top_df["max_cluster"] == cl])
    row_idx += n_cl
    if row_idx < n_rows:
        ax.axhline(y=row_idx - 0.5, color="white", linewidth=1.2)

ax.set_title("Eos Subcluster Unique Regulons", fontsize=11, fontweight="bold", pad=12)
ax.tick_params(length=0)

# Colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02, shrink=0.6)
cbar.set_label("Z-score", fontsize=8)
cbar.ax.tick_params(labelsize=7)

plt.tight_layout()

out_pdf = os.path.join(FIG_DIR, "Fig5_regulon_heatmap_eos_mac.pdf")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_pdf}")

# Save selection table
top_df.to_csv(os.path.join(SCENIC, "eos_unique_regulons_selected.csv"), index=False)
print(f"  Saved: {os.path.join(SCENIC, 'eos_unique_regulons_selected.csv')}")

print("\n16b complete.")
