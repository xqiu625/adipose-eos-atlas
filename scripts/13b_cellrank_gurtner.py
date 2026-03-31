#!/usr/bin/env python
"""
13b_cellrank_gurtner.py
CellRank trajectory analysis on Gurtner-integrated eosinophils
Uses CytoTRACE kernel (no RNA velocity needed for 10x Flex)
Root: Gurtner bone marrow eosinophils (most immature)

Input:  $BASE/integration/for_cellrank/ (from script 13)
Output: $BASE/integration/cellrank_*.pdf, cellrank_fate_probs.csv
SLURM:  batch, 64G, 8 CPUs, 4h, cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import mmread
from scipy.sparse import csc_matrix

import scanpy as sc
import cellrank as cr

# --- Publication figure settings ---------------------------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
    "pdf.fonttype": 42,       # TrueType
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# --- Color palettes (match 00_theme.R) ---------------------------------------
EOS_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4", "C6": "#2F3C7E",
}
DATASET_COLORS = {
    "Qiu_2025_Adipose": "#E63946",
    "Gurtner_2022": "#1D3557",
}
ORGAN_COLORS = {
    "Adipose": "#E63946",
    "Bone Marrow": "#1D3557",
    "Blood": "#457B9D",
    "Colon": "#2A9D8F",
    "Small Intestine": "#E9C46A",
    "Spleen": "#F4A261",
    "Stomach": "#264653",
    "Lung": "#A8DADC",
    "Uterus": "#CDB4DB",
    "Unknown": "#999999",
}

# --- Paths -------------------------------------------------------------------
BASE = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
EXPORT_DIR = os.path.join(BASE, "integration", "for_cellrank")
OUT_DIR = os.path.join(BASE, "integration")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# --- Load data exported from script 13 --------------------------------------
print("Loading exported data from script 13 ...")

counts = mmread(os.path.join(EXPORT_DIR, "counts.mtx")).T.tocsc()
genes = pd.read_csv(os.path.join(EXPORT_DIR, "genes.txt"), header=None)[0].values
barcodes = pd.read_csv(os.path.join(EXPORT_DIR, "barcodes.txt"), header=None)[0].values
metadata = pd.read_csv(os.path.join(EXPORT_DIR, "metadata.csv"), index_col=0)
harmony = pd.read_csv(os.path.join(EXPORT_DIR, "harmony_embeddings.csv"), index_col=0)
umap = pd.read_csv(os.path.join(EXPORT_DIR, "umap_embeddings.csv"), index_col=0)

# Build AnnData
adata = sc.AnnData(X=csc_matrix(counts))
adata.var_names = genes.astype(str)
adata.obs_names = barcodes.astype(str)
adata.obs = metadata.loc[adata.obs_names].copy()
adata.obsm["X_harmony"] = harmony.loc[adata.obs_names].values
adata.obsm["X_umap"] = umap.loc[adata.obs_names].values

print(f"  AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  Datasets: {adata.obs['dataset'].value_counts().to_dict()}")

# --- Preprocessing for CytoTRACE --------------------------------------------
print("Preprocessing ...")
# CytoTRACE needs normalized + log-transformed data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Use Harmony embeddings for neighbors
sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=30)

# CytoTRACE needs 'Ms' layer (moment-smoothed expression)
# Without RNA velocity data, use normalized log-transformed expression
adata.layers["Ms"] = adata.X.copy()

# --- CytoTRACE Kernel -------------------------------------------------------
# CytoTRACE estimates differentiation potential based on gene count
# Higher CytoTRACE score = more progenitor-like
print("Computing CytoTRACE kernel ...")
ctk = cr.kernels.CytoTRACEKernel(adata)
ctk.compute_cytotrace()
ctk.compute_transition_matrix()

print(f"  CytoTRACE score range: {adata.obs['ct_score'].min():.3f} - {adata.obs['ct_score'].max():.3f}")

# --- Ensure categorical columns for CellRank ---------------------------------
for col in ["original_cluster", "dataset", "organ"]:
    if col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")

# --- Identify initial and terminal states ------------------------------------
print("Computing macrostates with GPCCA ...")
g = cr.estimators.GPCCA(ctk)
g.compute_schur(n_components=20)
g.compute_macrostates(n_states=6, cluster_key="original_cluster")

print(f"  Macrostates: {g.macrostates.cat.categories.tolist()}")

# Identify initial states: should overlap with Gurtner BM cells
# Gurtner cells are more progenitor-like (higher CytoTRACE), so they should
# be detected as initial states automatically
g.predict_initial_states()
g.predict_terminal_states(allow_overlap=True)

print(f"  Initial states: {g.initial_states.cat.categories.tolist()}")
print(f"  Terminal states: {g.terminal_states.cat.categories.tolist()}")

# --- Compute fate probabilities ----------------------------------------------
print("Computing fate probabilities ...")
g.compute_fate_probabilities()

# --- Compute DPT (diffusion pseudotime) from BM root -----------------------
print("Computing diffusion pseudotime rooted on Gurtner bone marrow cells ...")

# Find root cell: Bone Marrow cell with highest CytoTRACE score (most immature)
bm_mask = adata.obs["organ"] == "Bone Marrow"
print(f"  Bone Marrow cells: {bm_mask.sum()}")
if bm_mask.sum() > 0:
    bm_ct = adata.obs.loc[bm_mask, "ct_score"]
    root_cell_idx = adata.obs_names.get_loc(bm_ct.idxmax())
    adata.uns["iroot"] = root_cell_idx
    sc.tl.diffmap(adata, n_comps=15)
    sc.tl.dpt(adata)
    print(f"  DPT range: {adata.obs['dpt_pseudotime'].min():.3f} - {adata.obs['dpt_pseudotime'].max():.3f}")
    has_dpt = True
else:
    print("  WARNING: No Bone Marrow cells found. Skipping DPT.")
    has_dpt = False

# --- Figures -----------------------------------------------------------------
print("Generating figures ...")

# Figure 1: CytoTRACE score on UMAP + organ labels
fig, axes = plt.subplots(1, 3, figsize=(20, 6))
sc.pl.embedding(adata, basis="umap", color="ct_score", cmap="viridis",
                title="CytoTRACE Score\n(high = progenitor-like)", ax=axes[0], show=False)
sc.pl.embedding(adata, basis="umap", color="organ", palette=ORGAN_COLORS,
                title="Organ / Tissue", ax=axes[1], show=False)
sc.pl.embedding(adata, basis="umap", color="dataset", palette=DATASET_COLORS,
                title="Dataset", ax=axes[2], show=False)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "13b_cytotrace_umap.pdf"))
plt.close()

# Figure 2: CytoTRACE score by dataset and cluster
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
# By dataset
df_ct = adata.obs[["ct_score", "dataset"]].copy()
df_ct.boxplot(column="ct_score", by="dataset", ax=axes[0])
axes[0].set_title("CytoTRACE by Dataset")
axes[0].set_xlabel("")
axes[0].set_ylabel("CytoTRACE Score")
plt.sca(axes[0])
plt.xticks(rotation=45, ha="right")

# By original cluster
df_ct2 = adata.obs[["ct_score", "original_cluster", "dataset"]].copy()
qiu_only = df_ct2[df_ct2["dataset"] == "Qiu_2025_Adipose"]
qiu_only.boxplot(column="ct_score", by="original_cluster", ax=axes[1])
axes[1].set_title("CytoTRACE by Eos Cluster (Qiu only)")
axes[1].set_xlabel("")
axes[1].set_ylabel("CytoTRACE Score")
plt.sca(axes[1])
plt.xticks(rotation=45, ha="right")
plt.suptitle("")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "13b_cytotrace_boxplots.pdf"))
plt.close()

# Figure 3: Macrostates on UMAP
fig, ax = plt.subplots(figsize=(8, 6))
g.plot_macrostates(which="all", basis="umap", legend_loc="right", ax=ax, show=False)
ax.set_title("CellRank Macrostates")
fig.savefig(os.path.join(FIG_DIR, "13b_cellrank_macrostates.pdf"))
plt.close()

# Figure 4: Fate probabilities on UMAP
try:
    g.plot_fate_probabilities(basis="umap", same_plot=False)
    plt.savefig(os.path.join(FIG_DIR, "13b_cellrank_fate_probabilities.pdf"))
    plt.close()
except Exception as e:
    print(f"  Warning: fate probability plot failed: {e}")

# Figure 5: DPT pseudotime on UMAP (if available)
if has_dpt:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    sc.pl.embedding(adata, basis="umap", color="dpt_pseudotime", cmap="viridis",
                    title="Diffusion Pseudotime\n(root = Gurtner BM)", ax=axes[0], show=False)
    sc.pl.embedding(adata, basis="umap", color="original_cluster",
                    title="Original Clusters", ax=axes[1], show=False)
    sc.pl.embedding(adata, basis="umap", color="dataset", palette=DATASET_COLORS,
                    title="Dataset", ax=axes[2], show=False)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "13b_dpt_pseudotime_umap.pdf"))
    plt.close()

    # Figure 6: DPT boxplot by condition (Qiu cells only)
    qiu_mask = adata.obs["dataset"] == "Qiu_2025_Adipose"
    if "condition" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 5))
        qiu_dpt = adata.obs.loc[qiu_mask, ["dpt_pseudotime", "condition"]].copy()
        qiu_dpt = qiu_dpt[qiu_dpt["dpt_pseudotime"] < 1.0]  # remove inf
        conditions = sorted(qiu_dpt["condition"].unique())
        data_by_cond = [qiu_dpt.loc[qiu_dpt["condition"] == c, "dpt_pseudotime"].values
                        for c in conditions]
        bp = ax.boxplot(data_by_cond, labels=conditions, patch_artist=True)
        ax.set_ylabel("Diffusion Pseudotime")
        ax.set_title("DPT by Condition (Qiu adipose eos)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "13b_dpt_by_condition.pdf"))
        plt.close()

# --- Save results ------------------------------------------------------------
print("Saving results ...")

# Fate probabilities
if g.fate_probabilities is not None:
    fate_df = pd.DataFrame(
        g.fate_probabilities.X,
        index=adata.obs_names,
        columns=g.fate_probabilities.names
    )
    fate_df.to_csv(os.path.join(OUT_DIR, "cellrank_fate_probabilities.csv"))

# CytoTRACE scores + DPT
out_df = adata.obs[["dataset", "original_cluster", "ct_score"]].copy()
if has_dpt:
    out_df["dpt_pseudotime"] = adata.obs["dpt_pseudotime"]
out_df.to_csv(os.path.join(OUT_DIR, "cellrank_cytotrace_scores.csv"))

# Save AnnData
adata.write_h5ad(os.path.join(OUT_DIR, "qiu_gurtner_cellrank.h5ad"))

print(f"\n========== CellRank analysis complete ==========")
print(f"Output directory: {OUT_DIR}")
