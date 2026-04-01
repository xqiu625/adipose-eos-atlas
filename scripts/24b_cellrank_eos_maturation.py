#!/usr/bin/env python
"""
24b_cellrank_eos_maturation.py
CellRank trajectory on three-way integrated eosinophils:
  Gurtner BM (progenitor) → Blood → Hu healthy adipose → Qiu HFD adipose

Uses CytoTRACE kernel to infer maturation trajectory without RNA velocity.
Root: Gurtner bone marrow (most progenitor-like)
Terminal: Qiu disease-specific clusters (C0, C3, C4)

Input:  $BASE/integration/three_way/for_cellrank/
Output: $BASE/integration/three_way/cellrank_*.pdf, .csv, .h5ad
SLURM:  epyc, 64G, 8 CPUs, 4h, cellrank_env
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
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# --- Color palettes -----------------------------------------------------------
DATASET_COLORS = {
    "Qiu_2025_Adipose": "#E63946",
    "Hu_2026_Healthy": "#4DAF4A",
    "Gurtner_2022_Blood_BM": "#1D3557",
}

EOS_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4", "C6": "#2F3C7E",
}

ORGAN_COLORS = {
    "Adipose_HFD": "#E63946",
    "Adipose_Healthy": "#4DAF4A",
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
EXPORT_DIR = os.path.join(BASE, "integration", "three_way", "for_cellrank")
OUT_DIR = os.path.join(BASE, "integration", "three_way")
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# --- Load exported data from script 24 ---------------------------------------
print("Loading exported data from three-way integration ...")

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
print(f"  Organs: {adata.obs['organ'].value_counts().to_dict()}")

# --- Preprocessing for CytoTRACE --------------------------------------------
print("\nPreprocessing ...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Use Harmony embeddings for neighbors
sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=30)

# CytoTRACE needs 'Ms' layer (moment-smoothed expression)
# Without RNA velocity data, use normalized log-transformed expression
adata.layers["Ms"] = adata.X.copy()

# --- CytoTRACE Kernel -------------------------------------------------------
print("Computing CytoTRACE kernel ...")
ctk = cr.kernels.CytoTRACEKernel(adata)
ctk.compute_cytotrace()
ctk.compute_transition_matrix()

print(f"  CytoTRACE score range: {adata.obs['ct_score'].min():.3f} - {adata.obs['ct_score'].max():.3f}")

# CytoTRACE by dataset summary
for ds in adata.obs["dataset"].unique():
    mask = adata.obs["dataset"] == ds
    scores = adata.obs.loc[mask, "ct_score"]
    print(f"  {ds}: mean={scores.mean():.3f}, median={scores.median():.3f}")

# --- Ensure categorical columns for CellRank ---------------------------------
# Create unique cluster labels (dataset prefix) to avoid duplicate names across datasets
ds_short = {"Qiu_2025_Adipose": "Qiu", "Hu_2026_Healthy": "Hu", "Gurtner_2022_Blood_BM": "Gurtner"}
adata.obs["dataset_cluster"] = (
    adata.obs["dataset"].map(ds_short).astype(str) + "_" +
    adata.obs["original_cluster"].astype(str)
)
for col in ["original_cluster", "dataset", "organ", "study", "dataset_cluster"]:
    if col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")

# --- GPCCA estimator --------------------------------------------------------
print("\nComputing macrostates with GPCCA ...")
g = cr.estimators.GPCCA(ctk)
g.compute_schur(n_components=20)
g.compute_macrostates(n_states=5, cluster_key="dataset_cluster")

print(f"  Macrostates: {g.macrostates.cat.categories.tolist()}")

g.predict_initial_states()
g.predict_terminal_states(allow_overlap=True)

print(f"  Initial states: {g.initial_states.cat.categories.tolist()}")
print(f"  Terminal states: {g.terminal_states.cat.categories.tolist()}")

# --- Fate probabilities ------------------------------------------------------
print("Computing fate probabilities ...")
g.compute_fate_probabilities()

# --- Diffusion pseudotime from BM root --------------------------------------
print("\nComputing diffusion pseudotime rooted on Gurtner BM cells ...")

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

# =============================================================================
# Figures
# =============================================================================
print("\nGenerating figures ...")

# --- Fig 1: CytoTRACE + Dataset + Organ on UMAP ---
fig, axes = plt.subplots(1, 3, figsize=(21, 6))
sc.pl.embedding(adata, basis="umap", color="ct_score", cmap="viridis",
                title="CytoTRACE Score\n(high = progenitor-like)", ax=axes[0], show=False)
sc.pl.embedding(adata, basis="umap", color="dataset", palette=DATASET_COLORS,
                title="Dataset", ax=axes[1], show=False)
sc.pl.embedding(adata, basis="umap", color="organ",
                title="Organ / Tissue", ax=axes[2], show=False)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "cellrank_cytotrace_overview.pdf"))
plt.close()

# --- Fig 2: CytoTRACE boxplots by dataset and organ ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# By dataset
datasets = sorted(adata.obs["dataset"].unique())
data_by_ds = [adata.obs.loc[adata.obs["dataset"] == ds, "ct_score"].values for ds in datasets]
bp1 = axes[0].boxplot(data_by_ds, labels=[d.split("_")[0] for d in datasets], patch_artist=True)
for patch, color in zip(bp1["boxes"], [DATASET_COLORS.get(d, "#999") for d in datasets]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
axes[0].set_ylabel("CytoTRACE Score")
axes[0].set_title("CytoTRACE by Dataset")

# By organ (top organs)
organ_counts = adata.obs["organ"].value_counts()
top_organs = organ_counts[organ_counts >= 50].index.tolist()
data_by_organ = [adata.obs.loc[adata.obs["organ"] == o, "ct_score"].values for o in top_organs]
bp2 = axes[1].boxplot(data_by_organ, labels=top_organs, patch_artist=True)
for patch, organ in zip(bp2["boxes"], top_organs):
    patch.set_facecolor(ORGAN_COLORS.get(organ, "#999"))
    patch.set_alpha(0.7)
axes[1].set_ylabel("CytoTRACE Score")
axes[1].set_title("CytoTRACE by Organ")
plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=45, ha="right")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "cellrank_cytotrace_boxplots.pdf"))
plt.close()

# --- Fig 3: CytoTRACE by Qiu cluster ---
qiu_mask = adata.obs["dataset"] == "Qiu_2025_Adipose"
if qiu_mask.sum() > 0:
    fig, ax = plt.subplots(figsize=(8, 5))
    qiu_data = adata.obs.loc[qiu_mask].copy()
    clusters = sorted(qiu_data["original_cluster"].unique())
    data_by_cl = [qiu_data.loc[qiu_data["original_cluster"] == c, "ct_score"].values for c in clusters]
    bp = ax.boxplot(data_by_cl, labels=clusters, patch_artist=True)
    for patch, cl in zip(bp["boxes"], clusters):
        cl_short = cl.replace("Eos_", "")
        patch.set_facecolor(EOS_COLORS.get(cl_short, "#999"))
        patch.set_alpha(0.7)
    ax.set_ylabel("CytoTRACE Score")
    ax.set_title("CytoTRACE by Eos Cluster (Qiu)")
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "cellrank_cytotrace_by_cluster.pdf"))
    plt.close()

# --- Fig 4: Macrostates on UMAP ---
fig, ax = plt.subplots(figsize=(10, 7))
g.plot_macrostates(which="all", basis="umap", legend_loc="right", ax=ax, show=False)
ax.set_title("CellRank Macrostates (Three-way)")
fig.savefig(os.path.join(FIG_DIR, "cellrank_macrostates.pdf"))
plt.close()

# --- Fig 5: Initial and terminal states ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
try:
    g.plot_macrostates(which="initial", basis="umap", ax=axes[0], show=False)
    axes[0].set_title("Initial States")
except Exception as e:
    print(f"  Initial states plot failed: {e}")
try:
    g.plot_macrostates(which="terminal", basis="umap", ax=axes[1], show=False)
    axes[1].set_title("Terminal States")
except Exception as e:
    print(f"  Terminal states plot failed: {e}")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "cellrank_initial_terminal.pdf"))
plt.close()

# --- Fig 6: Fate probabilities ---
try:
    g.plot_fate_probabilities(basis="umap", same_plot=False)
    plt.savefig(os.path.join(FIG_DIR, "cellrank_fate_probabilities.pdf"))
    plt.close()
except Exception as e:
    print(f"  Fate probability plot failed: {e}")

# --- Fig 7: DPT pseudotime ---
if has_dpt:
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    sc.pl.embedding(adata, basis="umap", color="dpt_pseudotime", cmap="viridis",
                    title="Diffusion Pseudotime\n(root = Gurtner BM)",
                    ax=axes[0], show=False)
    sc.pl.embedding(adata, basis="umap", color="dataset", palette=DATASET_COLORS,
                    title="Dataset", ax=axes[1], show=False)
    sc.pl.embedding(adata, basis="umap", color="original_cluster",
                    title="Original Clusters", ax=axes[2], show=False)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "cellrank_dpt_pseudotime.pdf"))
    plt.close()

    # DPT by dataset boxplot
    fig, ax = plt.subplots(figsize=(8, 5))
    dpt_data = adata.obs[["dpt_pseudotime", "dataset"]].copy()
    dpt_data = dpt_data[dpt_data["dpt_pseudotime"] < 1.0]  # remove inf
    ds_list = sorted(dpt_data["dataset"].unique())
    data_dpt = [dpt_data.loc[dpt_data["dataset"] == d, "dpt_pseudotime"].values for d in ds_list]
    bp = ax.boxplot(data_dpt, labels=[d.split("_")[0] for d in ds_list], patch_artist=True)
    for patch, d in zip(bp["boxes"], ds_list):
        patch.set_facecolor(DATASET_COLORS.get(d, "#999"))
        patch.set_alpha(0.7)
    ax.set_ylabel("Diffusion Pseudotime")
    ax.set_title("DPT by Dataset\n(BM = 0, Disease adipose = 1)")
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "cellrank_dpt_by_dataset.pdf"))
    plt.close()

    # DPT by Qiu cluster
    if qiu_mask.sum() > 0:
        fig, ax = plt.subplots(figsize=(8, 5))
        qiu_dpt = adata.obs.loc[qiu_mask, ["dpt_pseudotime", "original_cluster"]].copy()
        qiu_dpt = qiu_dpt[qiu_dpt["dpt_pseudotime"] < 1.0]
        clusters = sorted(qiu_dpt["original_cluster"].unique())
        data_cl = [qiu_dpt.loc[qiu_dpt["original_cluster"] == c, "dpt_pseudotime"].values
                    for c in clusters]
        bp = ax.boxplot(data_cl, labels=clusters, patch_artist=True)
        for patch, cl in zip(bp["boxes"], clusters):
            cl_short = cl.replace("Eos_", "")
            patch.set_facecolor(EOS_COLORS.get(cl_short, "#999"))
            patch.set_alpha(0.7)
        ax.set_ylabel("Diffusion Pseudotime")
        ax.set_title("DPT by Eos Cluster (Qiu)\n(BM root → tissue maturation)")
        plt.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "cellrank_dpt_by_cluster.pdf"))
        plt.close()

# =============================================================================
# Save results
# =============================================================================
print("\nSaving results ...")

# Fate probabilities
if g.fate_probabilities is not None:
    fate_df = pd.DataFrame(
        g.fate_probabilities.X,
        index=adata.obs_names,
        columns=g.fate_probabilities.names
    )
    fate_df.to_csv(os.path.join(OUT_DIR, "cellrank_fate_probabilities.csv"))

# CytoTRACE scores + DPT
out_cols = ["dataset", "organ", "original_cluster", "ct_score"]
if has_dpt:
    out_cols.append("dpt_pseudotime")
if "condition" in adata.obs.columns:
    out_cols.append("condition")
adata.obs[out_cols].to_csv(os.path.join(OUT_DIR, "cellrank_scores.csv"))

# Save AnnData — convert bool/object columns to strings for h5py compatibility
for col in adata.obs.columns:
    if adata.obs[col].dtype == object or adata.obs[col].dtype == bool:
        adata.obs[col] = adata.obs[col].astype(str)
adata.write_h5ad(os.path.join(OUT_DIR, "three_way_cellrank.h5ad"))

print(f"\n========== Three-way CellRank analysis complete ==========")
print(f"Output directory: {OUT_DIR}")
