#!/usr/bin/env python
# ==============================================================================
# 62b_foundation_combat_recluster.py   (Integration+decont redo — STAGE 2)
# On the decontaminated, C0-trimmed eosinophils (62a):
#   - lognorm the decontaminated counts; batch-aware HVGs
#   - ComBat batch correction on gem_reaction (the chosen primary integration)
#   - PCA on the ComBat-corrected matrix; neighbours; Leiden tuned to ~6 clusters
#   - UMAP; map new clusters to the original C0-C5 by max overlap (Jaccard)
#   - VERIFICATION GATE: ARI / overlap of new clustering vs original C0-C5
# Saves the foundation object that all downstream re-analyses read.
# Run in sc_scvi env.
#   X            = raw counts ; layers[decontX]=decont ; layers[lognorm]=lognorm decont
#   obs[new_cluster] (new Leiden), obs[new_C] (mapped C0-C5), obsm[X_combat], obsm[X_umap]
# ==============================================================================
import argparse, os, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--n_hvg", type=int, default=2000)
ap.add_argument("--n_pcs", type=int, default=30)
ap.add_argument("--target_k", type=int, default=6)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

ad = sc.read_h5ad(a.h5ad)
ad.obs["gem_reaction"] = ad.obs["gem_reaction"].astype("category")
ad.layers["counts"] = ad.X.copy()
print(f"loaded {ad.n_obs} cells x {ad.n_vars} genes", flush=True)

# lognorm of the DECONTAMINATED counts
ad.X = ad.layers["decontX"].copy()
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
ad.layers["lognorm"] = ad.X.copy()
sc.pp.highly_variable_genes(ad, n_top_genes=a.n_hvg, batch_key="gem_reaction",
                            flavor="seurat_v3", layer="counts")
adh = ad[:, ad.var["highly_variable"]].copy()

# ComBat on the (decontaminated) log-normalised HVG matrix -> primary integration
sc.pp.combat(adh, key="gem_reaction")
sc.pp.scale(adh, max_value=10)
sc.tl.pca(adh, n_comps=a.n_pcs, svd_solver="arpack")
ad.obsm["X_combat"] = adh.obsm["X_pca"].copy()

# neighbours + Leiden tuned to ~target_k clusters
sc.pp.neighbors(ad, use_rep="X_combat", n_neighbors=30, random_state=a.seed)
best = None
for res in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]:
    sc.tl.leiden(ad, resolution=res, random_state=a.seed, key_added="new_cluster")
    k = ad.obs["new_cluster"].nunique()
    if best is None or abs(k - a.target_k) < abs(best[1] - a.target_k):
        best = (res, k)
    if k >= a.target_k:
        break
sc.tl.leiden(ad, resolution=best[0], random_state=a.seed, key_added="new_cluster")
sc.tl.umap(ad, random_state=a.seed)
print(f"clustering: res={best[0]} -> {ad.obs['new_cluster'].nunique()} clusters", flush=True)

# map new clusters to original C0-C5 by maximum overlap
ct = pd.crosstab(ad.obs["new_cluster"], ad.obs["old_eos_cluster"])
mapping = ct.idxmax(axis=1).to_dict()          # new leiden id -> old C label
ad.obs["new_C"] = ad.obs["new_cluster"].map(mapping).astype(str)
ct.to_csv(os.path.join(a.outdir, "stage2_new_vs_oldC0C5_crosstab.csv"))

# verification gate: ARI + how many old clusters are recovered as a majority
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
ari = adjusted_rand_score(ad.obs["old_eos_cluster"], ad.obs["new_cluster"])
nmi = normalized_mutual_info_score(ad.obs["old_eos_cluster"], ad.obs["new_cluster"])
recovered = sorted(set(mapping.values()))
summary = pd.DataFrame({
    "metric": ["n_cells", "n_new_clusters", "resolution", "ARI_vs_oldC0C5",
               "NMI_vs_oldC0C5", "old_clusters_recovered"],
    "value": [ad.n_obs, ad.obs["new_cluster"].nunique(), best[0],
              round(ari, 3), round(nmi, 3), ",".join(recovered)]})
summary.to_csv(os.path.join(a.outdir, "stage2_verification.csv"), index=False)
print("\n=== VERIFICATION (new clustering vs original C0-C5) ===")
print(summary.to_string(index=False), flush=True)
print("\ncrosstab (new rows x old cols):\n", ct.to_string(), flush=True)
print(f"\nmapping new->old: {mapping}", flush=True)

# composition shift on the NEW clustering (sex/genotype/infection) — quick read
ad.obs["male"] = (ad.obs["sex"].astype(str) == "Male")
ad.obs["inf"] = ad.obs["condition"].astype(str).str.contains("Inf", case=False)
comp = ad.obs.groupby("new_C").agg(n=("male", "size"), pct_male=("male", "mean"),
                                   pct_inf=("inf", "mean")).round(3)
comp.to_csv(os.path.join(a.outdir, "stage2_new_cluster_composition.csv"))
print("\nnew-cluster composition:\n", comp.to_string(), flush=True)

# UMAP coords for figures
um = pd.DataFrame(ad.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"])
for c in ["new_cluster", "new_C", "old_eos_cluster", "condition", "sex", "gem_reaction"]:
    um[c] = ad.obs[c].astype(str).values
um.to_csv(os.path.join(a.outdir, "stage2_umap_cells.csv"), index=False)

ad.write_h5ad(os.path.join(a.outdir, "eos_foundation_combat.h5ad"))
print("\nDone 62b. Wrote eos_foundation_combat.h5ad", flush=True)
