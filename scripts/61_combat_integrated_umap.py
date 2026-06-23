#!/usr/bin/env python
# ==============================================================================
# 61_combat_integrated_umap.py   (Reviewer #1 — actually perform the integration)
# Produce the ComBat-integrated eosinophil UMAP for the revised figures: neighbors
# + UMAP on the ComBat embedding already computed in the benchmark, colored by the
# original C0-C5 labels, condition, and sample. Exports coordinates for figure-
# making. Demonstrates the integration was carried out (not merely benchmarked).
# Run in sc_scvi env.
# ==============================================================================
import argparse, os, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

ad = sc.read_h5ad(a.h5ad)
assert "X_combat" in ad.obsm, "X_combat embedding not found in h5ad"
ad.obs["cluster"] = ad.obs["cluster"].astype(str)

# UMAP on the ComBat embedding (the integration the reviewer asked us to perform)
sc.pp.neighbors(ad, use_rep="X_combat", n_neighbors=30, random_state=a.seed)
sc.tl.umap(ad, random_state=a.seed)
um = pd.DataFrame(ad.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"])
for col in ["cluster", "condition", "sample", "sex", "genotype"]:
    if col in ad.obs:
        um[col] = ad.obs[col].astype(str).values
um.to_csv(os.path.join(a.outdir, "combat_umap_cells.csv"), index=False)

# per-cluster purity already known; also export cluster centroids for labels
cent = um.groupby("cluster")[["UMAP_1", "UMAP_2"]].median()
cent.to_csv(os.path.join(a.outdir, "combat_umap_centroids.csv"))
print("cells:", um.shape[0], "clusters:", sorted(um["cluster"].unique()), flush=True)
print("exported combat_umap_cells.csv", flush=True)
