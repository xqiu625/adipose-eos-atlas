#!/usr/bin/env python
# Export ComBat embedding + UMAP from eos_combat010.h5ad to CSV (avoids the
# anndata-R uns/log1p incompatibility when attaching reductions in Seurat).
import argparse, os
import numpy as np, pandas as pd, scanpy as sc
ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True); ap.add_argument("--outdir", required=True)
a = ap.parse_args()
ad = sc.read_h5ad(a.h5ad)
if "X_umap" not in ad.obsm:                      # 67 saved X_combat but not UMAP
    sc.pp.neighbors(ad, use_rep="X_combat", n_neighbors=30, random_state=0)
    sc.tl.umap(ad, random_state=0)
emb = pd.DataFrame(ad.obsm["X_combat"], index=ad.obs_names,
                   columns=[f"combat_{i+1}" for i in range(ad.obsm["X_combat"].shape[1])])
emb["UMAP_1"] = ad.obsm["X_umap"][:, 0]
emb["UMAP_2"] = ad.obsm["X_umap"][:, 1]
emb.insert(0, "barcode", ad.obs_names)
emb.to_csv(os.path.join(a.outdir, "combat010_embedding.csv"), index=False)
print("wrote combat010_embedding.csv", emb.shape, flush=True)
