#!/usr/bin/env python
# ==============================================================================
# 67_combat010_recluster.py   (execute ComBat on decont + >0.10 data)
# ComBat integration as the PRIMARY representation, on the decontaminated,
# contamination->0.10-filtered eosinophils. Re-cluster and compare to the
# published C0-C5 so the team can decide on data (un-integrated vs ComBat).
#   - filter eos_decont_trimmed.h5ad to decontx_contam <= 0.10
#   - lognorm decont -> batch HVG -> ComBat(gem_reaction) -> PCA -> Leiden(~6) -> UMAP
#   - map to old C0-C5; ARI/crosstab; per-cluster sex/infection composition
#   - export new_C labels (barcode -> cluster) for the R downstream (68)
# Run in sc_scvi env.
# ==============================================================================
import argparse, os, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True); ap.add_argument("--outdir", required=True)
ap.add_argument("--target_k", type=int, default=6); ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)

ad = sc.read_h5ad(a.h5ad)
ad = ad[ad.obs["decontx_contam"].astype(float) <= 0.10].copy()      # >0.10 filter
ad.obs["gem_reaction"] = ad.obs["gem_reaction"].astype("category")
ad.layers["counts"] = ad.X.copy()
print(f"cells after >0.10 filter: {ad.n_obs}", flush=True)

ad.X = ad.layers["decontX"].copy()
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad); ad.layers["lognorm"] = ad.X.copy()
sc.pp.highly_variable_genes(ad, n_top_genes=2000, batch_key="gem_reaction",
                            flavor="seurat_v3", layer="counts")
adh = ad[:, ad.var["highly_variable"]].copy()
sc.pp.combat(adh, key="gem_reaction")
sc.pp.scale(adh, max_value=10); sc.tl.pca(adh, n_comps=30, svd_solver="arpack")
ad.obsm["X_combat"] = adh.obsm["X_pca"].copy()

sc.pp.neighbors(ad, use_rep="X_combat", n_neighbors=30, random_state=a.seed)
best = None
for res in [0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2]:
    sc.tl.leiden(ad, resolution=res, random_state=a.seed, key_added="new_cluster")
    k = ad.obs["new_cluster"].nunique()
    if best is None or abs(k-a.target_k) < abs(best[1]-a.target_k): best=(res,k)
    if k >= a.target_k: break
sc.tl.leiden(ad, resolution=best[0], random_state=a.seed, key_added="new_cluster")
ct = pd.crosstab(ad.obs["new_cluster"], ad.obs["old_eos_cluster"])
mapping = ct.idxmax(axis=1).to_dict()
ad.obs["new_C"] = ad.obs["new_cluster"].map(mapping).astype(str)
ct.to_csv(os.path.join(a.outdir, "combat010_new_vs_oldC0C5.csv"))

from sklearn.metrics import adjusted_rand_score
ari = adjusted_rand_score(ad.obs["old_eos_cluster"], ad.obs["new_cluster"])
recovered = sorted(set(mapping.values()))
print(f"\nres={best[0]} -> {ad.obs['new_cluster'].nunique()} clusters; ARI vs C0-C5={ari:.3f}; recovered={recovered}", flush=True)
print("crosstab (new rows x old cols):\n", ct.to_string(), flush=True)

ad.obs["male"] = (ad.obs["sex"].astype(str)=="Male")
ad.obs["inf"]  = ad.obs["condition"].astype(str).str.contains("Inf", case=False)
comp = ad.obs.groupby("new_C").agg(n=("male","size"), pct_male=("male","mean"), pct_inf=("inf","mean")).round(3)
comp.to_csv(os.path.join(a.outdir, "combat010_composition.csv"))
print("\nComBat new-cluster composition:\n", comp.to_string(), flush=True)
pd.DataFrame({"barcode": ad.obs_names, "new_cluster": ad.obs["new_cluster"].values,
              "new_C": ad.obs["new_C"].values}).to_csv(
    os.path.join(a.outdir, "combat010_labels.csv"), index=False)
pd.DataFrame({"metric":["n_cells","n_clusters","resolution","ARI_vs_C0C5","recovered"],
              "value":[ad.n_obs, ad.obs["new_cluster"].nunique(), best[0], round(ari,3), ",".join(recovered)]}
            ).to_csv(os.path.join(a.outdir,"combat010_verification.csv"), index=False)
ad.write_h5ad(os.path.join(a.outdir, "eos_combat010.h5ad"))
print("\nDone 67.", flush=True)
