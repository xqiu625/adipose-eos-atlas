#!/usr/bin/env python
# ==============================================================================
# 54b_integrated_conclusions_matched.py   (Reviewer #1, Concern #1 — fair test)
# Like 54, but tune the integrated-clustering resolution so it yields ~6 clusters
# (matched granularity to the un-integrated C0-C5), then test whether the
# sex / RELMa / infection composition shifts re-appear. This removes the
# resolution-mismatch confound (res 0.4 on the integrated graph gave only 3).
# Run in sc_scvi env.
# ==============================================================================
import argparse, os, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--target", type=int, default=6)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

ad = sc.read_h5ad(a.h5ad)
ob = ad.obs
ob["cluster"] = ob["cluster"].astype(str)
ob["sex"] = (ob["sex"].astype(str) if "sex" in ob else ob["condition"].astype(str).str.extract(r"(Male|Female)")[0]).values
ob["genotype"] = (ob["genotype"].astype(str) if "genotype" in ob else ob["condition"].astype(str).str.extract(r"(WT|KO)")[0]).values
ob["infected"] = ob["condition"].astype(str).str.contains("Inf", case=False).values
base_male = (ob["sex"] == "Male").mean(); base_inf = ob["infected"].mean()
fem = ob["sex"] == "Female"; base_koF = (ob.loc[fem, "genotype"] == "KO").mean()

def tune(emb):
    sc.pp.neighbors(ad, use_rep=emb, n_neighbors=30, random_state=a.seed)
    best = None
    for res in [0.4,0.6,0.8,1.0,1.2,1.5,2.0]:
        sc.tl.leiden(ad, resolution=res, random_state=a.seed, key_added="intg")
        nk = ad.obs["intg"].nunique()
        if best is None or abs(nk-a.target) < abs(best[1]-a.target):
            best = (res, nk)
        if nk >= a.target:
            break
    sc.tl.leiden(ad, resolution=best[0], random_state=a.seed, key_added="intg")
    return best

for emb in [e for e in ["X_pca_harmony","X_combat"] if e in ad.obsm]:
    res, nk = tune(emb)
    ct = pd.crosstab(ad.obs["cluster"], ad.obs["intg"])
    ct.to_csv(os.path.join(a.outdir, f"matched_crosstab_{emb}.csv"))
    rows = []
    for ic in sorted(ad.obs["intg"].unique(), key=int):
        sub = ob[ad.obs["intg"] == ic]; nF = (sub["sex"]=="Female").sum()
        rows.append({"embedding":emb,"res":res,"n_clusters":nk,"intg":ic,"n":len(sub),
            "dominant_origin":ct[ic].idxmax(),"dominant_frac":round(ct[ic].max()/ct[ic].sum(),2),
            "male_enrich":round((sub["sex"]=="Male").mean()-base_male,2),
            "inf_enrich":round(sub["infected"].mean()-base_inf,2),
            "KO_F_enrich":round((sub.loc[sub["sex"]=="Female","genotype"]=="KO").mean()-base_koF,2) if nF else np.nan})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.outdir, f"matched_composition_{emb}.csv"), index=False)
    print(f"\n=== [{emb}] tuned res={res} -> {nk} integrated clusters ===")
    print(df.to_string(index=False), flush=True)
    # which original clusters got a majority-dominated integrated cluster?
    recovered = df[df.dominant_frac >= 0.5]["dominant_origin"].unique().tolist()
    c0 = df[df.dominant_origin=="C0"]
    print(f"[{emb}] original clusters recovered as a majority integrated cluster: {sorted(recovered)}")
    if len(c0):
        print(f"[{emb}] C0-like cluster FOUND: male_enrich={c0.male_enrich.max():+.2f}, KO_F_enrich={c0.KO_F_enrich.max():+.2f}")
    else:
        print(f"[{emb}] NO C0-dominated integrated cluster even at res {res} (C0/MHC-II does NOT re-separate)")
print("\nDone 54b", flush=True)
