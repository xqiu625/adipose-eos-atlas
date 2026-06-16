#!/usr/bin/env python
# ==============================================================================
# 54_integrated_conclusions.py   (Reviewer #1, Concern #1 — robustness of CONCLUSIONS)
# Re-cluster on the integrated embedding, map integrated clusters back to the
# manuscript C0-C5, and test whether the key compositional conclusions still hold:
#   - a MHC-II (C0-like) cluster enriched in MALES and expanded in KO females
#   - a migratory (C1-like) cluster enriched in FEMALES
#   - infection-enriched cluster(s) (C4/C5-like)
# If the same sex/RELMa/infection shifts are recovered after integration, the
# biological conclusions do not depend on the un-integrated choice.
# Run in sc_scvi env (scanpy).
# ==============================================================================
import argparse, os, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--res", type=float, default=0.4)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

ad = sc.read_h5ad(a.h5ad)
ob = ad.obs
ob["cluster"] = ob["cluster"].astype(str)
# derive covariates
sx = ob["sex"].astype(str) if "sex" in ob else ob["condition"].astype(str).str.extract(r"(Male|Female)")[0]
gt = ob["genotype"].astype(str) if "genotype" in ob else ob["condition"].astype(str).str.extract(r"(WT|KO)")[0]
inf = ob["condition"].astype(str).str.contains("Inf", case=False)
ob["sex"], ob["genotype"], ob["infected"] = sx.values, gt.values, inf.values

base_male = (ob["sex"] == "Male").mean()
base_inf  = ob["infected"].mean()
fem = ob["sex"] == "Female"
base_ko_inF = (ob.loc[fem, "genotype"] == "KO").mean()
print(f"baselines: %male={base_male:.2f} %infected={base_inf:.2f} %KO(among F)={base_ko_inF:.2f}", flush=True)

for emb in [e for e in ["X_pca_harmony", "X_combat"] if e in ad.obsm]:
    sc.pp.neighbors(ad, use_rep=emb, n_neighbors=30, random_state=a.seed)
    sc.tl.leiden(ad, resolution=a.res, random_state=a.seed, key_added="intg")
    ct = pd.crosstab(ad.obs["cluster"], ad.obs["intg"])
    ct.to_csv(os.path.join(a.outdir, f"concl_crosstab_{emb}.csv"))

    rows = []
    for ic in sorted(ad.obs["intg"].unique(), key=int):
        m = ad.obs["intg"] == ic
        sub = ad.obs[m]
        dom = ct[ic].idxmax(); dom_frac = ct[ic].max() / ct[ic].sum()
        nF = (sub["sex"] == "Female").sum()
        rows.append({
            "embedding": emb, "intg_cluster": ic, "n": int(m.sum()),
            "dominant_origin": dom, "dominant_frac": round(dom_frac, 2),
            "pct_male": round(100*(sub["sex"] == "Male").mean(), 1),
            "male_enrich": round((sub["sex"] == "Male").mean() - base_male, 2),
            "pct_infected": round(100*sub["infected"].mean(), 1),
            "inf_enrich": round(sub["infected"].mean() - base_inf, 2),
            "pct_KO_amongF": round(100*(sub.loc[sub["sex"]=="Female","genotype"]=="KO").mean(), 1) if nF else np.nan,
            "KO_F_enrich": round((sub.loc[sub["sex"]=="Female","genotype"]=="KO").mean() - base_ko_inF, 2) if nF else np.nan,
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.outdir, f"concl_composition_{emb}.csv"), index=False)
    print(f"\n=== [{emb}] integrated clusters: origin + condition enrichment ===")
    print(df.to_string(index=False), flush=True)

    # conclusion-recovery check, tied to manuscript claims
    c0c = df[df.dominant_origin == "C0"]
    c1c = df[df.dominant_origin == "C1"]
    infc = df[df.dominant_origin.isin(["C4", "C5"])]
    print(f"\n[{emb}] CONCLUSION RECOVERY:")
    if len(c0c): print(f"  C0-like (MHC-II) cluster male_enrich={c0c.male_enrich.max():+.2f}, KO_F_enrich={c0c.KO_F_enrich.max():+.2f}  (manuscript: male+ / KO-female expanded)")
    if len(c1c): print(f"  C1-like (migratory) cluster male_enrich={c1c.male_enrich.min():+.2f}  (manuscript: female+, i.e. male_enrich negative)")
    if len(infc): print(f"  C4/C5-like cluster inf_enrich={infc.inf_enrich.max():+.2f}  (manuscript: infection+)")

print("\nDone 54_integrated_conclusions.py", flush=True)
