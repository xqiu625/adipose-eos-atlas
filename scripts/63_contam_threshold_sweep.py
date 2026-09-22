#!/usr/bin/env python
# ==============================================================================
# 63_contam_threshold_sweep.py   (raise the contamination bar — quantify)
# From the decontaminated eosinophil object, quantify how many cells a stricter
# contamination bar would remove, per cluster, using two complementary metrics:
#   (a) per-cell DecontX contamination fraction (obs[decontx_contam])
#   (b) a macrophage-associated module score (abundant ambient mac genes:
#       C1qa/C1qb/Apoe/Lyz2/Folr2/F13a1) on the decontaminated lognorm matrix
# Outputs per-cluster distributions + a removal sweep so we can pick a defensible
# stricter threshold before re-trimming.
# Run in sc_scvi env.
# ==============================================================================
import argparse, os, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

ad = sc.read_h5ad(a.h5ad)
cl = ad.obs["old_eos_cluster"].astype(str)
contam = ad.obs["decontx_contam"].astype(float).values

# mac-associated module score on decontaminated lognorm
ad.X = ad.layers["decontX"].copy()
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
mac_genes = [g for g in ["C1qa","C1qb","Apoe","Lyz2","Folr2","F13a1"] if g in ad.var_names]
sc.tl.score_genes(ad, mac_genes, score_name="mac_score")
mac = ad.obs["mac_score"].values

clusters = sorted(cl.unique())
# distribution of DecontX contamination per cluster
dist = []
for c in clusters:
    v = contam[cl.values == c]
    dist.append({"cluster": c, "n": int(v.size), "median": round(np.median(v),3),
                 "p90": round(np.percentile(v,90),3), "p95": round(np.percentile(v,95),3),
                 "p99": round(np.percentile(v,99),3), "max": round(v.max(),3)})
pd.DataFrame(dist).to_csv(os.path.join(a.outdir,"contam_distribution_by_cluster.csv"), index=False)
print("=== DecontX contamination distribution by cluster ===")
print(pd.DataFrame(dist).to_string(index=False), flush=True)

# removal sweep: DecontX-contamination thresholds
rows = []
for thr in [0.30, 0.25, 0.20, 0.15, 0.10]:
    rem = contam > thr
    row = {"metric": f"DecontX_contam>{thr}", "total_removed": int(rem.sum()),
           "pct_total": round(100*rem.mean(),1)}
    for c in clusters:
        row[c] = int((rem & (cl.values==c)).sum())
    rows.append(row)
# removal sweep: mac-score top percentiles (global)
for pct in [99, 97.5, 95, 90]:
    thr = np.percentile(mac, pct); rem = mac > thr
    row = {"metric": f"mac_score>p{pct}", "total_removed": int(rem.sum()),
           "pct_total": round(100*rem.mean(),1)}
    for c in clusters:
        row[c] = int((rem & (cl.values==c)).sum())
    rows.append(row)
sweep = pd.DataFrame(rows)
sweep.to_csv(os.path.join(a.outdir,"contam_removal_sweep.csv"), index=False)
print("\n=== removal sweep (cells removed at each bar; per-cluster columns) ===")
print(sweep.to_string(index=False), flush=True)

# combined bar candidate: DecontX>0.2 OR mac_score>p95 (a 'raised bar')
for label, mask in [("contam>0.2", contam>0.2),
                    ("contam>0.2 OR mac>p95", (contam>0.2) | (mac>np.percentile(mac,95))),
                    ("contam>0.15 OR mac>p90", (contam>0.15) | (mac>np.percentile(mac,90)))]:
    by = pd.Series(mask).groupby(cl.values).sum()
    print(f"\ncandidate [{label}] total={int(mask.sum())} ({100*mask.mean():.1f}%): "
          + ", ".join(f"{c}={int(by.get(c,0))}" for c in clusters), flush=True)
print("\nDone 63.", flush=True)
