#!/usr/bin/env python
# ==============================================================================
# 49_c2_integration_dispersion.py   (Reviewer #1, Concern #1)
# The reviewer's explicit question: when the data are integrated, do C2 cells
# stay together (a robust state) or disperse among other clusters (a male-
# specific transcriptional response)? Re-cluster on the integrated embeddings
# already computed in the benchmark and cross-tabulate against the original
# C0-C5 labels. Run in the sc_scvi env (scanpy).
#
#   python 49_c2_integration_dispersion.py --h5ad <integrated_eos_benchmark.h5ad> --outdir <dir>
# ==============================================================================
import argparse, os, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--label_key", default="cluster")
ap.add_argument("--res", type=float, default=0.4)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

adata = sc.read_h5ad(a.h5ad)
adata.obs[a.label_key] = adata.obs[a.label_key].astype(str)
# integration embeddings present, in order of bio-conservation preference
embs = [e for e in ["X_combat", "X_pca_harmony", "X_scVI"] if e in adata.obsm]
print("Embeddings to test:", embs, flush=True)

summary = []
for emb in embs:
    sc.pp.neighbors(adata, use_rep=emb, n_neighbors=30, random_state=a.seed)
    sc.tl.leiden(adata, resolution=a.res, random_state=a.seed, key_added="intg")
    ct = pd.crosstab(adata.obs[a.label_key], adata.obs["intg"])
    ct.to_csv(os.path.join(a.outdir, f"crosstab_origC0C5_vs_{emb}.csv"))
    # per original cluster: how concentrated it stays (purity = max share in one integrated cluster)
    for oc in sorted(adata.obs[a.label_key].unique()):
        row = ct.loc[oc]
        purity = row.max() / row.sum()
        # how many integrated clusters hold >=10% of this original cluster
        spread = int((row / row.sum() >= 0.10).sum())
        summary.append({"embedding": emb, "orig_cluster": oc, "n_cells": int(row.sum()),
                        "stays_together_purity": round(float(purity), 3),
                        "n_integrated_clusters_>=10pct": spread,
                        "top_integrated_cluster": row.idxmax()})
    print(f"\n[{emb}] crosstab (orig C0-C5 rows x integrated cols):\n{ct}\n", flush=True)

sm = pd.DataFrame(summary)
sm.to_csv(os.path.join(a.outdir, "integration_dispersion_summary.csv"), index=False)
print("=== dispersion summary (low purity / high spread = cluster dissolves on integration) ===")
print(sm.to_string(index=False), flush=True)
c2 = sm[sm.orig_cluster == "C2"]
print("\n=== C2 specifically ===")
print(c2.to_string(index=False), flush=True)
print("\nDone.", flush=True)
