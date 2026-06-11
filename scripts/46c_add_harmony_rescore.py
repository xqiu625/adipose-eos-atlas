#!/usr/bin/env python
# ==============================================================================
# 46c_add_harmony_rescore.py
# Inject the Harmony embedding into an existing benchmark result and re-score all
# embeddings together — used when a 46 run completed but Harmony failed (the
# scanpy 1.11 / harmonypy wrapper mis-shapes X_pca_harmony). Avoids re-training
# scVI: reads the saved integrated_eos_benchmark.h5ad, which already holds every
# other embedding plus X_pca (Harmony's input).
#
# Computes Harmony directly from X_pca, appends its Leiden-vs-C0C5 recovery, then
# re-runs scib-metrics over the full embedding set and rewrites the metrics CSVs
# and the scib figure.
#
# Usage:
#   python 46c_add_harmony_rescore.py \
#       --h5ad   $B/integrated_eos_benchmark.h5ad \
#       --outdir $B --batch_key gem_reaction --label_key cluster
# ==============================================================================

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--batch_key", default="gem_reaction")
ap.add_argument("--label_key", default="cluster")
ap.add_argument("--leiden_res", type=float, default=0.4)
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()

print(f"Loading {args.h5ad}", flush=True)
adata = sc.read_h5ad(args.h5ad)
adata.obs[args.batch_key] = adata.obs[args.batch_key].astype("category")
adata.obs[args.label_key] = adata.obs[args.label_key].astype("category")
if "X_pca" not in adata.obsm:
    raise KeyError("X_pca not in obsm; cannot compute Harmony")

# ── Harmony, direct (robust to the scanpy wrapper shape bug) ──────────────────
print("Running Harmony directly via harmonypy", flush=True)
import harmonypy
ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [args.batch_key])
adata.obsm["X_pca_harmony"] = np.asarray(ho.Z_corr).T.copy()
assert adata.obsm["X_pca_harmony"].shape[0] == adata.n_obs, "harmony shape wrong"
print(f"  X_pca_harmony shape {adata.obsm['X_pca_harmony'].shape}", flush=True)

# ── Harmony cluster recovery, appended to the existing CSV ─────────────────────
tmp = adata.copy()
sc.pp.neighbors(tmp, use_rep="X_pca_harmony", n_neighbors=30, random_state=args.seed)
sc.tl.leiden(tmp, resolution=args.leiden_res, random_state=args.seed, key_added="_l")
truth = tmp.obs[args.label_key].astype(str).values
pred = tmp.obs["_l"].astype(str).values
h_ari = adjusted_rand_score(truth, pred)
h_nmi = normalized_mutual_info_score(truth, pred)
h_ncl = tmp.obs["_l"].nunique()
print(f"  harmony recovery: ARI={h_ari:.3f} NMI={h_nmi:.3f} ({h_ncl} clusters)", flush=True)

rec_path = os.path.join(args.outdir, "integration_cluster_recovery.csv")
rec = pd.read_csv(rec_path)
rec = rec[rec["method"] != "harmony"]
rec = pd.concat([rec, pd.DataFrame([{
    "method": "harmony", "leiden_ARI_vs_C0C5": h_ari,
    "leiden_NMI_vs_C0C5": h_nmi, "n_leiden_clusters": h_ncl}])], ignore_index=True)
order = ["unintegrated", "combat", "harmony", "scanorama", "scvi", "scanvi"]
rec["__o"] = rec["method"].map({m: i for i, m in enumerate(order)}).fillna(99)
rec = rec.sort_values("__o").drop(columns="__o")
rec.to_csv(rec_path, index=False)

# ── Re-score every embedding together (consistent min-max scaling) ────────────
from scib_metrics.benchmark import Benchmarker
emb = [k for k in ["X_pca", "X_combat", "X_pca_harmony", "X_scanorama",
                   "X_scVI", "X_scANVI"] if k in adata.obsm]
print(f"Scoring {len(emb)} embeddings with scib-metrics: {emb}", flush=True)
bm = Benchmarker(adata, batch_key=args.batch_key, label_key=args.label_key,
                 embedding_obsm_keys=emb, n_jobs=-1)
bm.benchmark()
bm.get_results(min_max_scale=False).to_csv(
    os.path.join(args.outdir, "integration_benchmark_metrics.csv"))
bm.get_results(min_max_scale=True).to_csv(
    os.path.join(args.outdir, "integration_benchmark_metrics_scaled.csv"))
try:
    bm.plot_results_table(min_max_scale=False, save_dir=args.outdir)
    src = os.path.join(args.outdir, "scib_results.svg")
    if os.path.exists(src):
        os.replace(src, os.path.join(args.outdir, "integration_benchmark_scib.svg"))
except Exception as e:
    print(f"  ! scib plot failed (CSVs written): {e}", flush=True)

adata.write_h5ad(args.h5ad)
print("done", flush=True)
