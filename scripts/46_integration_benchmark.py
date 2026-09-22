#!/usr/bin/env python
# ==============================================================================
# 46_integration_benchmark.py
# Benchmark single-cell integration / batch-correction tools on the adipose
# eosinophil atlas and let the data decide whether integration helps.
#
# Methods span the categories in the AI4Bio batch-correction module:
#   - Unintegrated   : PCA on log-normalized HVGs (baseline, current primary)
#   - ComBat         : classical empirical-Bayes correction, then PCA
#   - Harmony        : iterative PCA-space correction (harmonypy)
#   - Scanorama      : mutual-nearest-neighbor panorama stitching
#   - scVI           : negative-binomial VAE, batch as covariate
#   - scANVI         : scVI + label supervision (semi-supervised)
#   - BBKNN          : batch-balanced kNN graph (graph-only; scored separately)
#
# Each embedding method writes an embedding to adata.obsm; scib-metrics scores
# every embedding on two axes:
#   - Batch correction   (iLISI, kBET, graph connectivity, PCR, silhouette-batch)
#   - Bio conservation   (ARI, NMI, cLISI, isolated labels, silhouette-label)
# The "Total" column is the published scIB weighted average (0.4*batch + 0.6*bio).
#
# For each method we ALSO run Leiden on its embedding and report ARI / NMI versus
# the manuscript C0-C5 labels — the cluster-recovery view the reviewers asked for
# (a method that scrambles C0-C5 is over-correcting regardless of its batch score).
#
# batch_key defaults to the technical multiplex batch (gem_reaction); pass
# --batch_key sample to score per-library mixing instead. label_key is the
# bio-conservation ground truth (manuscript clusters).
#
# Usage (see 46b_run_integration_benchmark.sh):
#   python 46_integration_benchmark.py \
#       --h5ad   $BASE/integration/benchmark/eos_for_integration.h5ad \
#       --outdir $BASE/integration/benchmark \
#       --batch_key gem_reaction --label_key cluster
#
# Outputs (in --outdir):
#   integration_benchmark_metrics.csv        full scib-metrics table (raw scores)
#   integration_benchmark_metrics_scaled.csv min-max scaled table (for the plot)
#   integration_benchmark_scib.svg/.png      scib-metrics summary table figure
#   integration_cluster_recovery.csv         Leiden-vs-C0C5 ARI/NMI per method
#   integration_umap_grid.png                UMAP per method (batch + label)
#   integrated_eos_benchmark.h5ad            adata with every method's embedding
#   benchmark_run.log                        method-by-method run log
# ==============================================================================

import argparse
import gc
import importlib.util
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Single-cell integration benchmark")
    p.add_argument("--h5ad", required=True, help="Input h5ad (raw counts in X)")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--batch_key", default="gem_reaction",
                   help="obs column used as the batch to correct")
    p.add_argument("--label_key", default="cluster",
                   help="obs column used as the bio-conservation ground truth")
    p.add_argument("--n_hvg", type=int, default=2000)
    p.add_argument("--n_pcs", type=int, default=30)
    p.add_argument("--n_latent", type=int, default=30, help="scVI/scANVI latent dim")
    p.add_argument("--scvi_epochs", type=int, default=400)
    p.add_argument("--leiden_res", type=float, default=0.4,
                   help="Leiden resolution for cluster-recovery (matches manuscript)")
    p.add_argument("--methods", default="all",
                   help="comma list, or 'all'. options: unintegrated,combat,"
                        "harmony,scanorama,scvi,scanvi,bbknn")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────
class Tee:
    """Mirror stdout to a log file."""
    def __init__(self, path):
        self.term = sys.stdout
        self.log = open(path, "w")

    def write(self, m):
        self.term.write(m); self.log.write(m); self.log.flush()

    def flush(self):
        self.term.flush(); self.log.flush()


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def leiden_recovery(adata, use_rep, label_key, res, seed):
    """Cluster on an embedding and score agreement with the manuscript labels."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    tmp = adata.copy()
    sc.pp.neighbors(tmp, use_rep=use_rep, n_neighbors=30, random_state=seed)
    sc.tl.leiden(tmp, resolution=res, random_state=seed, key_added="_bench_leiden")
    truth = tmp.obs[label_key].astype(str).values
    pred = tmp.obs["_bench_leiden"].astype(str).values
    ari = adjusted_rand_score(truth, pred)
    nmi = normalized_mutual_info_score(truth, pred)
    n_clusters = tmp.obs["_bench_leiden"].nunique()
    del tmp; gc.collect()
    return ari, nmi, n_clusters


# ──────────────────────────────────────────────────────────────────────────────
# Integration methods — each returns the obsm key it wrote, or None on failure
# ──────────────────────────────────────────────────────────────────────────────
def run_unintegrated(adata, args):
    # X_pca already computed in preprocessing
    return "X_pca"


def run_combat(adata, args):
    a = adata.copy()
    sc.pp.combat(a, key=args.batch_key)
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=args.n_pcs, svd_solver="arpack")
    adata.obsm["X_combat"] = a.obsm["X_pca"]
    del a; gc.collect()
    return "X_combat"


def run_harmony(adata, args):
    # Call harmonypy directly rather than sc.external.pp.harmony_integrate: the
    # scanpy wrapper mis-shapes the output under scanpy 1.11 / harmonypy >=0.2
    # (writes a (n_pcs,) vector instead of (n_cells, n_pcs)). Z_corr is
    # (n_pcs, n_cells); transpose to cells-by-PCs.
    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [args.batch_key])
    adata.obsm["X_pca_harmony"] = np.asarray(ho.Z_corr).T.copy()
    return "X_pca_harmony"


def run_scanorama(adata, args):
    # scanorama integrates per-batch; cells must be contiguous by batch
    order = np.argsort(adata.obs[args.batch_key].astype(str).values, kind="stable")
    a = adata[order].copy()
    # approx=False -> exact NN (sklearn) instead of annoy. annoy's SIMD kernels
    # SIGILL ("Illegal instruction") on some HPCC compute nodes; exact is fine
    # at this cell count and cannot crash the interpreter.
    sc.external.pp.scanorama_integrate(a, key=args.batch_key, basis="X_pca",
                                       adjusted_basis="X_scanorama", approx=False)
    # scatter back to original order
    emb = np.zeros((adata.n_obs, a.obsm["X_scanorama"].shape[1]), dtype=np.float32)
    emb[order] = a.obsm["X_scanorama"]
    adata.obsm["X_scanorama"] = emb
    del a; gc.collect()
    return "X_scanorama"


def run_scvi(adata, args):
    import scvi
    scvi.settings.seed = args.seed
    a = adata.copy()
    a.X = a.layers["counts"].copy()
    scvi.model.SCVI.setup_anndata(a, batch_key=args.batch_key)
    model = scvi.model.SCVI(a, n_layers=2, n_latent=args.n_latent,
                            gene_likelihood="nb")
    model.train(max_epochs=args.scvi_epochs, early_stopping=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.uns["_scvi_model"] = model  # reused by scANVI
    del a; gc.collect()
    return "X_scVI"


def run_scanvi(adata, args):
    import scvi
    model = adata.uns.get("_scvi_model")
    if model is None:
        run_scvi(adata, args)
        model = adata.uns["_scvi_model"]
    scanvi = scvi.model.SCANVI.from_scvi_model(
        model, adata=model.adata, labels_key=args.label_key, unlabeled_category="Unknown")
    scanvi.train(max_epochs=20, n_samples_per_label=100)
    adata.obsm["X_scANVI"] = scanvi.get_latent_representation()
    return "X_scANVI"


def run_bbknn(adata, args):
    # graph-only method: produce a 2-D embedding so it can be scored alongside
    # the embedding methods (PCA-UMAP on the batch-balanced graph).
    import bbknn  # noqa: F401
    a = adata.copy()
    # use_annoy=False + approx=False -> exact NN, avoids the annoy SIGILL (see scanorama)
    sc.external.pp.bbknn(a, batch_key=args.batch_key, use_rep="X_pca",
                         approx=False, use_annoy=False)
    sc.tl.umap(a, random_state=args.seed)
    adata.obsm["X_bbknn"] = a.obsm["X_umap"]
    del a; gc.collect()
    return "X_bbknn"


METHODS = {
    "unintegrated": run_unintegrated,
    "combat": run_combat,
    "harmony": run_harmony,
    "scanorama": run_scanorama,
    "scvi": run_scvi,
    "scanvi": run_scanvi,
    "bbknn": run_bbknn,
}

# Top-level package each method needs, for a preflight availability check.
# environment.yml only covers the core (unintegrated/combat/harmony) path;
# scanorama, scvi-tools (scvi/scanvi), and bbknn need separate environments
# (scvi-tools in particular conflicts with the core env's anndata pin -- see
# git history). None of that is a reason to let a partial run report success:
# this preflight check, plus the exit-code check at the end of main(), make a
# silently incomplete benchmark impossible to mistake for a complete one.
METHOD_REQUIRES = {
    "unintegrated": None,
    "combat": None,          # part of scanpy itself
    "harmony": "harmonypy",
    "scanorama": "scanorama",
    "scvi": "scvi",
    "scanvi": "scvi",
    "bbknn": "bbknn",
}


def check_method_availability(method_list):
    """Report which requested methods have their required package importable.
    Returns (available, unavailable) method-name lists."""
    available, unavailable = [], []
    for m in method_list:
        pkg = METHOD_REQUIRES.get(m)
        if pkg is None or importlib.util.find_spec(pkg) is not None:
            available.append(m)
        else:
            unavailable.append((m, pkg))
    return available, unavailable


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sys.stdout = Tee(outdir / "benchmark_run.log")
    sc.settings.figdir = str(outdir)

    if args.methods == "all":
        method_list = list(METHODS.keys())
    else:
        method_list = [m.strip() for m in args.methods.split(",") if m.strip()]

    # Preflight: report method availability BEFORE the expensive preprocessing
    # runs, rather than discovering a missing package deep into a long job.
    # This does not abort by itself (this script's purpose is comparing
    # whatever methods are available), but every unavailable method is
    # recorded and reflected in the exit code at the end of main() -- a
    # partial run can no longer look identical to a complete one.
    available, unavailable = check_method_availability(method_list)
    failed_methods = []
    if unavailable:
        log("Preflight: the following requested methods are UNAVAILABLE "
            "(package not importable in this environment) and will be skipped:")
        for m, pkg in unavailable:
            log(f"  ! {m}: requires '{pkg}' (not installed)")
        failed_methods.extend(m for m, _ in unavailable)
    log(f"Preflight: {len(available)}/{len(method_list)} requested methods available: "
        f"{available}")
    method_list = available

    log(f"Loading {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    log(f"  {adata.n_obs} cells x {adata.n_vars} genes")
    for k in (args.batch_key, args.label_key):
        if k not in adata.obs:
            raise KeyError(f"obs column '{k}' not found. have: {list(adata.obs.columns)}")
    adata.obs[args.batch_key] = adata.obs[args.batch_key].astype("category")
    adata.obs[args.label_key] = adata.obs[args.label_key].astype("category")
    log(f"  batch_key={args.batch_key} ({adata.obs[args.batch_key].nunique()} levels)"
        f" | label_key={args.label_key} ({adata.obs[args.label_key].nunique()} levels)")

    # ── Preprocessing shared by all methods ──────────────────────────────────
    log("Preprocessing: counts -> lognorm -> HVG -> scale -> PCA")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg,
                                batch_key=args.batch_key, flavor="seurat_v3",
                                layer="counts", subset=True)
    log(f"  kept {adata.n_vars} HVGs")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=args.n_pcs, svd_solver="arpack")
    adata.X = adata.layers["lognorm"].copy()  # restore lognorm for ComBat/scanorama

    # ── Run each method ──────────────────────────────────────────────────────
    embeddings = {}        # display name -> obsm key
    for m in method_list:
        if m not in METHODS:
            log(f"  ! unknown method '{m}', skipping")
            failed_methods.append(m)
            continue
        log(f"Running: {m}")
        t0 = time.time()
        try:
            key = METHODS[m](adata, args)
            embeddings[m] = key
            log(f"  {m} -> obsm['{key}']  ({time.time() - t0:.0f}s)")
        except Exception as e:  # one bad method must not kill the run
            log(f"  ! {m} FAILED: {type(e).__name__}: {e}")
            failed_methods.append(m)

    if "_scvi_model" in adata.uns:
        del adata.uns["_scvi_model"]

    # ── Cluster recovery (Leiden vs manuscript C0-C5) ────────────────────────
    log("Cluster recovery: Leiden on each embedding vs manuscript labels")
    rec_rows = []
    for name, key in embeddings.items():
        if adata.obsm[key].shape[1] < 3:  # bbknn 2-D umap: skip recovery
            log(f"  {name}: 2-D embedding, skipping Leiden recovery")
            continue
        ari, nmi, ncl = leiden_recovery(adata, key, args.label_key,
                                        args.leiden_res, args.seed)
        rec_rows.append({"method": name, "leiden_ARI_vs_C0C5": ari,
                         "leiden_NMI_vs_C0C5": nmi, "n_leiden_clusters": ncl})
        log(f"  {name}: ARI={ari:.3f} NMI={nmi:.3f} ({ncl} clusters)")
    if rec_rows:
        pd.DataFrame(rec_rows).to_csv(outdir / "integration_cluster_recovery.csv",
                                      index=False)

    # ── scib-metrics on the >2-D embeddings ──────────────────────────────────
    emb_keys = [k for k in embeddings.values() if adata.obsm[k].shape[1] >= 3]
    log(f"Scoring {len(emb_keys)} embeddings with scib-metrics: {emb_keys}")
    try:
        from scib_metrics.benchmark import Benchmarker
        bm = Benchmarker(adata, batch_key=args.batch_key, label_key=args.label_key,
                         embedding_obsm_keys=emb_keys, n_jobs=-1)
        bm.benchmark()
        res_raw = bm.get_results(min_max_scale=False)
        res_scaled = bm.get_results(min_max_scale=True)
        res_raw.to_csv(outdir / "integration_benchmark_metrics.csv")
        res_scaled.to_csv(outdir / "integration_benchmark_metrics_scaled.csv")
        try:
            bm.plot_results_table(min_max_scale=False,
                                  save_dir=str(outdir))
            # scib-metrics writes scib_results.svg; give it our name
            src = outdir / "scib_results.svg"
            if src.exists():
                src.rename(outdir / "integration_benchmark_scib.svg")
        except Exception as e:
            log(f"  ! scib plot failed (metrics CSVs still written): {e}")
        log("  scib-metrics done. Summary (Total):")
        if "Total" in res_raw.columns:
            # last row is a non-numeric "Metric Type" header; coerce + drop it
            tot = pd.to_numeric(res_raw["Total"], errors="coerce").dropna()
            print(tot.sort_values(ascending=False).to_string())
    except Exception as e:
        log(f"  ! scib-metrics FAILED: {type(e).__name__}: {e}")
        log("    (per-method embeddings + cluster-recovery CSV are still saved)")
        # scib-metrics produces this script's principal output (the batch-
        # correction/bio-conservation scores) -- its failure must not let the
        # run exit 0 with integration_benchmark_metrics*.csv missing.
        failed_methods.append("scib-metrics (batch-correction/bio-conservation scoring)")

    # ── UMAP grid per method ─────────────────────────────────────────────────
    log("Building per-method UMAP grid")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = len(embeddings)
        fig, axes = plt.subplots(2, n, figsize=(4 * n, 8), squeeze=False)
        for j, (name, key) in enumerate(embeddings.items()):
            if adata.obsm[key].shape[1] >= 3:
                sc.pp.neighbors(adata, use_rep=key, n_neighbors=30,
                                random_state=args.seed)
                sc.tl.umap(adata, random_state=args.seed)
                umap = adata.obsm["X_umap"]
            else:
                umap = adata.obsm[key]  # bbknn already 2-D
            for i, color in enumerate([args.batch_key, args.label_key]):
                ax = axes[i][j]
                cats = adata.obs[color].astype("category")
                for c in cats.cat.categories:
                    m = (cats == c).values
                    ax.scatter(umap[m, 0], umap[m, 1], s=2, label=str(c))
                ax.set_title(f"{name} — {color}", fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
            axes[0][j].legend(markerscale=4, fontsize=5, loc="best", frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / "integration_umap_grid.png", dpi=150)
        plt.close(fig)
    except Exception as e:
        log(f"  ! UMAP grid failed: {e}")

    # ── Persist the integrated object ────────────────────────────────────────
    log("Saving integrated_eos_benchmark.h5ad")
    adata.X = adata.layers["lognorm"]
    for lyr in ("counts", "lognorm"):
        if lyr in adata.layers and lyr != "counts":
            del adata.layers[lyr]
    adata.write_h5ad(outdir / "integrated_eos_benchmark.h5ad")

    log("Done. Methods run: " + ", ".join(embeddings.keys()))

    if failed_methods:
        log(f"INCOMPLETE: {len(failed_methods)} requested method(s) did not "
            f"produce an embedding: {failed_methods}")
        log("Exiting nonzero -- do not treat this as a completed benchmark.")
        sys.exit(1)


if __name__ == "__main__":
    main()
