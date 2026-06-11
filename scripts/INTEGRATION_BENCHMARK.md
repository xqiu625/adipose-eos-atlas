# Integration-Method Benchmark

Collaborator-requested benchmark of single-cell batch-correction / integration
tools on the adipose eosinophil atlas. Rather than commit to one method, we run
the standard panel and let `scib-metrics` quantify the batch-mixing vs.
bio-conservation trade-off for each. Method selection follows the AI4Bio
[batch-correction module](https://xqiu625.github.io/ai4bio-concepts/) categories.

## What is benchmarked

| Method | Category (module) | Output scored | Notes |
|--------|-------------------|---------------|-------|
| **Unintegrated** | baseline | PCA on log-norm HVGs | current primary representation |
| **ComBat** | classical / empirical Bayes | PCA on corrected matrix | robust for few samples/batch |
| **Harmony** | embedding correction | corrected PCs | fast CPU; the method tested in scripts 19/37 |
| **Scanorama** | MNN / panorama stitching | corrected embedding | avoids over-correcting non-shared types |
| **scVI** | deep generative (NB-VAE) | latent z | batch as covariate, raw counts |
| **scANVI** | semi-supervised VAE | latent z | scVI + C0–C5 label supervision |
| **BBKNN** | batch-balanced graph | UMAP (graph-only) | scored on graph metrics only |

Foundation models (scGPT, UCE) and cross-species tools (SATURN) from the module
are intentionally excluded: they are human-pretrained and bring their own batch
axis, which is the wrong tool for a 16-sample mouse eosinophil study. They are
noted for completeness only.

## Two evaluation axes

1. **scib-metrics** (`46_integration_benchmark.py`) — the published scIB panel.
   - *Batch correction*: iLISI, kBET, graph connectivity, PCR, silhouette-batch.
   - *Bio conservation*: ARI, NMI, cLISI, isolated-label F1, silhouette-label.
   - **Total** = 0.4·batch + 0.6·bio (scIB weighting). Higher = better overall,
     but read the two sub-scores separately — a high batch score with a collapsed
     bio score is over-correction.

2. **Cluster recovery** — Leiden (res 0.4, matching the manuscript) on each
   method's embedding, scored by ARI/NMI against the manuscript C0–C5 labels.
   This is the reviewer-facing number: a method that scrambles C0–C5 is erasing
   the biology the paper is about.

### Caveat on the bio label

There is no integration-independent ground-truth annotation within this study, so
the bio-conservation `label_key` is the un-integrated C0–C5 assignment. This
mildly favors the unintegrated baseline, so the headline comparison is the
**batch-mixing gain** each method buys versus the **bio-conservation it spends**,
not the absolute Total alone.

## Batch key

Default `batch_key = gem_reaction` — the 10x Flex multiplex reaction (GEM1–GEM4),
the only true technical batch (4 samples pooled per GEM). Because the multiplex
design partially confounds GEM reaction with biological condition, aggressive
integration is expected to remove real signal; the benchmark makes that visible
instead of assuming it. A `--batch_key sample` sensitivity run (per-library) is
provided commented-out in `46b`.

## How to run

```bash
# one-time, in the project conda env
pip install harmonypy scanorama bbknn scvi-tools scib-metrics

sbatch scripts/46b_run_integration_benchmark.sh
```

`46b` runs the export (`45`) then the benchmark (`46`). scVI/scANVI use the GPU
partition; for a CPU-only subset use `--methods unintegrated,combat,harmony,scanorama,bbknn`.

## Outputs (`$BASE/integration/benchmark/`)

| File | Contents |
|------|----------|
| `integration_benchmark_metrics.csv` | raw scib-metrics, one row per method |
| `integration_benchmark_metrics_scaled.csv` | min-max scaled (drives the figure) |
| `integration_benchmark_scib.svg` | scib-metrics summary table figure |
| `integration_cluster_recovery.csv` | Leiden-vs-C0C5 ARI/NMI per method |
| `integration_umap_grid.png` | per-method UMAP, colored by batch and by cluster |
| `integrated_eos_benchmark.h5ad` | object carrying every method's embedding |
| `benchmark_run.log` | per-method run log (timings, failures) |

## Interpreting against the manuscript narrative

The manuscript's primary representation is unintegrated; Harmony was reported as
over-correcting (low ARI vs C0–C5). This benchmark generalizes that single
comparison to the full method panel so the claim is "we tested N integration
tools, here is the batch-vs-bio trade-off for each," rather than "we tried
Harmony." Keep the response narrative consistent: report the trade-off table and
the cluster-recovery ARIs; do not also claim integration *recovers* subclusters
elsewhere.
