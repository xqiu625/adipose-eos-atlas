#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --job-name=int_bench
#SBATCH --output=int_bench_%j.out

# ==============================================================================
# 46b_run_integration_benchmark.sh
# Driver for the integration-method benchmark.
#
#   Step 1 (45, CPU): export eos_subclustered.rds -> eos_for_integration.h5ad
#   Step 2 (46, GPU): run + score every integration tool with scib-metrics
#
# GPU is requested for scVI / scANVI. Harmony, Scanorama, ComBat, BBKNN are CPU;
# to run only those, drop scvi,scanvi from --methods and use --partition=batch.
#
# One-time env setup (matches the AI4Bio batch-correction module):
#   conda activate /bigdata/godziklab/shared/Xinru/cellrank_env
#   pip install harmonypy scanorama bbknn scvi-tools scib-metrics
#
# Submit:  sbatch scripts/46b_run_integration_benchmark.sh
# ==============================================================================

set -euo pipefail

BASE=/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas
SCRIPTS=/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts
ENV=/bigdata/godziklab/shared/Xinru/cellrank_env
BENCH=$BASE/integration/benchmark
H5AD=$BENCH/eos_for_integration.h5ad

source activate "$ENV" || conda activate "$ENV"

# ── Step 1: export (skip if h5ad already present) ────────────────────────────
if [ ! -f "$H5AD" ]; then
  echo "[45] Exporting eos object to h5ad..."
  Rscript "$SCRIPTS/45_integration_export_h5ad.R"
else
  echo "[45] $H5AD exists, skipping export."
fi

# ── Step 2: benchmark. Default batch = GEM reaction (technical multiplex) ─────
echo "[46] Running integration benchmark (batch=gem_reaction)..."
python "$SCRIPTS/46_integration_benchmark.py" \
  --h5ad "$H5AD" --outdir "$BENCH" \
  --batch_key gem_reaction --label_key cluster \
  --methods all

# ── Optional sensitivity: per-library batch (orig.ident) ─────────────────────
# echo "[46] Re-running with batch=sample..."
# python "$SCRIPTS/46_integration_benchmark.py" \
#   --h5ad "$H5AD" --outdir "$BENCH/by_sample" \
#   --batch_key sample --label_key cluster --methods all

echo "Done. Results in $BENCH"
