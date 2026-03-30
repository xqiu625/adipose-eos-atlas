#!/bin/bash
# ==============================================================================
# 15_scenic_run.sh
# pySCENIC pipeline: GRN inference -> cisTarget -> AUCell
# Submits 3 SLURM jobs with dependency chaining (--dependency=afterok)
# ==============================================================================
#
# Usage:
#   bash 15_scenic_run.sh
#
# Prerequisites:
#   - 14_scenic_export.R must have completed (eos_for_scenic.h5ad exists)
#   - pySCENIC, arboreto installed in conda env
#   - Ranking databases and TF list available
# ==============================================================================

set -euo pipefail

# -- Configuration -------------------------------------------------------------
BASE="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
LOGS="/bigdata/godziklab/shared/Xinru/302025/06_LOGS"
CONDA="/bigdata/godziklab/shared/Xinru/cellrank_env"

# Input/output paths
H5AD="${BASE}/scenic/eos_for_scenic.h5ad"
ADJ="${BASE}/scenic/adj.csv"
REGULONS="${BASE}/scenic/regulons.csv"
LOOM_OUT="${BASE}/scenic/eos_scenic_output.loom"

# Reference databases
TF_LIST="/bigdata/godziklab/shared/Xinru/scenic_data/allTFs_mm.txt"
RANK_10KB="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
RANK_500BP="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
MOTIF_ANNOT="/bigdata/godziklab/shared/Xinru/scenic_data/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

# -- Verify input exists -------------------------------------------------------
if [[ ! -f "${H5AD}" ]]; then
    echo "ERROR: Input file not found: ${H5AD}"
    echo "       Run 14_scenic_export.R first."
    exit 1
fi

mkdir -p "${BASE}/scenic"
mkdir -p "${LOGS}"

echo "========================================"
echo "15_scenic_run.sh — pySCENIC pipeline"
echo "========================================"
echo ""
echo "Input h5ad:  ${H5AD}"
echo "Output loom: ${LOOM_OUT}"
echo "Log dir:     ${LOGS}"
echo ""

# -- Job 1: GRN inference (GRNBoost2 via arboreto) ----------------------------
echo "Submitting Job 1: GRN inference (GRNBoost2)..."
JOB1=$(sbatch --parsable \
    --partition=batch \
    --mem=64G \
    --cpus-per-task=16 \
    --time=8:00:00 \
    --job-name=scenic_grn \
    --output="${LOGS}/scenic_grn_%j.out" \
    --error="${LOGS}/scenic_grn_%j.err" \
    --wrap="source activate ${CONDA} && \
export DASK_DATAFRAME__QUERY_PLANNING=False && \
pyscenic grn \
  ${H5AD} \
  ${TF_LIST} \
  -o ${ADJ} \
  --num_workers 16 && \
echo 'GRN inference complete.'"
)
echo "  Job 1 (GRN):     ${JOB1}"

# -- Job 2: cisTarget motif enrichment -----------------------------------------
echo "Submitting Job 2: cisTarget (depends on Job 1)..."
JOB2=$(sbatch --parsable \
    --dependency=afterok:${JOB1} \
    --partition=batch \
    --mem=64G \
    --cpus-per-task=16 \
    --time=4:00:00 \
    --job-name=scenic_ctx \
    --output="${LOGS}/scenic_ctx_%j.out" \
    --error="${LOGS}/scenic_ctx_%j.err" \
    --wrap="source activate ${CONDA} && \
pyscenic ctx \
  ${ADJ} \
  ${RANK_10KB} \
  ${RANK_500BP} \
  --annotations_fname ${MOTIF_ANNOT} \
  --expression_mtx_fname ${H5AD} \
  --output ${REGULONS} \
  --num_workers 16 && \
echo 'cisTarget complete.'"
)
echo "  Job 2 (ctx):     ${JOB2} (depends on ${JOB1})"

# -- Job 3: AUCell scoring -----------------------------------------------------
echo "Submitting Job 3: AUCell (depends on Job 2)..."
JOB3=$(sbatch --parsable \
    --dependency=afterok:${JOB2} \
    --partition=batch \
    --mem=64G \
    --cpus-per-task=16 \
    --time=2:00:00 \
    --job-name=scenic_aucell \
    --output="${LOGS}/scenic_aucell_%j.out" \
    --error="${LOGS}/scenic_aucell_%j.err" \
    --wrap="source activate ${CONDA} && \
pyscenic aucell \
  ${H5AD} \
  ${REGULONS} \
  --output ${LOOM_OUT} \
  --num_workers 16 && \
echo 'AUCell scoring complete.'"
)
echo "  Job 3 (aucell):  ${JOB3} (depends on ${JOB2})"

# -- Summary -------------------------------------------------------------------
echo ""
echo "========================================"
echo "All pySCENIC jobs submitted"
echo "========================================"
echo "  Job 1 (GRN):     ${JOB1}  [8h, 64G, 16 CPUs]"
echo "  Job 2 (ctx):     ${JOB2}  [4h, 64G, 16 CPUs] -> after ${JOB1}"
echo "  Job 3 (aucell):  ${JOB3}  [2h, 64G, 16 CPUs] -> after ${JOB2}"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${JOB1},${JOB2},${JOB3} --format=JobID,JobName,State,Elapsed,MaxRSS"
echo ""
echo "Output files:"
echo "  ${ADJ}"
echo "  ${REGULONS}"
echo "  ${LOOM_OUT}"
