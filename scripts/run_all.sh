#!/bin/bash
# ==============================================================================
# run_all.sh
# Master SLURM submission script for the adipose eosinophil subtype atlas
# Submits all pipeline steps with dependency chaining via sbatch --parsable
# ==============================================================================
#
# Usage:
#   cd /path/to/adipose-eos-subtype-atlas
#   bash scripts/run_all.sh
#
# Pipeline DAG:
#
#   00_decontx
#       |
#   01_qc_doublet
#       |
#   02_all_cells_umap
#       |
#   03_annotate_celltypes
#       |
#   +--------+----------+
#   |        |          |
#  04_eos  06_mac    (parallel)
#   |        |
#  05_pseudo 09_combined <-- (needs both 04 and 06)
#   |        |
#  07_deg   10_cellchat <-- (needs 09)
#   |
#  08_pathway
#   |
#   +--------+----------+
#   |        |          |
#  11_monocle 12_hu  13_gurtner   (parallel, all need 04)
#
#  14_scenic_export  (needs 04)
#       |
#  15_scenic_run     (needs 14)
#       |
#  16_scenic_viz     (needs 15)
#
# ==============================================================================

set -euo pipefail

# -- Configuration -------------------------------------------------------------
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
LOGS="/bigdata/godziklab/shared/Xinru/302025/06_LOGS"
CONDA="/bigdata/godziklab/shared/Xinru/cellrank_env"

mkdir -p "${LOGS}"
mkdir -p "${BASE}"

echo "========================================"
echo "Adipose Eosinophil Atlas Pipeline"
echo "========================================"
echo "Scripts: ${SCRIPTS}"
echo "Output:  ${BASE}"
echo "Logs:    ${LOGS}"
echo "Conda:   ${CONDA}"
echo "========================================"
echo ""

# -- Helper function -----------------------------------------------------------
# submit_r  JOB_NAME TIME MEM CPUS SCRIPT [DEP_JOBID...]
submit_r() {
    local job_name="$1" time="$2" mem="$3" cpus="$4" script="$5"
    shift 5
    local dep_flag=""
    if [[ $# -gt 0 ]]; then
        local dep_ids=$(IFS=:; echo "$*")
        dep_flag="--dependency=afterok:${dep_ids}"
    fi
    sbatch --parsable \
        --partition=batch \
        --mem="${mem}" \
        --cpus-per-task="${cpus}" \
        --time="${time}" \
        --job-name="${job_name}" \
        --output="${LOGS}/${job_name}_%j.out" \
        --error="${LOGS}/${job_name}_%j.err" \
        ${dep_flag} \
        --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/${script}"
}

# submit_py  JOB_NAME TIME MEM CPUS SCRIPT [DEP_JOBID...]
submit_py() {
    local job_name="$1" time="$2" mem="$3" cpus="$4" script="$5"
    shift 5
    local dep_flag=""
    if [[ $# -gt 0 ]]; then
        local dep_ids=$(IFS=:; echo "$*")
        dep_flag="--dependency=afterok:${dep_ids}"
    fi
    sbatch --parsable \
        --partition=batch \
        --mem="${mem}" \
        --cpus-per-task="${cpus}" \
        --time="${time}" \
        --job-name="${job_name}" \
        --output="${LOGS}/${job_name}_%j.out" \
        --error="${LOGS}/${job_name}_%j.err" \
        ${dep_flag} \
        --wrap="source activate ${CONDA}; python ${SCRIPTS}/${script}"
}

# ==============================================================================
# Phase 0: Preprocessing
# ==============================================================================
echo "--- Phase 0: Preprocessing ---"

JOB_00=$(submit_r "00_decontx" "4:00:00" "128G" 8 "00_decontx_per_sample.R")
echo "  00_decontx:           ${JOB_00}"

JOB_01=$(submit_r "01_qc_doublet" "4:00:00" "128G" 8 "01_qc_doublet_filter.R" "${JOB_00}")
echo "  01_qc_doublet:        ${JOB_01}  -> after ${JOB_00}"

JOB_02=$(submit_r "02_umap" "2:00:00" "64G" 4 "02_all_cells_umap.R" "${JOB_01}")
echo "  02_all_cells_umap:    ${JOB_02}  -> after ${JOB_01}"

JOB_03=$(submit_r "03_annotate" "2:00:00" "64G" 4 "03_annotate_celltypes.R" "${JOB_02}")
echo "  03_annotate:          ${JOB_03}  -> after ${JOB_02}"

echo ""

# ==============================================================================
# Phase 1: Cell-type-specific analyses (parallel from 03)
# ==============================================================================
echo "--- Phase 1: Cell-type subsets (parallel) ---"

JOB_04=$(submit_r "04_eos" "2:00:00" "64G" 4 "04_eos_subcluster.R" "${JOB_03}")
echo "  04_eos_subcluster:    ${JOB_04}  -> after ${JOB_03}"

JOB_06=$(submit_r "06_mac" "2:00:00" "64G" 4 "06_mac_subcluster.R" "${JOB_03}")
echo "  06_mac_subcluster:    ${JOB_06}  -> after ${JOB_03}"

echo ""

# ==============================================================================
# Phase 2: Downstream of eosinophils (04)
# ==============================================================================
echo "--- Phase 2: Eosinophil downstream ---"

JOB_05=$(submit_r "05_pseudo" "2:00:00" "64G" 4 "05_eos_pseudobulk.R" "${JOB_04}")
echo "  05_eos_pseudobulk:    ${JOB_05}  -> after ${JOB_04}"

JOB_07=$(submit_r "07_deg" "2:00:00" "64G" 4 "07_deg.R" "${JOB_05}")
echo "  07_deg:               ${JOB_07}  -> after ${JOB_05}"

JOB_08=$(submit_r "08_pathway" "2:00:00" "64G" 4 "08_pathway.R" "${JOB_07}")
echo "  08_pathway:           ${JOB_08}  -> after ${JOB_07}"

echo ""

# ==============================================================================
# Phase 3: Combined eos + mac analyses
# ==============================================================================
echo "--- Phase 3: Combined (eos + mac) ---"

JOB_09=$(submit_r "09_combined" "2:00:00" "64G" 4 "09_combined.R" "${JOB_04}" "${JOB_06}")
echo "  09_combined:          ${JOB_09}  -> after ${JOB_04},${JOB_06}"

JOB_10=$(submit_r "10_cellchat" "2:00:00" "64G" 4 "10_cellchat.R" "${JOB_09}")
echo "  10_cellchat:          ${JOB_10}  -> after ${JOB_09}"

echo ""

# ==============================================================================
# Phase 4: Trajectory & integration (parallel, all need 04)
# ==============================================================================
echo "--- Phase 4: Trajectory & integration (parallel) ---"

JOB_11=$(submit_r "11_monocle" "2:00:00" "64G" 4 "11_monocle.R" "${JOB_04}")
echo "  11_monocle:           ${JOB_11}  -> after ${JOB_04}"

JOB_12=$(submit_r "12_hu" "2:00:00" "64G" 4 "12_hu.R" "${JOB_04}")
echo "  12_hu:                ${JOB_12}  -> after ${JOB_04}"

JOB_13=$(submit_r "13_gurtner" "2:00:00" "64G" 4 "13_gurtner.R" "${JOB_04}")
echo "  13_gurtner:           ${JOB_13}  -> after ${JOB_04}"

echo ""

# ==============================================================================
# Phase 5: SCENIC (needs 04)
# ==============================================================================
echo "--- Phase 5: SCENIC ---"

JOB_14=$(submit_r "14_scenic_export" "1:00:00" "64G" 4 "14_scenic_export.R" "${JOB_04}")
echo "  14_scenic_export:     ${JOB_14}  -> after ${JOB_04}"

# 15_scenic_run: inline the 3 pySCENIC steps here so dependency chain is exact.
# (15_scenic_run.sh can still be used standalone.)
SCENIC_TF="/bigdata/godziklab/shared/Xinru/scenic_data/allTFs_mm.txt"
SCENIC_R1="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
SCENIC_R2="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
SCENIC_MA="/bigdata/godziklab/shared/Xinru/scenic_data/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

JOB_15_GRN=$(sbatch --parsable \
    --partition=batch --mem=64G --cpus-per-task=16 --time=8:00:00 \
    --job-name="15_scenic_grn" \
    --output="${LOGS}/15_scenic_grn_%j.out" \
    --error="${LOGS}/15_scenic_grn_%j.err" \
    --dependency=afterok:${JOB_14} \
    --wrap="source activate ${CONDA} && pyscenic grn ${BASE}/scenic/eos_for_scenic.h5ad ${SCENIC_TF} -o ${BASE}/scenic/adj.csv --num_workers 16"
)
echo "  15a_scenic_grn:       ${JOB_15_GRN}  -> after ${JOB_14}  [8h]"

JOB_15_CTX=$(sbatch --parsable \
    --partition=batch --mem=64G --cpus-per-task=16 --time=4:00:00 \
    --job-name="15_scenic_ctx" \
    --output="${LOGS}/15_scenic_ctx_%j.out" \
    --error="${LOGS}/15_scenic_ctx_%j.err" \
    --dependency=afterok:${JOB_15_GRN} \
    --wrap="source activate ${CONDA} && pyscenic ctx ${BASE}/scenic/adj.csv ${SCENIC_R1} ${SCENIC_R2} --annotations_fname ${SCENIC_MA} --expression_mtx_fname ${BASE}/scenic/eos_for_scenic.h5ad --output ${BASE}/scenic/regulons.csv --num_workers 16"
)
echo "  15b_scenic_ctx:       ${JOB_15_CTX}  -> after ${JOB_15_GRN}  [4h]"

JOB_15_AUC=$(sbatch --parsable \
    --partition=batch --mem=64G --cpus-per-task=16 --time=2:00:00 \
    --job-name="15_scenic_aucell" \
    --output="${LOGS}/15_scenic_aucell_%j.out" \
    --error="${LOGS}/15_scenic_aucell_%j.err" \
    --dependency=afterok:${JOB_15_CTX} \
    --wrap="source activate ${CONDA} && pyscenic aucell ${BASE}/scenic/eos_for_scenic.h5ad ${BASE}/scenic/regulons.csv --output ${BASE}/scenic/eos_scenic_output.loom --num_workers 16"
)
echo "  15c_scenic_aucell:    ${JOB_15_AUC}  -> after ${JOB_15_CTX}  [2h]"

JOB_16=$(submit_py "16_scenic_viz" "2:00:00" "64G" 4 "16_scenic_visualize.py" "${JOB_15_AUC}")
echo "  16_scenic_visualize:  ${JOB_16}  -> after ${JOB_15_AUC}"

echo ""

# ==============================================================================
# Summary
# ==============================================================================
echo "========================================"
echo "All jobs submitted"
echo "========================================"
echo ""
echo "Job ID Summary:"
echo "  00_decontx:           ${JOB_00}"
echo "  01_qc_doublet:        ${JOB_01}"
echo "  02_all_cells_umap:    ${JOB_02}"
echo "  03_annotate:          ${JOB_03}"
echo "  04_eos_subcluster:    ${JOB_04}"
echo "  05_eos_pseudobulk:    ${JOB_05}"
echo "  06_mac_subcluster:    ${JOB_06}"
echo "  07_deg:               ${JOB_07}"
echo "  08_pathway:           ${JOB_08}"
echo "  09_combined:          ${JOB_09}"
echo "  10_cellchat:          ${JOB_10}"
echo "  11_monocle:           ${JOB_11}"
echo "  12_hu:                ${JOB_12}"
echo "  13_gurtner:           ${JOB_13}"
echo "  14_scenic_export:     ${JOB_14}"
echo "  15a_scenic_grn:       ${JOB_15_GRN}"
echo "  15b_scenic_ctx:       ${JOB_15_CTX}"
echo "  15c_scenic_aucell:    ${JOB_15_AUC}"
echo "  16_scenic_visualize:  ${JOB_16}"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct --format=JobID,JobName,State,Elapsed,MaxRSS -j ${JOB_00},${JOB_01},${JOB_02},${JOB_03},${JOB_04},${JOB_05},${JOB_06},${JOB_07},${JOB_08},${JOB_09},${JOB_10},${JOB_11},${JOB_12},${JOB_13},${JOB_14},${JOB_15_GRN},${JOB_15_CTX},${JOB_15_AUC},${JOB_16}"
echo ""
echo "Logs: ${LOGS}/"
