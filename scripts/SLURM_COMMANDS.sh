#!/bin/bash
# ==============================================================================
# SLURM_COMMANDS.sh
# Step-by-step SLURM submission commands for the adipose eosinophil atlas
#
# Run each block manually. Steps within a block can run in parallel.
# Wait for each block to complete (squeue -u $USER) before proceeding.
#
# Pipeline DAG:
#   00 → 01 → 02 → 03 → {04, 06} parallel
#   04 → {05, 07, 11, 12, 13, 14} parallel
#   13 → 13b (CellRank on Gurtner integration)
#   07 → 08
#   06 + 04 → 09 → 10
#   14 → 15 (grn → ctx → aucell) → 16
#   04 → {17, 18, 19, 20, 21} parallel
#   15 → 22
# ==============================================================================

SCRIPTS="/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
LOGS="/bigdata/godziklab/shared/Xinru/302025/06_LOGS"
CONDA="/bigdata/godziklab/shared/Xinru/cellrank_env"
BASE="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"

# SCENIC reference files
SCENIC_TF="/bigdata/godziklab/shared/Xinru/scenic_data/allTFs_mm.txt"
SCENIC_R1="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
SCENIC_R2="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
SCENIC_MA="/bigdata/godziklab/shared/Xinru/scenic_data/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

mkdir -p "$LOGS" "$BASE"

# ==============================================================================
# BLOCK 1: DecontX (start here)
# ==============================================================================
sbatch --job-name=00_decontx --partition=batch \
  --mem=128G --cpus-per-task=8 --time=4:00:00 \
  --output=${LOGS}/00_decontx_%j.out --error=${LOGS}/00_decontx_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/00_decontx_per_sample.R"

# >>> Wait for 00 to complete <<<

# ==============================================================================
# BLOCK 2: QC + doublet filter
# ==============================================================================
sbatch --job-name=01_qc --partition=batch \
  --mem=128G --cpus-per-task=8 --time=4:00:00 \
  --output=${LOGS}/01_qc_%j.out --error=${LOGS}/01_qc_%j.err \
  --wrap="source activate ${CONDA} && module load R/4.5.0 && Rscript ${SCRIPTS}/01_qc_doublet_filter.R"

# >>> Wait for 01 to complete <<<

# ==============================================================================
# BLOCK 3: All cells UMAP
# ==============================================================================
sbatch --job-name=02_umap --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/02_umap_%j.out --error=${LOGS}/02_umap_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/02_all_cells_umap.R"

# >>> Wait for 02 to complete <<<

# ==============================================================================
# BLOCK 4: Cell type annotation
# ==============================================================================
sbatch --job-name=03_annotate --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/03_annotate_%j.out --error=${LOGS}/03_annotate_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/03_annotate_celltypes.R"

# >>> Wait for 03 to complete <<<

# ==============================================================================
# BLOCK 5: Eos subcluster + Mac subtype (PARALLEL)
# ==============================================================================
sbatch --job-name=04_eos --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/04_eos_%j.out --error=${LOGS}/04_eos_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/04_eos_subcluster.R"

sbatch --job-name=06_mac --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/06_mac_%j.out --error=${LOGS}/06_mac_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/06_mac_subtype.R"

# >>> Wait for BOTH 04 and 06 to complete <<<

# ==============================================================================
# BLOCK 6: Eos downstream + combined UMAP (PARALLEL)
# All depend on 04; 09 also depends on 06 (which is done)
# ==============================================================================
sbatch --job-name=05_pseudo --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/05_pseudo_%j.out --error=${LOGS}/05_pseudo_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/05_eos_pseudobulk.R"

sbatch --job-name=07_deg --partition=batch \
  --mem=64G --cpus-per-task=4 --time=8:00:00 \
  --output=${LOGS}/07_deg_%j.out --error=${LOGS}/07_deg_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/07_deg_eos.R"

sbatch --job-name=09_combined --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/09_combined_%j.out --error=${LOGS}/09_combined_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/09_eos_mac_combined_umap.R"

sbatch --job-name=11_monocle --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/11_monocle_%j.out --error=${LOGS}/11_monocle_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/11_monocle3.R"

sbatch --job-name=12_hu --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/12_hu_%j.out --error=${LOGS}/12_hu_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/12_integrate_hu.R"

sbatch --job-name=13_gurtner --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/13_gurtner_%j.out --error=${LOGS}/13_gurtner_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/13_integrate_gurtner.R"

sbatch --job-name=14_scenic_exp --partition=batch \
  --mem=64G --cpus-per-task=4 --time=1:00:00 \
  --output=${LOGS}/14_scenic_exp_%j.out --error=${LOGS}/14_scenic_exp_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/14_scenic_export.R"

# Reviewer supplement scripts (all depend on 04)
sbatch --job-name=17_qc_supp --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/17_qc_supp_%j.out --error=${LOGS}/17_qc_supp_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/17_qc_supplement.R"

sbatch --job-name=18_clustree --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/18_clustree_%j.out --error=${LOGS}/18_clustree_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/18_clustree_resolution.R"

sbatch --job-name=19_batch --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/19_batch_%j.out --error=${LOGS}/19_batch_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/19_batch_effect_evidence.R"

sbatch --job-name=20_compos --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/20_compos_%j.out --error=${LOGS}/20_compos_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/20_compositional_analysis.R"

sbatch --job-name=21_ambient --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/21_ambient_%j.out --error=${LOGS}/21_ambient_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/21_ambient_rna_evidence.R"

# >>> Wait for 07 to complete, then: <<<

# ==============================================================================
# BLOCK 7: Pathway (depends on 07_deg)
# ==============================================================================
sbatch --job-name=08_pathway --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/08_pathway_%j.out --error=${LOGS}/08_pathway_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/08_pathway_eos.R"

# >>> Wait for 13 to complete, then: <<<

# ==============================================================================
# BLOCK 7b: CellRank on Gurtner integration (depends on 13)
# ==============================================================================
sbatch --job-name=13b_cellrank --partition=batch \
  --mem=64G --cpus-per-task=8 --time=4:00:00 \
  --output=${LOGS}/13b_cellrank_%j.out --error=${LOGS}/13b_cellrank_%j.err \
  --wrap="source activate ${CONDA} && python ${SCRIPTS}/13b_cellrank_gurtner.py"

# >>> Wait for 09 to complete, then: <<<

# ==============================================================================
# BLOCK 8: CellChat (depends on 09_combined)
# ==============================================================================
sbatch --job-name=10_cellchat --partition=batch \
  --mem=64G --cpus-per-task=4 --time=6:00:00 \
  --output=${LOGS}/10_cellchat_%j.out --error=${LOGS}/10_cellchat_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/10_cellchat.R"

# >>> Wait for 14 to complete, then: <<<

# ==============================================================================
# BLOCK 9: SCENIC pipeline (sequential: grn → ctx → aucell)
# ==============================================================================
# Step 15a: GRN inference (~8h)
sbatch --job-name=15a_grn --partition=batch \
  --mem=64G --cpus-per-task=16 --time=8:00:00 \
  --output=${LOGS}/15a_grn_%j.out --error=${LOGS}/15a_grn_%j.err \
  --wrap="source activate ${CONDA} && pyscenic grn ${BASE}/scenic/eos_for_scenic.h5ad ${SCENIC_TF} -o ${BASE}/scenic/adj.csv --num_workers 16"

# >>> Wait for 15a to complete <<<

# Step 15b: cisTarget (~4h)
sbatch --job-name=15b_ctx --partition=batch \
  --mem=64G --cpus-per-task=16 --time=4:00:00 \
  --output=${LOGS}/15b_ctx_%j.out --error=${LOGS}/15b_ctx_%j.err \
  --wrap="source activate ${CONDA} && pyscenic ctx ${BASE}/scenic/adj.csv ${SCENIC_R1} ${SCENIC_R2} --annotations_fname ${SCENIC_MA} --expression_mtx_fname ${BASE}/scenic/eos_for_scenic.h5ad --output ${BASE}/scenic/regulons.csv --num_workers 16"

# >>> Wait for 15b to complete <<<

# Step 15c: AUCell (~2h)
sbatch --job-name=15c_auc --partition=batch \
  --mem=64G --cpus-per-task=16 --time=2:00:00 \
  --output=${LOGS}/15c_auc_%j.out --error=${LOGS}/15c_auc_%j.err \
  --wrap="source activate ${CONDA} && pyscenic aucell ${BASE}/scenic/eos_for_scenic.h5ad ${BASE}/scenic/regulons.csv --output ${BASE}/scenic/eos_scenic_output.loom --num_workers 16"

# >>> Wait for 15c to complete <<<

# ==============================================================================
# BLOCK 10: SCENIC visualization (PARALLEL, both depend on 15c)
# ==============================================================================
sbatch --job-name=16_scenic_viz --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/16_scenic_viz_%j.out --error=${LOGS}/16_scenic_viz_%j.err \
  --wrap="source activate ${CONDA} && python ${SCRIPTS}/16_scenic_visualize.py"

sbatch --job-name=22_scenic_cond --partition=batch \
  --mem=64G --cpus-per-task=8 --time=2:00:00 \
  --output=${LOGS}/22_scenic_cond_%j.out --error=${LOGS}/22_scenic_cond_%j.err \
  --wrap="source activate ${CONDA} && python ${SCRIPTS}/22_scenic_by_condition.py"

# ==============================================================================
# BLOCK 11: Hu mast marker validation (depends on 04 + Hu data)
# ==============================================================================
sbatch --job-name=23_mast --partition=batch \
  --mem=64G --cpus-per-task=4 --time=2:00:00 \
  --output=${LOGS}/23_mast_%j.out --error=${LOGS}/23_mast_%j.err \
  --wrap="module load R/4.5.0; Rscript ${SCRIPTS}/23_hu_mast_marker_validation.R"

# ==============================================================================
# Monitor jobs:
#   squeue -u $USER
#   sacct --format=JobID,JobName,State,Elapsed,MaxRSS -j <JOB_IDs>
#   tail -f ${LOGS}/<step>_*.out
# ==============================================================================
