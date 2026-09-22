#!/bin/bash
# pySCENIC on the ComBat-integrated eosinophils (export loom -> GRN -> ctx -> AUCell).
set -euo pipefail
ENV=/bigdata/godziklab/shared/Xinru/cellrank_env
PY=$ENV/bin/python
PYS=$ENV/bin/pyscenic
S=/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts
RES=/bigdata/nairlab/shared/XQiu/EOS_project/R02_draft_v1/results/pySCENIC/resources
REDO=/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo
OUT=$REDO/scenic; mkdir -p "$OUT"
LOOM=$OUT/integrated_eos.loom
NW=16

echo "[0/4] strip uns (sc_scvi wrote the h5ad; cellrank anndata can't read uns/log1p)"
SCVI=/bigdata/godziklab/shared/Xinru/sc_scvi/bin/python
$SCVI -c "import scanpy as sc; a=sc.read_h5ad('$REDO/eos_combat010.h5ad'); a.uns.clear(); a.write_h5ad('$REDO/eos_combat010_clean.h5ad')"

echo "[1/4] export loom"
$PY $S/73_scenic_export_loom.py --h5ad $REDO/eos_combat010_clean.h5ad --out $LOOM

echo "[2/4] GRN (GRNBoost2)"
$PYS grn --num_workers $NW --seed 1 -o $OUT/adj.tsv --method grnboost2 $LOOM $RES/mm_mgi_tfs.txt

echo "[3/4] ctx (cisTarget -> regulons)"
$PYS ctx $OUT/adj.tsv \
  $RES/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  $RES/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  --annotations_fname $RES/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
  --expression_mtx_fname $LOOM --mode dask_multiprocessing \
  --output $OUT/regulons.csv --num_workers $NW

echo "[4/4] AUCell"
$PYS aucell $LOOM $OUT/regulons.csv --output $OUT/aucell.loom --num_workers $NW

echo "DONE pySCENIC integrated -> $OUT"
