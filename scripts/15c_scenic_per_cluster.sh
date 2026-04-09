#!/bin/bash
#SBATCH --job-name=scenic_per_cluster
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/scenic_per_cluster_%a.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/scenic_per_cluster_%a.err
#SBATCH --array=0-5
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=epyc

CLUSTERS=(C0 C1 C2 C3 C4 C5)
CLUSTER=${CLUSTERS[$SLURM_ARRAY_TASK_ID]}

BASE="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/scenic_per_cluster"
CLUSTER_DIR="${BASE}/${CLUSTER}"
PYTHON="/bigdata/godziklab/shared/Xinru/cellrank_env/bin/python"
PYSCENIC="/bigdata/godziklab/shared/Xinru/cellrank_env/bin/pyscenic"

# SCENIC resources
TF_LIST="/bigdata/godziklab/shared/Xinru/scenic_data/allTFs_mm.txt"
RANKING_DB="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
MOTIF_ANN="/bigdata/godziklab/shared/Xinru/scenic_data/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

INPUT_LOOM="${CLUSTER_DIR}/input.loom"
ADJ_CSV="${CLUSTER_DIR}/adj.csv"
REG_CSV="${CLUSTER_DIR}/regulons.csv"
AUC_CSV="${CLUSTER_DIR}/auc_matrix.csv"

echo "========================================"
echo "pySCENIC per-cluster: ${CLUSTER}"
echo "========================================"
echo "Input: ${INPUT_LOOM}"

# --- Step 1: GRN (GRNBoost2) ---
echo "Step 1: GRN..."
if [ ! -f "$ADJ_CSV" ]; then
    $PYSCENIC grn \
        --num_workers 8 \
        --output "$ADJ_CSV" \
        --method grnboost2 \
        "$INPUT_LOOM" "$TF_LIST"
    echo "  GRN done: $ADJ_CSV"
else
    echo "  GRN exists, skipping"
fi

# --- Step 2: cisTarget ---
echo "Step 2: cisTarget..."
if [ ! -f "$REG_CSV" ]; then
    $PYSCENIC ctx \
        "$ADJ_CSV" "$RANKING_DB" \
        --annotations_fname "$MOTIF_ANN" \
        --expression_mtx_fname "$INPUT_LOOM" \
        --output "$REG_CSV" \
        --num_workers 8 \
        --mode "dask_multiprocessing"
    echo "  ctx done: $REG_CSV"
else
    echo "  ctx exists, skipping"
fi

# --- Step 3: AUCell ---
echo "Step 3: AUCell..."
export CLUSTER="${CLUSTER}"
if [ ! -f "$AUC_CSV" ]; then
    $PYTHON << 'PYEOF'
import os, sys
import pandas as pd
from pyscenic.export import export2loom
from pyscenic.aucell import aucell
from pyscenic.utils import load_motifs

cluster = os.environ.get("CLUSTER", "")
base = f"/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/scenic_per_cluster/{cluster}"
reg_csv = f"{base}/regulons.csv"
input_loom = f"{base}/input.loom"
auc_csv = f"{base}/auc_matrix.csv"

import anndata as ad
import numpy as np

# Load expression matrix directly via h5py (bypasses loompy numpy2.0 np.string_ bug)
import h5py, scipy.sparse as _sp
with h5py.File(input_loom, "r") as f:
    mat      = f["matrix"][:]                          # genes x cells
    genes    = [g.decode() for g in f["row_attrs"]["Gene"][:]]
    barcodes = [b.decode() for b in f["col_attrs"]["CellID"][:]]
ex_mtx = pd.DataFrame(mat.T.astype(np.float32), index=barcodes, columns=genes)

print(f"  Expression: {ex_mtx.shape}")

# Load regulons
from ctxcore.genesig import GeneSignature
from pyscenic.utils import load_motifs
import csv

regulons = {}
try:
    motifs = load_motifs(reg_csv)
    for motif in motifs.itertuples():
        tf = motif.TF
        genes = set(g for g, _ in motif.TargetGenes)
        if tf not in regulons:
            regulons[tf] = set()
        regulons[tf].update(genes)
except Exception as e:
    print(f"  load_motifs failed: {e}, trying manual parse")
    with open(reg_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            tf = row.get("TF", row.get("", ""))
            tg = row.get("TargetGenes", "")
            if tf and tg:
                try:
                    import ast
                    genes = [g[0] for g in ast.literal_eval(tg)]
                    regulons.setdefault(tf, set()).update(genes)
                except:
                    pass

gene_sigs = [GeneSignature(name=tf, gene2weight={g: 1.0 for g in genes})
             for tf, genes in regulons.items() if len(genes) >= 5]
print(f"  Regulons: {len(gene_sigs)}")

# AUCell
from pyscenic.aucell import aucell
auc_mtx = aucell(ex_mtx, gene_sigs, num_workers=8)
auc_mtx.to_csv(auc_csv)
print(f"  Saved AUCell: {auc_csv}")
PYEOF
    echo "  AUCell done: $AUC_CSV"
else
    echo "  AUCell exists, skipping"
fi

echo ""
echo "Done: ${CLUSTER}"
