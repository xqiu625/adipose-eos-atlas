#!/bin/bash
#SBATCH --job-name=scenic_bat
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/scenic_bat_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/scenic_bat_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=epyc

DEPOT="BAT"
BASE="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
DEPOT_DIR="${BASE}/scenic/hu_scenic/${DEPOT}"
mkdir -p "${DEPOT_DIR}"

PYTHON="/bigdata/godziklab/shared/Xinru/cellrank_env/bin/python"
PYSCENIC="/bigdata/godziklab/shared/Xinru/cellrank_env/bin/pyscenic"

TF_LIST="/bigdata/godziklab/shared/Xinru/scenic_data/allTFs_mm.txt"
RANKING_DB="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
MOTIF_ANN="/bigdata/godziklab/shared/Xinru/scenic_data/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

INPUT_LOOM="${DEPOT_DIR}/input.loom"
ADJ_CSV="${DEPOT_DIR}/adj.csv"
REG_CSV="${DEPOT_DIR}/regulons.csv"
AUC_CSV="${DEPOT_DIR}/auc_matrix.csv"

echo "========================================"
echo "pySCENIC: BAT (standalone re-run)"
echo "========================================"

export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/cellrank_env/lib:/opt/linux/rocky/8.x/x86_64/pkgs/gcc/12.2.0/lib64:$LD_LIBRARY_PATH

# Loom should already exist from array job
if [ ! -f "$INPUT_LOOM" ]; then
    echo "ERROR: input.loom not found at $INPUT_LOOM"
    exit 1
fi

echo "Step 1: GRN..."
if [ ! -f "$ADJ_CSV" ] || [ ! -s "$ADJ_CSV" ]; then
    rm -f "$ADJ_CSV"
    $PYSCENIC grn --num_workers 8 --output "$ADJ_CSV" --method grnboost2 "$INPUT_LOOM" "$TF_LIST"
    echo "  GRN done: $(wc -l < $ADJ_CSV) lines"
else
    echo "  GRN exists ($(wc -l < $ADJ_CSV) lines), skipping"
fi

echo "Step 2: cisTarget..."
if [ ! -f "$REG_CSV" ] || [ ! -s "$REG_CSV" ]; then
    rm -f "$REG_CSV"
    $PYSCENIC ctx "$ADJ_CSV" "$RANKING_DB" --annotations_fname "$MOTIF_ANN" --expression_mtx_fname "$INPUT_LOOM" --output "$REG_CSV" --num_workers 8 --mode "dask_multiprocessing"
    echo "  cisTarget done: $(wc -l < $REG_CSV) lines"
else
    echo "  cisTarget exists, skipping"
fi

echo "Step 3: AUCell..."
if [ ! -f "$AUC_CSV" ]; then
    export DEPOT="${DEPOT}"
    $PYTHON << 'PYEOF'
import os, sys, numpy as np, pandas as pd, h5py

d = f"/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/scenic/hu_scenic/{os.environ['DEPOT']}"
with h5py.File(os.path.join(d, "input.loom"), "r") as f:
    mat = f["matrix"][:]
    genes = [g.decode() for g in f["row_attrs"]["Gene"][:]]
    barcodes = [b.decode() for b in f["col_attrs"]["CellID"][:]]
ex_mtx = pd.DataFrame(mat.T.astype(np.float32), index=barcodes, columns=genes)
print(f"  Expression: {ex_mtx.shape}")

from pyscenic.utils import load_motifs
from pyscenic.prune import df2regulons
try:
    motifs = load_motifs(os.path.join(d, "regulons.csv"))
    regulons = df2regulons(motifs)
    print(f"  Regulons: {len(regulons)}")
except Exception as e:
    print(f"  Failed: {e}")
    sys.exit(1)

from pyscenic.aucell import aucell
auc_mtx = aucell(ex_mtx, regulons, num_workers=8)
auc_mtx.to_csv(os.path.join(d, "auc_matrix.csv"))
print(f"  Saved: {auc_mtx.shape[0]} cells x {auc_mtx.shape[1]} regulons")
PYEOF
    echo "  AUCell done"
else
    echo "  AUCell exists, skipping"
fi

echo "Done: BAT"
