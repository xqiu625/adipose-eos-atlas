#!/bin/bash
#SBATCH --job-name=scenic_hu_all
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/scenic_hu_all_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/scenic_hu_all_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=4
#SBATCH --partition=highmem

BASE="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
HU_DIR="${BASE}/scenic/hu_scenic/combined"
mkdir -p "${HU_DIR}"

PYTHON="/bigdata/godziklab/shared/Xinru/cellrank_env/bin/python"
PYSCENIC="/bigdata/godziklab/shared/Xinru/cellrank_env/bin/pyscenic"

TF_LIST="/bigdata/godziklab/shared/Xinru/scenic_data/allTFs_mm.txt"
RANKING_DB="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
MOTIF_ANN="/bigdata/godziklab/shared/Xinru/scenic_data/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

INPUT_LOOM="${HU_DIR}/input.loom"
ADJ_CSV="${HU_DIR}/adj.csv"
REG_CSV="${HU_DIR}/regulons.csv"
AUC_CSV="${HU_DIR}/auc_matrix.csv"

echo "========================================"
echo "pySCENIC on ALL Hu eosinophils (combined)"
echo "========================================"

export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/cellrank_env/lib:/opt/linux/rocky/8.x/x86_64/pkgs/gcc/12.2.0/lib64:$LD_LIBRARY_PATH

# --- Step 0: Loom from full Hu expression ---
echo "Step 0: Creating combined loom..."
if [ ! -f "$INPUT_LOOM" ]; then
    $PYTHON << 'PYEOF'
import os, numpy as np, pandas as pd, h5py

base = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
expr = pd.read_csv(os.path.join(base, "scenic", "hu_expr_matrix.csv"), index_col=0)
print(f"  All Hu: {expr.shape[0]} genes x {expr.shape[1]} cells")

loom_path = os.path.join(base, "scenic", "hu_scenic", "combined", "input.loom")
with h5py.File(loom_path, "w") as f:
    f.create_dataset("matrix", data=expr.values.astype(np.float32))
    dt = h5py.string_dtype()
    f.create_group("row_attrs").create_dataset("Gene", data=np.array(expr.index.tolist(), dtype=object), dtype=dt)
    f.create_group("col_attrs").create_dataset("CellID", data=np.array(expr.columns.tolist(), dtype=object), dtype=dt)
print(f"  Created: {loom_path}")
PYEOF
    echo "  Loom done"
else
    echo "  Loom exists, skipping"
fi

# --- Step 1: GRN ---
echo "Step 1: GRN..."
if [ ! -f "$ADJ_CSV" ]; then
    $PYSCENIC grn --num_workers 8 --output "$ADJ_CSV" --method grnboost2 "$INPUT_LOOM" "$TF_LIST"
    echo "  GRN done"
else
    echo "  GRN exists, skipping"
fi

# --- Step 2: cisTarget ---
echo "Step 2: cisTarget..."
if [ ! -f "$REG_CSV" ]; then
    $PYSCENIC ctx "$ADJ_CSV" "$RANKING_DB" --annotations_fname "$MOTIF_ANN" --expression_mtx_fname "$INPUT_LOOM" --output "$REG_CSV" --num_workers 8 --mode "dask_multiprocessing"
    echo "  cisTarget done"
else
    echo "  cisTarget exists, skipping"
fi

# --- Step 3: AUCell ---
echo "Step 3: AUCell..."
if [ ! -f "$AUC_CSV" ]; then
    $PYTHON << 'PYEOF'
import os, sys, numpy as np, pandas as pd, h5py

d = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/scenic/hu_scenic/combined"
with h5py.File(os.path.join(d, "input.loom"), "r") as f:
    mat = f["matrix"][:]
    genes = [g.decode() for g in f["row_attrs"]["Gene"][:]]
    barcodes = [b.decode() for b in f["col_attrs"]["CellID"][:]]
ex_mtx = pd.DataFrame(mat.T.astype(np.float32), index=barcodes, columns=genes)
print(f"  Expression: {ex_mtx.shape}")

from pyscenic.utils import load_motifs
from pyscenic.prune import df2regulons
motifs = load_motifs(os.path.join(d, "regulons.csv"))
regulons = df2regulons(motifs)
print(f"  Regulons: {len(regulons)}")

from pyscenic.aucell import aucell
auc_mtx = aucell(ex_mtx, regulons, num_workers=8)
auc_mtx.to_csv(os.path.join(d, "auc_matrix.csv"))
print(f"  Saved: {auc_mtx.shape[0]} cells x {auc_mtx.shape[1]} regulons")
PYEOF
    echo "  AUCell done"
else
    echo "  AUCell exists, skipping"
fi

echo "Done: combined Hu SCENIC"
