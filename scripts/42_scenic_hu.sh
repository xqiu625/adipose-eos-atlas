#!/bin/bash
#SBATCH --job-name=scenic_hu
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/scenic_hu_%a.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/scenic_hu_%a.err
#SBATCH --array=0-2
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=epyc

DEPOTS=(eWAT scWAT BAT)
DEPOT=${DEPOTS[$SLURM_ARRAY_TASK_ID]}

BASE="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
HU_BASE="${BASE}/scenic/hu_scenic"
DEPOT_DIR="${HU_BASE}/${DEPOT}"
mkdir -p "${DEPOT_DIR}"

PYTHON="/bigdata/godziklab/shared/Xinru/cellrank_env/bin/python"
PYSCENIC="/bigdata/godziklab/shared/Xinru/cellrank_env/bin/pyscenic"

# SCENIC resources
TF_LIST="/bigdata/godziklab/shared/Xinru/scenic_data/allTFs_mm.txt"
RANKING_DB="/bigdata/godziklab/shared/Xinru/scenic_data/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
MOTIF_ANN="/bigdata/godziklab/shared/Xinru/scenic_data/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

INPUT_LOOM="${DEPOT_DIR}/input.loom"
ADJ_CSV="${DEPOT_DIR}/adj.csv"
REG_CSV="${DEPOT_DIR}/regulons.csv"
AUC_CSV="${DEPOT_DIR}/auc_matrix.csv"

echo "========================================"
echo "pySCENIC per-depot: ${DEPOT}"
echo "========================================"

export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/cellrank_env/lib:/opt/linux/rocky/8.x/x86_64/pkgs/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
export DEPOT="${DEPOT}"

# --- Step 0: Export depot-specific loom ---
echo "Step 0: Creating loom for ${DEPOT}..."
if [ ! -f "$INPUT_LOOM" ]; then
    $PYTHON << 'PYEOF'
import os
import numpy as np
import pandas as pd
import h5py

depot = os.environ["DEPOT"]
base = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
depot_dir = os.path.join(base, "scenic", "hu_scenic", depot)
expr_csv = os.path.join(base, "scenic", "hu_expr_matrix.csv")
meta_csv = os.path.join(base, "scenic", "hu_depot_metadata.csv")

if not os.path.exists(expr_csv):
    print("ERROR: hu_expr_matrix.csv not found. Run script 41 first.")
    exit(1)

print(f"  Loading expression matrix...")
expr = pd.read_csv(expr_csv, index_col=0)
meta = pd.read_csv(meta_csv).set_index("barcode")

# Filter to this depot
depot_cells = meta.index[meta["depot"] == depot].tolist()
shared = [c for c in depot_cells if c in expr.columns]
expr_depot = expr[shared]
print(f"  {depot}: {expr_depot.shape[0]} genes x {expr_depot.shape[1]} cells")

# Create loom
loom_path = os.path.join(depot_dir, "input.loom")
mat = expr_depot.values.astype(np.float32)
genes = expr_depot.index.tolist()
barcodes = expr_depot.columns.tolist()

with h5py.File(loom_path, "w") as f:
    f.create_dataset("matrix", data=mat)
    dt = h5py.string_dtype()
    ra = f.create_group("row_attrs")
    ra.create_dataset("Gene", data=np.array(genes, dtype=object), dtype=dt)
    ca = f.create_group("col_attrs")
    ca.create_dataset("CellID", data=np.array(barcodes, dtype=object), dtype=dt)

print(f"  Created: {loom_path}")
PYEOF
    echo "  Loom done: $INPUT_LOOM"
else
    echo "  Loom exists, skipping"
fi

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
    echo "  cisTarget done: $REG_CSV"
else
    echo "  cisTarget exists, skipping"
fi

# --- Step 3: AUCell ---
echo "Step 3: AUCell..."
if [ ! -f "$AUC_CSV" ]; then
    $PYTHON << 'PYEOF'
import os
import sys
import numpy as np
import pandas as pd

depot = os.environ["DEPOT"]
depot_dir = f"/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/scenic/hu_scenic/{depot}"
reg_csv = os.path.join(depot_dir, "regulons.csv")
input_loom = os.path.join(depot_dir, "input.loom")
auc_csv = os.path.join(depot_dir, "auc_matrix.csv")

import h5py
with h5py.File(input_loom, "r") as f:
    mat = f["matrix"][:]
    genes = [g.decode() for g in f["row_attrs"]["Gene"][:]]
    barcodes = [b.decode() for b in f["col_attrs"]["CellID"][:]]
ex_mtx = pd.DataFrame(mat.T.astype(np.float32), index=barcodes, columns=genes)
print(f"  Expression: {ex_mtx.shape}")

from pyscenic.utils import load_motifs
from pyscenic.prune import df2regulons

try:
    motifs = load_motifs(reg_csv)
    regulons = df2regulons(motifs)
    print(f"  Regulons: {len(regulons)}")
except Exception as e:
    print(f"  Failed to load regulons: {e}")
    sys.exit(1)

if not regulons:
    print("  ERROR: No regulons found")
    sys.exit(1)

from pyscenic.aucell import aucell
auc_mtx = aucell(ex_mtx, regulons, num_workers=8)
auc_mtx.to_csv(auc_csv)
print(f"  Saved AUCell: {auc_csv}")
print(f"  {auc_mtx.shape[0]} cells x {auc_mtx.shape[1]} regulons")
PYEOF
    echo "  AUCell done: $AUC_CSV"
else
    echo "  AUCell exists, skipping"
fi

echo ""
echo "Done: ${DEPOT}"
