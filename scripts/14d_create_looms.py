#!/usr/bin/env python3
"""
14d_create_looms.py
Convert per-cluster mtx exports from 14c to loom files for pySCENIC.

SLURM: --partition=epyc --mem=32G --cpus-per-task=2 --time=00:30:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import anndata as ad

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
CLUSTERS   = ["C0", "C1", "C2", "C3", "C4", "C5"]

print("14d_create_looms.py")
print("=" * 50)

for cl in CLUSTERS:
    cl_dir   = os.path.join(SCENIC_DIR, cl)
    mtx_path = os.path.join(cl_dir, "matrix.mtx")
    loom_path = os.path.join(cl_dir, "input.loom")

    if not os.path.exists(mtx_path):
        print(f"  {cl}: matrix.mtx not found, skipping")
        continue

    if os.path.exists(loom_path):
        print(f"  {cl}: loom already exists, skipping")
        continue

    print(f"\n{cl}: loading matrix...")

    # Read sparse matrix (genes x cells, COO/CSC from R writeMM)
    mat = scipy.io.mmread(mtx_path).tocsc()   # genes x cells

    genes = pd.read_csv(os.path.join(cl_dir, "features.tsv"),
                        header=None, names=["gene"])["gene"].tolist()
    barcodes = pd.read_csv(os.path.join(cl_dir, "barcodes.tsv"),
                           header=None, names=["barcode"])["barcode"].tolist()
    meta = pd.read_csv(os.path.join(cl_dir, "metadata.csv"), index_col="barcode")

    print(f"  {mat.shape[0]} genes x {mat.shape[1]} cells")

    # Transpose to cells x genes for AnnData
    X = mat.T.tocsr()

    obs = meta.loc[barcodes]   # align order
    var = pd.DataFrame(index=genes)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = barcodes
    adata.var_names = genes

    adata.write_loom(loom_path, write_obsm_varm=False)
    print(f"  Saved: {loom_path}")

print("\n14d complete.")
