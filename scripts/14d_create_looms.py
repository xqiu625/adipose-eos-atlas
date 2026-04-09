#!/usr/bin/env python3
"""
14d_create_looms.py
Convert per-cluster mtx exports from 14c to loom files for pySCENIC.
Uses h5py directly to bypass loompy (which has np.string_ numpy 2.0 bug).

Loom format: HDF5 with /matrix (genes x cells), /row_attrs/Gene, /col_attrs/CellID

SLURM: --partition=epyc --mem=32G --cpus-per-task=2 --time=00:30:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import h5py

BASE       = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR = os.path.join(BASE, "scenic_per_cluster")
CLUSTERS   = ["C0", "C1", "C2", "C3", "C4", "C5"]

print("14d_create_looms.py")
print("=" * 50)


def write_loom(path, X_cells_genes, obs_names, var_names):
    """Write a loom file using h5py (genes x cells matrix), bypassing loompy."""
    # loom convention: matrix is genes x cells
    if scipy.sparse.issparse(X_cells_genes):
        mat = X_cells_genes.T.toarray().astype(np.float32)   # genes x cells
    else:
        mat = X_cells_genes.T.astype(np.float32)

    n_genes, n_cells = mat.shape
    chunk_r = min(64, n_genes)
    chunk_c = min(64, n_cells)

    with h5py.File(path, "w") as f:
        f.create_dataset("matrix", data=mat,
                         chunks=(chunk_r, chunk_c),
                         compression="gzip", compression_opts=2)
        ra = f.create_group("row_attrs")
        ra.create_dataset("Gene",   data=np.array(var_names,  dtype="S"))
        ca = f.create_group("col_attrs")
        ca.create_dataset("CellID", data=np.array(obs_names,  dtype="S"))
        f.create_group("layers")
        f.create_group("col_graphs")
        f.create_group("row_graphs")
        f.attrs["CreationDate"]      = np.bytes_("2025-01-01")
        f.attrs["LOOM_SPEC_VERSION"] = np.bytes_("3.0.0")


for cl in CLUSTERS:
    cl_dir    = os.path.join(SCENIC_DIR, cl)
    mtx_path  = os.path.join(cl_dir, "matrix.mtx")
    loom_path = os.path.join(cl_dir, "input.loom")

    if not os.path.exists(mtx_path):
        print(f"  {cl}: matrix.mtx not found, skipping")
        continue

    if os.path.exists(loom_path):
        print(f"  {cl}: loom already exists, skipping")
        continue

    print(f"\n{cl}: loading matrix...")

    mat = scipy.io.mmread(mtx_path).tocsr()   # genes x cells (R writeMM is genes x cells)

    genes    = pd.read_csv(os.path.join(cl_dir, "features.tsv"),
                           header=None, names=["gene"])["gene"].tolist()
    barcodes = pd.read_csv(os.path.join(cl_dir, "barcodes.tsv"),
                           header=None, names=["barcode"])["barcode"].tolist()

    print(f"  {mat.shape[0]} genes x {mat.shape[1]} cells")

    # mat is genes x cells; write_loom expects cells x genes (it transposes internally)
    X_cells_genes = mat.T   # cells x genes

    write_loom(loom_path, X_cells_genes, obs_names=barcodes, var_names=genes)
    print(f"  Saved: {loom_path}")

print("\n14d complete.")
