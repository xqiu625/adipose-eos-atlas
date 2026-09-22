#!/usr/bin/env python
# Export the ComBat-integrated eosinophils to a pySCENIC loom (raw counts,
# genes x cells, with the integrated cluster label). Run in cellrank_env (loompy).
import argparse, numpy as np, scanpy as sc, loompy as lp
from scipy.sparse import issparse
ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
ad = sc.read_h5ad(a.h5ad)
M = ad.layers["counts"] if "counts" in ad.layers else ad.X
M = M.toarray() if issparse(M) else np.asarray(M)        # cells x genes
M = M.T                                                  # genes x cells
row_attrs = {"Gene": np.asarray(ad.var_names, dtype=object)}
col_attrs = {"CellID": np.asarray(ad.obs_names, dtype=object),
             "cluster": np.asarray(ad.obs["new_cluster"].astype(str), dtype=object)}
lp.create(a.out, M, row_attrs, col_attrs)
print("wrote", a.out, "genes x cells:", M.shape, flush=True)
