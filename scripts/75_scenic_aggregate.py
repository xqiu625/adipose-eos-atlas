#!/usr/bin/env python
# Aggregate the integrated-atlas SCENIC AUCell matrix by cluster: mean regulon AUC
# per cluster + regulon specificity (RSS-like), for the Fig 5 regulon heatmap.
import h5py, numpy as np, pandas as pd
REDO = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
f = h5py.File(f"{REDO}/scenic/aucell.loom", "r")
auc = f["col_attrs/RegulonsAUC"][:]           # structured array: cells x regulons
regs = list(auc.dtype.names)
cellid = [c.decode() if isinstance(c, bytes) else str(c) for c in f["col_attrs/CellID"][:]]
A = pd.DataFrame({r: auc[r] for r in regs}, index=cellid)   # cells x regulons
f.close()

meta = pd.read_csv(f"{REDO}/fig2/fig2_percell.csv").set_index("barcode")
common = A.index.intersection(meta.index)
A = A.loc[common]; cl = meta.loc[common, "cluster"]
print("AUC:", A.shape, "regulons; clusters:", sorted(cl.unique()), flush=True)

# mean AUC per cluster (regulons x clusters)
mean_auc = A.groupby(cl.values).mean().T
mean_auc.to_csv(f"{REDO}/scenic/regulon_mean_auc_by_cluster.csv")

# regulon specificity score (RSS) per cluster via z-score across clusters
z = (mean_auc.sub(mean_auc.mean(1), axis=0)).div(mean_auc.std(1) + 1e-9, axis=0)
z.to_csv(f"{REDO}/scenic/regulon_zscore_by_cluster.csv")
# top specific regulons per cluster
print("\nTop cluster-specific regulons (z-score):")
for c in z.columns:
    top = z[c].sort_values(ascending=False).head(6).index.tolist()
    print(f"  {c}: {', '.join(r.replace('(+)','') for r in top)}", flush=True)
print("\nwrote regulon_mean_auc_by_cluster.csv + regulon_zscore_by_cluster.csv", flush=True)
