#!/usr/bin/env python
"""
step2_regulon_generate.py
Prune GRNBoost2 adjacencies with cisTarget (new prune2df API) to get regulons.
"""

import argparse, logging, sys, time, pathlib
import pandas as pd
import scanpy as sc
from pyscenic.prune import prune2df
from pyscenic.utils import modules_from_adjacencies      # 关键
try:
    from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
except ModuleNotFoundError:
    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase  # type: ignore

# 有些版本没有 load_motif_annotations，不用它，prune2df 只要文件名
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_adata",  required=True)
    ap.add_argument("--input_adj",    required=True)
    ap.add_argument("--motif_tbl",    required=True)
    ap.add_argument("--rankings_db",  required=True, nargs="+")
    ap.add_argument("--output",       required=True)
    ap.add_argument("--threads",      type=int, default=8)
    ap.add_argument("--rank_threshold", type=int,   default=1500)
    ap.add_argument("--auc_threshold",  type=float, default=0.05)
    ap.add_argument("--nes_threshold",  type=float, default=3.0)
    ap.add_argument("--motif_similarity_fdr", type=float, default=0.001)
    args = ap.parse_args()

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )

    # 1) load data
    logging.info(f"Reading AnnData: {args.input_adata}")
    adata = sc.read_h5ad(args.input_adata)
    expr_df = pd.DataFrame(
        data=adata.X.A if hasattr(adata.X, "A") else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )

    logging.info(f"Reading adjacencies: {args.input_adj}")
    adj = pd.read_csv(args.input_adj)

    # 2) build modules
    logging.info("Building modules from adjacencies …")
    modules = list(modules_from_adjacencies(adj, expr_df))

    # 3) load ranking DBs (objects list)
    logging.info("Loading ranking DBs …")
    rnkdbs = [RankingDatabase(f, pathlib.Path(f).stem) for f in args.rankings_db]

    # 4) prune2df
    logging.info("Running prune2df …")
    t0 = time.time()
    pruned = prune2df(
        rnkdbs,
        modules,
        args.motif_tbl,
        rank_threshold=args.rank_threshold,
        auc_threshold=args.auc_threshold,
        nes_threshold=args.nes_threshold,
        motif_similarity_fdr=args.motif_similarity_fdr,
        client_or_address="dask_multiprocessing",
        num_workers=args.threads
    )
    logging.info(f"prune2df done in {(time.time()-t0)/60:.1f} min, rows={pruned.shape[0]}")

    # 5) turn pruned df → regulon table
    lc = {c.lower(): c for c in pruned.columns}
    tf_col  = lc.get("tf") or lc.get("transcription_factor") or lc.get("regulator") or lc.get("motif")
    tgt_col = lc.get("target") or lc.get("gene") or lc.get("target_gene") or lc.get("genes")
    if tf_col is None or tgt_col is None:
        raise ValueError(f"Cannot find TF/target columns in pruned_df: {pruned.columns.tolist()}")

    rows = []
    for tf, sub in pruned.groupby(tf_col):
        tgts = sorted(set(sub[tgt_col].astype(str)))
        rows.append({"regulon": f"{tf}({len(tgts)})", "tf": tf, "target_genes": ";".join(tgts)})

    pd.DataFrame(rows).sort_values("regulon").to_csv(out_path, index=False)
    logging.info(f"[✓] regulons written → {out_path}")

if __name__ == "__main__":
    main()