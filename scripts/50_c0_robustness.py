#!/usr/bin/env python
# ==============================================================================
# 50_c0_robustness.py   (Reviewer #1, Concern #1)
# The reviewer's exact question for the condition-skewed clusters: when they merge
# into another cluster under integration, do they RETAIN their own signature
# (= a real state that integration over-corrected), or do they just upregulate a
# handful of genes (= a transcriptional response)?  We test C0 (the least stable,
# integration purity ~0.5) and C2. For the integrated cluster each one mostly joins,
# we DE the original-cluster cells vs the rest of that integrated cluster.
# Run in the sc_scvi env (scanpy).
# ==============================================================================
import argparse, os, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--label_key", default="cluster")
ap.add_argument("--res", type=float, default=0.4)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

adata = sc.read_h5ad(a.h5ad)
adata.obs[a.label_key] = adata.obs[a.label_key].astype(str)
MHC2 = [g for g in ["H2-Aa","H2-Ab1","H2-Eb1","H2-DMa","H2-DMb1","Cd74","Ciita"] if g in adata.var_names]

rows, topgenes = [], {}
for emb in [e for e in ["X_pca_harmony","X_combat"] if e in adata.obsm]:
    sc.pp.neighbors(adata, use_rep=emb, n_neighbors=30, random_state=a.seed)
    sc.tl.leiden(adata, resolution=a.res, random_state=a.seed, key_added="intg")
    for t in ["C0", "C2"]:
        mask_t = (adata.obs[a.label_key] == t).values
        if mask_t.sum() == 0:
            continue
        top_intg = adata.obs.loc[mask_t, "intg"].value_counts().idxmax()
        sub = adata[adata.obs["intg"] == top_intg].copy()
        sub.obs["grp"] = np.where(sub.obs[a.label_key].values == t, t, "other")
        n_t = int((sub.obs["grp"] == t).sum()); n_o = int((sub.obs["grp"] == "other").sum())
        rec = {"embedding": emb, "orig_cluster": t, "joins_integrated_cluster": top_intg,
               "n_origin_here": n_t, "n_other_here": n_o}
        if n_t >= 20 and n_o >= 20:
            sc.tl.rank_genes_groups(sub, "grp", groups=[t], reference="other", method="wilcoxon")
            de = sc.get.rank_genes_groups_df(sub, group=t)
            sig = de[(de.pvals_adj < 0.05) & (de.logfoldchanges.abs() > 0.5)]
            rec["n_DE_retained"] = int(sig.shape[0])
            rec["n_DE_up"] = int((sig.logfoldchanges > 0).sum())
            rec["top_up_genes"] = ";".join(sig.sort_values("logfoldchanges", ascending=False).names.head(15))
            topgenes[f"{emb}_{t}"] = sig
            if t == "C0" and MHC2:
                m = de[de.names.isin(MHC2)]
                rec["MHC2_meanLFC"] = round(float(m.logfoldchanges.mean()), 2)
                rec["MHC2_genes_sig_up"] = int(((m.pvals_adj < 0.05) & (m.logfoldchanges > 0.5)).sum())
        else:
            rec["n_DE_retained"] = np.nan
        rows.append(rec)
        print(f"[{emb}] {t} -> integrated {top_intg}: origin {n_t} vs other {n_o}, "
              f"DE retained = {rec.get('n_DE_retained')}", flush=True)

summ = pd.DataFrame(rows)
summ.to_csv(os.path.join(a.outdir, "c0_c2_signature_retention.csv"), index=False)
for k, df in topgenes.items():
    df.to_csv(os.path.join(a.outdir, f"retained_DE_{k}.csv"), index=False)

print("\n=== signature retention (many DE retained = a real state integration over-merged;\n"
      "    ~0 DE = not a distinct state) ===")
print(summ.to_string(index=False), flush=True)
print("\nInterpretation key: for C0, MHC2_genes_sig_up > 0 means the MHC-II / antigen-"
      "presentation program is still differentially expressed even inside the integrated "
      "cluster — i.e., a genuine state that integration blurred.", flush=True)
