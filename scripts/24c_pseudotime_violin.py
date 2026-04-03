#!/usr/bin/env python
"""
24c_pseudotime_violin.py
Fig 2H-style violin plot: per-group DPT ranked by median pseudotime.
Matches original manuscript approach (scripts 09+10):
  1. Load integrated Qiu+Gurtner h5ad from script 13b (global DPT already computed)
  2. Re-compute DPT per group (each group gets own neighbors, diffmap, DPT)
  3. Violin plot ranked by median per-group DPT

Groups:
  - Gurtner Bone Marrow, Gurtner Blood
  - Qiu eos clusters C0-C6

Input:  13b h5ad (two-way integration with global DPT)
Output: violin plot PDF
SLURM:  batch, 32G, 4 CPUs, 1h, cellrank_env
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Publication figure settings ---------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# --- Color palettes (match 00_theme.R) ---------------------------------------
EOS_COLORS = {
    "C0": "#6CA6CD", "C1": "#E84A5F", "C2": "#FF6F1F",
    "C3": "#CBAACB", "C4": "#B03A48", "C5": "#5BA4A4", "C6": "#2F3C7E",
}

EXTRA_COLORS = {
    "Gurtner\nBone Marrow": "#2E8B57",
    "Gurtner\nBlood": "#DAA520",
}

# --- Paths -------------------------------------------------------------------
BASE = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
INTEG_DIR = os.path.join(BASE, "integration")
H5AD_PATH = os.path.join(INTEG_DIR, "qiu_gurtner_cellrank.h5ad")
META_PATH = os.path.join(INTEG_DIR, "for_cellrank", "metadata.csv")
SCORES_PATH = os.path.join(INTEG_DIR, "cellrank_cytotrace_scores.csv")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# --- Load integrated AnnData from script 13b ---------------------------------
print("Loading integrated h5ad from script 13b ...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  h5ad obs columns: {adata.obs.columns.tolist()}")

# Load metadata CSV (has dataset, organ, original_cluster columns)
meta = pd.read_csv(META_PATH, index_col=0)
print(f"  metadata.csv: {len(meta)} cells, columns: {meta.columns.tolist()}")

# Load scores CSV (has ct_score, dpt_pseudotime)
scores = pd.read_csv(SCORES_PATH, index_col=0)
print(f"  scores.csv: {len(scores)} cells, columns: {scores.columns.tolist()}")

# Merge metadata + scores into adata.obs
for col in ["dataset", "organ", "original_cluster"]:
    if col in meta.columns:
        adata.obs[col] = meta.loc[adata.obs_names, col].values
        print(f"  Added '{col}' from metadata.csv")

for col in ["ct_score", "dpt_pseudotime"]:
    if col in scores.columns:
        adata.obs[col] = scores.loc[adata.obs_names, col].values
        print(f"  Added '{col}' from scores.csv")

print(f"  Final obs columns: {adata.obs.columns.tolist()}")

# --- Assign groups -----------------------------------------------------------
print("Assigning groups ...")

def assign_group(row):
    ds = str(row.get("dataset", ""))
    organ = str(row.get("organ", ""))
    cluster = str(row.get("original_cluster", ""))

    if "Gurtner" in ds:
        if organ == "Bone Marrow":
            return "Gurtner\nBone Marrow"
        elif organ == "Blood":
            return "Gurtner\nBlood"
        else:
            return None
    elif "Qiu" in ds:
        return cluster  # C0, C1, ..., C6
    else:
        return None

adata.obs["violin_group"] = adata.obs.apply(assign_group, axis=1)
adata = adata[adata.obs["violin_group"].notna()].copy()

groups = sorted(adata.obs["violin_group"].unique())
print(f"  Groups: {adata.obs['violin_group'].value_counts().to_dict()}")

# --- Per-group DPT (matching original manuscript script 10) ------------------
print("\nComputing per-group DPT (like original manuscript) ...")

group_dpt = {}  # group_name -> array of DPT values

for grp in groups:
    print(f"\n  --- {grp} ---")
    mask = adata.obs["violin_group"] == grp
    adata_grp = adata[mask].copy()
    n_cells = adata_grp.n_obs
    print(f"    n_cells = {n_cells}")

    if n_cells < 30:
        print(f"    WARNING: too few cells ({n_cells}), skipping per-group DPT")
        # Use global DPT if available
        if "dpt_pseudotime" in adata_grp.obs.columns:
            vals = adata_grp.obs["dpt_pseudotime"].astype(float).values
            vals = vals[~np.isinf(vals) & ~np.isnan(vals)]
            group_dpt[grp] = vals
            print(f"    Using global DPT (fallback): median={np.median(vals):.4f}")
        continue

    # Re-compute neighbors and diffusion map per group
    n_neighbors = min(30, n_cells - 1)
    n_pcs = min(30, adata_grp.obsm["X_harmony"].shape[1])

    sc.pp.neighbors(adata_grp, use_rep="X_harmony",
                    n_neighbors=n_neighbors, n_pcs=n_pcs)
    n_comps = min(15, n_cells - 2)
    sc.tl.diffmap(adata_grp, n_comps=n_comps)

    # Root cell: cell with highest CytoTRACE score (most progenitor-like)
    if "ct_score" in adata_grp.obs.columns:
        ct = adata_grp.obs["ct_score"].astype(float)
        root_cell = ct.idxmax()
        root_idx = adata_grp.obs_names.get_loc(root_cell)
        print(f"    Root cell: {root_cell} (ct_score={ct[root_cell]:.3f})")
    else:
        root_idx = 0
        print(f"    Root cell: index 0 (no ct_score available)")

    adata_grp.uns["iroot"] = root_idx
    n_dcs = min(10, n_comps)
    sc.tl.dpt(adata_grp, n_dcs=n_dcs)

    vals = adata_grp.obs["dpt_pseudotime"].values.astype(float)
    vals = vals[~np.isinf(vals) & ~np.isnan(vals)]
    group_dpt[grp] = vals
    print(f"    DPT range: {vals.min():.4f} - {vals.max():.4f}, "
          f"median={np.median(vals):.4f}, n_valid={len(vals)}")

# --- Order by median DPT (ascending) ----------------------------------------
medians = {g: np.median(v) for g, v in group_dpt.items()}
ordered_groups = sorted(medians, key=lambda g: medians[g])

print(f"\n  Group order by median per-group DPT:")
for g in ordered_groups:
    print(f"    {g:20s}  median={medians[g]:.4f}  n={len(group_dpt[g])}")

# --- Build color map ---------------------------------------------------------
all_colors = {}
all_colors.update(EXTRA_COLORS)
all_colors.update(EOS_COLORS)
colors = [all_colors.get(g, "#999999") for g in ordered_groups]

# --- Violin plot (Fig 2H style) ----------------------------------------------
print("\nGenerating violin plot ...")

fig, ax = plt.subplots(figsize=(max(10, len(ordered_groups) * 0.9 + 2), 6))

data_list = [group_dpt[g] for g in ordered_groups]

parts = ax.violinplot(data_list, positions=range(len(ordered_groups)),
                      showmeans=False, showmedians=False, showextrema=False)

for i, pc in enumerate(parts["bodies"]):
    pc.set_facecolor(colors[i])
    pc.set_edgecolor("black")
    pc.set_linewidth(0.8)
    pc.set_alpha(0.85)

bp = ax.boxplot(data_list, positions=range(len(ordered_groups)),
                widths=0.15, showfliers=False, patch_artist=True,
                medianprops=dict(color="white", linewidth=2),
                boxprops=dict(facecolor="black", linewidth=0.8),
                whiskerprops=dict(color="black", linewidth=0.8),
                capprops=dict(color="black", linewidth=0.8))

ax.set_xticks(range(len(ordered_groups)))
ax.set_xticklabels(ordered_groups, rotation=45, ha="right", fontsize=10,
                   fontweight="bold")
ax.set_ylabel("Diffusion Pseudotime", fontsize=13, fontweight="bold")
ax.set_title("Diffusion Pseudotime by Group\n(per-group DPT, root = most progenitor-like cell)",
             fontsize=14, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

out_path = os.path.join(FIG_DIR, "24c_dpt_violin.pdf")
fig.savefig(out_path)
plt.close()
print(f"\n  Saved: {out_path}")

# Save summary
summary_rows = []
for g in ordered_groups:
    v = group_dpt[g]
    summary_rows.append({
        "group": g.replace("\n", " "),
        "n": len(v),
        "mean": np.mean(v),
        "median": np.median(v),
        "std": np.std(v),
        "min": np.min(v),
        "max": np.max(v),
    })
summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(INTEG_DIR, "dpt_violin_summary.csv"), index=False)
print(f"  Saved: dpt_violin_summary.csv")

print("\nDone.")
