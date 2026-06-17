#!/usr/bin/env python
"""Nature-style broad-reference figure (Reviewer #1 — identity + subcluster reality).
Co-clustering our cells with published eosinophils (Hu, Gurtner) and our own
macrophages/mast/neutrophils/DC/NK/B as negative-control lineages."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.8, "legend.frameon": False,
})

um = pd.read_csv(HERE / "umap_cells_for_fig.csv")
q1 = pd.read_csv(HERE / "q1_identity_refeos_neighbor_fraction.csv")
cor = pd.read_csv(HERE / "q2_pseudobulk_correlation_C0C5_vs_Hu.csv", index_col=0)

CATCOL = {"Our eosinophils": "#C0392B", "Our macrophages": "#2c5aa0",
          "Our other lineages": "#BBBBBB", "Ref eos (Hu)": "#E69F00", "Ref eos (Gurtner)": "#009E73"}

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.42, wspace=0.34)

# A — joint UMAP by category
axA = fig.add_subplot(gs[0, :])
order = ["Our other lineages", "Our macrophages", "Ref eos (Gurtner)", "Ref eos (Hu)", "Our eosinophils"]
for cat in order:
    d = um[um.category == cat]
    axA.scatter(d.UMAP_1, d.UMAP_2, s=2.2, c=CATCOL[cat], label=f"{cat} (n={len(d)})",
                alpha=0.65 if cat != "Our other lineages" else 0.4, linewidths=0, rasterized=True)
axA.set_xticks([]); axA.set_yticks([]); axA.set_xlabel("UMAP 1"); axA.set_ylabel("UMAP 2")
axA.set_title("A  Joint embedding: our cells + published eosinophils (Hu, Gurtner)",
              loc="left", fontweight="bold", fontsize=8)
axA.legend(markerscale=4, fontsize=5.6, loc="best", ncol=1)

# B — Q1 identity bar
axB = fig.add_subplot(gs[1, 0])
q1 = q1.sort_values("mean_frac_ref_eos")
iseos = q1.lineage.str.match(r"^C[0-5]$")
axB.barh(q1.lineage, q1.mean_frac_ref_eos,
         color=[("#C0392B" if e else "#999999") for e in iseos], edgecolor="white", lw=0.4)
for y, (v, n) in enumerate(zip(q1.mean_frac_ref_eos, q1.n)):
    axB.text(v + 0.004, y, f"{v:.2f}", va="center", fontsize=5.2)
axB.set_xlabel("fraction of neighbours that are\nreference eosinophils")
axB.set_title("B  Identity: C0–C5 cluster with reference\neosinophils; other lineages do not",
              loc="left", fontweight="bold", fontsize=7.6)
axB.tick_params(axis="y", labelsize=6)

# C — Q2 correlation heatmap (C0-C5 x Hu Eos1-8)
axC = fig.add_subplot(gs[1, 1])
M = cor.values
im = axC.imshow(M, cmap="YlOrRd", vmin=0.3, vmax=0.85, aspect="auto")
axC.set_xticks(range(cor.shape[1])); axC.set_xticklabels(cor.columns, fontsize=6, rotation=45, ha="right")
axC.set_yticks(range(cor.shape[0])); axC.set_yticklabels(cor.index, fontsize=6)
for i in range(M.shape[0]):
    j = int(np.argmax(M[i]));
    for jj in range(M.shape[1]):
        axC.text(jj, i, f"{M[i,jj]:.2f}", ha="center", va="center", fontsize=4.3,
                 color="white" if M[i,jj] > 0.65 else "#333",
                 fontweight="bold" if jj == j else "normal")
    axC.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False, edgecolor="#0F4D92", lw=1.3))
axC.set_title("C  Subcluster: C0–C5 vs Hu Eos1–8\n(pseudobulk Spearman; best match boxed)",
              loc="left", fontweight="bold", fontsize=7.6)
cb = fig.colorbar(im, ax=axC, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=5)
cb.set_label("Spearman r", fontsize=5.5)

for ext in ("svg", "pdf", "png"):
    fig.savefig(HERE / f"fig_broad_reference.{ext}", dpi=500 if ext == "png" else None, bbox_inches="tight")
print("wrote fig_broad_reference.{svg,pdf,png}")
