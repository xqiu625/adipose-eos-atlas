#!/usr/bin/env python
"""Fig 6 — eosinophil clusters are not macrophages (dendrogram + correlation) and
the ComBat-integrated UMAP recovers C0-C5. Reviewer #1 / Meera's requested panels."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import linkage, dendrogram

HERE = Path(__file__).resolve().parent
UM = HERE.parent / "combat_umap"
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.8, "legend.frameon": False,
})
CLCOL = {"C0": "#E41A1C", "C1": "#377EB8", "C2": "#4DAF4A", "C3": "#984EA3", "C4": "#FF7F00", "C5": "#A65628"}
EOS_C = "#C0392B"; MAC_C = "#2c5aa0"

cor = pd.read_csv(HERE / "eos_mac_pseudobulk_correlation.csv", index_col=0)
um = pd.read_csv(UM / "combat_umap_cells.csv")
cent = pd.read_csv(UM / "combat_umap_centroids.csv", index_col=0)

fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.5, wspace=0.4)

# A — dendrogram (recompute linkage from 1-corr, Ward) so it renders cleanly
axA = fig.add_subplot(gs[0, 0])
labels = list(cor.index)
from scipy.spatial.distance import squareform
D = squareform(1 - cor.values, checks=False)
Z = linkage(D, method="ward")
dd = dendrogram(Z, labels=labels, ax=axA, color_threshold=0,
                above_threshold_color="#555", leaf_rotation=90)
for lbl in axA.get_xticklabels():
    t = lbl.get_text()
    lbl.set_color(EOS_C if t.startswith("Eos") else MAC_C)
    lbl.set_fontsize(6)
axA.set_title("A  Eosinophil clusters vs macrophage subtypes", loc="left", fontweight="bold", fontsize=8)
axA.set_ylabel("1 − Spearman (Ward)")
axA.text(0.02, 0.97, "eosinophils and macrophages split at the top node",
         transform=axA.transAxes, fontsize=5.4, color="#555", va="top")

# B — correlation heatmap, honest 0-1 scale, eos block vs mac block
axB = fig.add_subplot(gs[0, 1])
M = cor.values
im = axB.imshow(M, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
axB.set_xticks(range(len(labels))); axB.set_xticklabels(labels, fontsize=5, rotation=90)
axB.set_yticks(range(len(labels))); axB.set_yticklabels(labels, fontsize=5)
for i in range(len(labels)):
    for j in range(len(labels)):
        axB.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=3.6,
                 color="white" if M[i,j] > 0.6 else "#333")
# separate eos (6) from mac (3) blocks
axB.axhline(5.5, color="k", lw=1.0); axB.axvline(5.5, color="k", lw=1.0)
axB.set_title("B  Pseudobulk correlation (0–1)", loc="left", fontweight="bold", fontsize=8)
cb = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=5)
cb.set_label("Spearman r", fontsize=5.5)

# C — ComBat-integrated UMAP colored by cluster
axC = fig.add_subplot(gs[1, 0])
for c in ["C0", "C1", "C2", "C3", "C4", "C5"]:
    d = um[um.cluster == c]
    axC.scatter(d.UMAP_1, d.UMAP_2, s=2, c=CLCOL[c], label=c, alpha=0.7, linewidths=0, rasterized=True)
for c in cent.index:
    axC.text(cent.loc[c, "UMAP_1"], cent.loc[c, "UMAP_2"], c, fontsize=6, fontweight="bold", ha="center", va="center")
axC.set_xticks([]); axC.set_yticks([]); axC.set_xlabel("UMAP 1"); axC.set_ylabel("UMAP 2")
axC.set_title("C  ComBat-integrated UMAP — clusters", loc="left", fontweight="bold", fontsize=8)
axC.legend(markerscale=4, fontsize=5.4, loc="best", ncol=2)

# D — ComBat-integrated UMAP colored by condition (batch mixing)
axD = fig.add_subplot(gs[1, 1])
colcol = "condition" if "condition" in um.columns else "sample"
cats = sorted(um[colcol].astype(str).unique())
pal = plt.cm.tab10(np.linspace(0, 1, len(cats)))
for k, cc in enumerate(cats):
    d = um[um[colcol].astype(str) == cc]
    axD.scatter(d.UMAP_1, d.UMAP_2, s=2, color=pal[k], label=cc, alpha=0.6, linewidths=0, rasterized=True)
axD.set_xticks([]); axD.set_yticks([]); axD.set_xlabel("UMAP 1"); axD.set_ylabel("UMAP 2")
axD.set_title(f"D  ComBat-integrated UMAP — {colcol}", loc="left", fontweight="bold", fontsize=8)
axD.legend(markerscale=4, fontsize=4.4, loc="best", ncol=1)

for ext in ("svg", "pdf", "png"):
    fig.savefig(HERE / f"fig_eosmac_combat.{ext}", dpi=500 if ext == "png" else None, bbox_inches="tight")
print("wrote fig_eosmac_combat.{svg,pdf,png}")
