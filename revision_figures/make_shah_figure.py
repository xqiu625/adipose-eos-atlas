#!/usr/bin/env python
"""Nature-style Shah-concordance figure (Reviewer #1 — fixation/representativeness).
Our scRNA eosinophil pseudobulk vs published bulk-sorted eosinophils (Shah 2023,
GSE217551). High, eosinophil-specific concordance = our FLEX/fixed transcriptome
is representative, not a fixation artifact."""
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
gl = pd.read_csv(HERE / "shah_concordance_genelevel.csv")     # gene, ours, shah_gonadal, shah_blood
sp = pd.read_csv(HERE / "shah_specificity_control.csv")

fig = plt.figure(figsize=(7.2, 3.1))
gs = fig.add_gridspec(1, 2, wspace=0.34, width_ratios=[1, 1.05])

# A — scatter: our eos pseudobulk vs Shah gonadal-adipose bulk eos
axA = fig.add_subplot(gs[0, 0])
axA.scatter(gl.shah_gonadal, gl.ours, s=3, c="#C0392B", alpha=0.25, linewidths=0, rasterized=True)
r = gl[["ours", "shah_gonadal"]].corr(method="spearman").iloc[0, 1]
lim = [0, max(gl.ours.max(), gl.shah_gonadal.max()) * 1.02]
axA.plot(lim, lim, "--", color="#999", lw=0.7)
axA.set_xlim(lim); axA.set_ylim(lim)
axA.set_xlabel("Shah 2023 bulk — gonadal adipose eos\nlog(CPM+1)")
axA.set_ylabel("Our scRNA eosinophils (pseudobulk)\nlog(CPM+1)")
axA.set_title("A  Our eosinophil transcriptome vs\nindependent bulk eosinophils", loc="left", fontweight="bold", fontsize=7.8)
axA.text(0.04, 0.93, f"Spearman r = {r:.2f}\n(all common genes)", transform=axA.transAxes, fontsize=6.5, va="top")

# B — specificity bars: eos vs mac vs neutro, all-gene & variable-gene
axB = fig.add_subplot(gs[0, 1])
labs = sp.our_pseudobulk.str.replace(" (C0-C5)", "", regex=False)
x = np.arange(len(sp)); w = 0.38
cols = ["#C0392B", "#2c5aa0", "#E69F00"]
axB.bar(x - w/2, sp.vs_Shah_eos_allgenes, w, color=cols, edgecolor="white", lw=0.4, label="all genes")
axB.bar(x + w/2, sp.vs_Shah_eos_top3000var, w, color=cols, alpha=0.5, hatch="///", edgecolor="white", lw=0.4, label="variable genes")
for i in range(len(sp)):
    axB.text(i - w/2, sp.vs_Shah_eos_allgenes[i] + 0.01, f"{sp.vs_Shah_eos_allgenes[i]:.2f}", ha="center", fontsize=5.2)
    axB.text(i + w/2, sp.vs_Shah_eos_top3000var[i] + 0.01, f"{sp.vs_Shah_eos_top3000var[i]:.2f}", ha="center", fontsize=5.2)
axB.set_xticks(x); axB.set_xticklabels(labs, fontsize=6.5)
axB.set_ylabel("Spearman r vs Shah bulk eosinophils"); axB.set_ylim(0, 1)
axB.set_title("B  Concordance is eosinophil-specific", loc="left", fontweight="bold", fontsize=7.8)
from matplotlib.patches import Patch
axB.legend(handles=[Patch(facecolor="#777", label="all genes"),
                    Patch(facecolor="#777", alpha=0.5, hatch="///", label="variable genes (housekeeping removed)")],
           fontsize=5.2, loc="upper right")
axB.text(0.5, -0.30, "on variable genes our eosinophils (0.66) >> macrophages (0.24)",
         transform=axB.transAxes, ha="center", fontsize=5.4, color="#666")

for ext in ("svg", "pdf", "png"):
    fig.savefig(HERE / f"fig_shah_concordance.{ext}", dpi=500 if ext == "png" else None, bbox_inches="tight")
print("wrote fig_shah_concordance.{svg,pdf,png}")
