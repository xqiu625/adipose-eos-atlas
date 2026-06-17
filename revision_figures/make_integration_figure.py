#!/usr/bin/env python
"""Figure 3 — integration recovers the sex-dependent compositional conclusions.
At matched granularity (res 1.0 -> 6 clusters), ComBat (the higher-bio-conservation
method the reviewer recommended) re-separates the C0-C5 clusters and reproduces the
male-enriched MHC-II (C0) and female-enriched migratory (C1) shifts; Harmony loses
C0's male enrichment. Magnitudes are attenuated (confounded Flex design)."""
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
CLCOL = {"C0": "#E41A1C", "C1": "#377EB8", "C2": "#4DAF4A", "C3": "#984EA3", "C4": "#FF7F00", "C5": "#A65628"}
RED, BLUE, GREY = "#C0392B", "#2c5aa0", "#767676"

cb = pd.read_csv(HERE / "matched_composition_X_combat.csv")
hm = pd.read_csv(HERE / "matched_composition_X_pca_harmony.csv")

fig = plt.figure(figsize=(7.2, 3.0))
gs = fig.add_gridspec(1, 2, wspace=0.34, width_ratios=[1.25, 1])

# A — ComBat: per integrated cluster male enrichment (female = negative)
axA = fig.add_subplot(gs[0, 0])
d = cb.copy()
lab = [f"{r.intg}\n({r.dominant_origin})" for r in d.itertuples()]
cols = [CLCOL[o] for o in d["dominant_origin"]]
axA.bar(range(len(d)), d["male_enrich"], color=cols, edgecolor="white", lw=0.5, width=0.7)
axA.axhline(0, color="#333", lw=0.7)
axA.set_xticks(range(len(d))); axA.set_xticklabels(lab, fontsize=6)
axA.set_ylabel("male enrichment\n(− = female-enriched)")
axA.set_title("A  ComBat integration recovers the sex shifts", loc="left", fontweight="bold", fontsize=8)
for i, r in enumerate(d.itertuples()):
    if r.dominant_origin == "C0":
        axA.annotate("MHC-II (C0)\nmale+", (i, r.male_enrich), xytext=(0, 6), textcoords="offset points",
                     ha="center", fontsize=5.4, color=CLCOL["C0"])
    if r.dominant_origin == "C1" and r.male_enrich < -0.1:
        axA.annotate("migratory (C1)\nfemale+", (i, r.male_enrich), xytext=(0, -16), textcoords="offset points",
                     ha="center", fontsize=5.4, color=CLCOL["C1"])

# B — Harmony vs ComBat for the key clusters
axB = fig.add_subplot(gs[0, 1])
def pick(df, origin, want_female=False):
    sub = df[df.dominant_origin == origin]
    return (sub["male_enrich"].min() if want_female else sub["male_enrich"].max())
items = [("C0 (MHC-II)\nmale enrich", pick(hm, "C0"), pick(cb, "C0")),
        ("C1 (migratory)\nfemale enrich", pick(hm, "C1", True), pick(cb, "C1", True))]
x = np.arange(len(items)); w = 0.36
axB.bar(x - w/2, [it[1] for it in items], w, color=GREY, edgecolor="white", lw=0.4, label="Harmony")
axB.bar(x + w/2, [it[2] for it in items], w, color=BLUE, edgecolor="white", lw=0.4, label="ComBat")
axB.axhline(0, color="#333", lw=0.7)
for i, it in enumerate(items):
    axB.text(i - w/2, it[1] + (0.005 if it[1] >= 0 else -0.012), f"{it[1]:+.2f}", ha="center", fontsize=5.2)
    axB.text(i + w/2, it[2] + (0.005 if it[2] >= 0 else -0.012), f"{it[2]:+.2f}", ha="center", fontsize=5.2)
axB.set_xticks(x); axB.set_xticklabels([it[0] for it in items], fontsize=6)
axB.set_ylabel("condition enrichment")
axB.set_title("B  ComBat recovers C0 male+; Harmony does not", loc="left", fontweight="bold", fontsize=8)
axB.legend(fontsize=5.6, loc="lower left")

for ext in ("svg", "pdf", "png"):
    fig.savefig(HERE / f"fig_integration_conclusions.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
print("wrote fig_integration_conclusions.{svg,pdf,png}")
