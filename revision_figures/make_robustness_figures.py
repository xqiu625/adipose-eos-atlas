#!/usr/bin/env python
"""Nature-style figures for the Reviewer #1 response (cluster robustness + identity).

Conclusion the figures defend: the six adipose eosinophil clusters are robust
eosinophil states, NOT misclassified macrophages. C0/C2 are not contamination
(MerTK/CD64/CD68 low, 0% doublets); their elevated C1q/Apoe is ambient RNA from
the macrophage-rich sample; C0 is defined by a genuine MHC-II program.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.spines.right": False, "axes.spines.top": False,
    "axes.linewidth": 0.8, "legend.frameon": False,
})
P = {"blue": "#0F4D92", "blue2": "#3775BA", "teal": "#42949E", "violet": "#9A4D8E",
     "gold": "#C9A227", "green": "#2E9E44", "red": "#C0392B", "grey": "#767676",
     "lgrey": "#CFCECE", "ink": "#2c2a22"}
CL = ["C0", "C1", "C2", "C3", "C4", "C5"]
CLCOL = {"C0": "#E41A1C", "C1": "#377EB8", "C2": "#4DAF4A", "C3": "#984EA3", "C4": "#FF7F00", "C5": "#A65628"}


def save(fig, name):
    for ext in ("svg", "pdf", "png"):
        fig.savefig(HERE / f"{name}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


# ── Load ──────────────────────────────────────────────────────────────────────
con = pd.read_csv(HERE / "contamination_by_cluster.csv").set_index("cluster").loc[CL]
con["contam"] = pd.read_csv(HERE / "ambient_correction_by_cluster.csv").set_index("cluster")["contam"].reindex(CL).values
ieg_frac = pd.read_csv(HERE / "c2_ieg_fraction_by_cluster.csv").set_index("cluster").loc[CL]
ieg_mod = pd.read_csv(HERE / "ieg_module_score_by_cluster.csv").set_index("cluster").loc[CL]
disp = pd.read_csv(HERE / "integration_dispersion_summary.csv")
strict = pd.read_csv(HERE / "siglecf_strict_macmarker_by_celltype.csv")
roles = pd.read_csv(HERE / "marker_roles_eos_vs_other.csv")
coord = pd.read_csv(HERE / "c0_coordination_summary.csv")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Identity: real eosinophils, not macrophage contamination
# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.42, height_ratios=[1, 1])

# A — per-cluster identity / contamination metrics (heatmap, each row own scale)
axA = fig.add_subplot(gs[0, :])
rows = [("Eosinophil score", con["eos_score"]),
        ("MHC-II score", con["mhc2_score"]),
        ("Mac ambient (C1q/Apoe/Lyz2)", con["mac_abundant"]),
        ("Mac intrinsic (MerTK/CD64/CD68)", con["mac_intrinsic"]),
        ("MerTK+ (%)", con["pct_Mertk"]),
        ("CD64+ (%)", con["pct_Fcgr1"]),
        ("Doublet (%)", con["pct_doublet"]),
        ("Contamination", con["contam"])]
M = np.vstack([r[1].values.astype(float) for r in rows])
Mz = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)
im = axA.imshow(Mz, cmap="RdBu_r", aspect="auto", vmin=-1.8, vmax=1.8)
axA.set_xticks(range(6)); axA.set_xticklabels(CL)
axA.set_yticks(range(len(rows))); axA.set_yticklabels([r[0] for r in rows], fontsize=6.3)
for i in range(len(rows)):
    for j in range(6):
        axA.text(j, i, f"{M[i, j]:g}", ha="center", va="center", fontsize=5,
                 color="white" if abs(Mz[i, j]) > 1.1 else "#222")
axA.set_title("A  Per-cluster identity & contamination metrics (row-scaled colour, raw values shown)",
              loc="left", fontweight="bold", fontsize=8)
cb = fig.colorbar(im, ax=axA, fraction=0.018, pad=0.01); cb.ax.tick_params(labelsize=5)
cb.set_label("row z-score", fontsize=5.5)

# B — strict macrophage co-expression by cell type
axB = fig.add_subplot(gs[1, 0])
st = strict.sort_values("pct_within_celltype", ascending=True)
cols = [P["red"] if c == "Eosinophils" else P["lgrey"] for c in st["celltype"]]
axB.barh(st["celltype"], st["pct_within_celltype"], color=cols, edgecolor="white", lw=0.4)
for y, v in enumerate(st["pct_within_celltype"]):
    axB.text(v + 0.05, y, f"{v}%", va="center", fontsize=5.4)
axB.set_xlabel("SiglecF+ AND (MerTK or CD64)+  (% of cell type)")
axB.set_title("B  Strict macrophage co-expression", loc="left", fontweight="bold", fontsize=8)
axB.set_xlim(0, max(st["pct_within_celltype"]) * 1.25)
axB.tick_params(axis="y", labelsize=6)

# C — marker roles: eos markers, discriminating, non-discriminating
axC = fig.add_subplot(gs[1, 1])
order = roles.copy()
rolecol = {"eosinophil": P["blue"], "mac_discriminating": P["red"], "mac_nondiscriminating": P["grey"]}
yb = np.arange(len(order)); w = 0.4
axC.barh(yb - w/2, order["pct_pos_eos"], w, color=[rolecol[r] for r in order["role"]], label="eosinophils", edgecolor="white", lw=0.3)
axC.barh(yb + w/2, order["pct_pos_nonEos"], w, color=[rolecol[r] for r in order["role"]], alpha=0.45, hatch="///", label="other cells", edgecolor="white", lw=0.3)
axC.set_yticks(yb); axC.set_yticklabels(order["gene"], fontsize=6)
axC.invert_yaxis()
axC.set_xlabel("% positive")
axC.set_title("C  Marker roles (eos vs other)", loc="left", fontweight="bold", fontsize=8)
from matplotlib.patches import Patch
axC.legend(handles=[Patch(facecolor=P["grey"], label="eosinophils"),
                    Patch(facecolor=P["grey"], alpha=0.45, hatch="///", label="other cells")],
           fontsize=5.2, loc="lower right")
axC.text(1.0, 1.02, "blue=eos marker · red=discriminating mac · grey=non-discriminating",
         transform=axC.transAxes, ha="right", fontsize=4.8, color=P["grey"])
save(fig, "fig_identity")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Robustness: clusters survive integration; C2 not an IEG artifact
# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(7.2, 3.2))
gs = fig.add_gridspec(1, 3, wspace=0.42)

# A — IEG module score + IEG marker fraction by cluster
axA = fig.add_subplot(gs[0, 0])
yb = np.arange(6); w = 0.4
axA.bar(yb - w/2, ieg_mod["mean"].loc[CL], w, color=[CLCOL[c] for c in CL], edgecolor="white", lw=0.3, label="IEG module score")
axA.bar(yb + w/2, ieg_frac["ieg_frac_top"].loc[CL], w, color=[CLCOL[c] for c in CL], alpha=0.45, hatch="///", edgecolor="white", lw=0.3, label="IEG frac of markers")
axA.set_xticks(yb); axA.set_xticklabels(CL, fontsize=6)
axA.set_title("A  C2 is not an IEG/stress artifact", loc="left", fontweight="bold", fontsize=8)
axA.set_ylabel("score / fraction")
axA.legend(handles=[Patch(facecolor=P["grey"], label="IEG module score"),
                    Patch(facecolor=P["grey"], alpha=0.45, hatch="///", label="IEG frac. of markers")],
           fontsize=5, loc="upper left")
axA.annotate("C2", (2, ieg_mod["mean"].loc["C2"]), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=6, color=CLCOL["C2"])

# B — integration purity by cluster x method
axB = fig.add_subplot(gs[0, 1])
methods = ["X_pca_harmony", "X_combat", "X_scVI"]
mlab = {"X_pca_harmony": "Harmony", "X_combat": "ComBat", "X_scVI": "scVI"}
mhatch = {"X_pca_harmony": "", "X_combat": "///", "X_scVI": ".."}
yb = np.arange(6); w = 0.26
for j, m in enumerate(methods):
    sub = disp[disp.embedding == m].set_index("orig_cluster").loc[CL]
    axB.bar(yb + (j-1)*w, sub["stays_together_purity"], w, color=[CLCOL[c] for c in CL],
            alpha=[1.0, 0.6, 0.8][j], hatch=mhatch[m], edgecolor="white", lw=0.3)
axB.axhline(0.5, color=P["lgrey"], ls="--", lw=0.7)
axB.set_xticks(yb); axB.set_xticklabels(CL, fontsize=6)
axB.set_ylabel("stays-together purity"); axB.set_ylim(0, 1)
axB.set_title("B  Cluster purity under integration", loc="left", fontweight="bold", fontsize=8)
axB.legend(handles=[Patch(facecolor=P["grey"], hatch=mhatch[m], alpha=[1.0,0.6,0.8][j], label=mlab[m]) for j, m in enumerate(methods)],
           fontsize=5, loc="lower right")

# C — C0 MHC-II coordination vs ambient
axC = fig.add_subplot(gs[0, 2])
lab = ["within\nMHC-II", "within\nmac-ambient", "MHC-II ×\nmac (cross)"]
vals = coord["mean_spearman"].values
cols = [P["blue"], P["grey"], P["red"]]
axC.bar(range(3), vals, color=cols, edgecolor="white", lw=0.4, width=0.62)
for i, v in enumerate(vals):
    axC.text(i, v + 0.002, f"{v:.3f}", ha="center", fontsize=6)
axC.set_xticks(range(3)); axC.set_xticklabels(lab, fontsize=5.6)
axC.set_ylabel("mean Spearman (C0 cells)")
axC.set_title("C  C0 MHC-II is a coordinated program", loc="left", fontweight="bold", fontsize=8)
axC.text(0.5, 0.93, "within-MHC-II > cross\n→ intrinsic, not ambient", transform=axC.transAxes,
         ha="center", va="top", fontsize=5, color=P["grey"])
save(fig, "fig_robustness")
print("wrote fig_identity.{svg,pdf,png} and fig_robustness.{svg,pdf,png}")
