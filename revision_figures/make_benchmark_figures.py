#!/usr/bin/env python
"""Nature-style figures for the eosinophil integration-method benchmark.

Conclusion the figure defends: across seven integration tools, no UNSUPERVISED
method improves biological conservation over the un-integrated representation;
integration only trades bio-conservation for batch mixing. scANVI's apparent lead
is supervised on C0-C5 (circular) and is flagged as such.
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
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

P = {"blue_main": "#0F4D92", "blue_secondary": "#3775BA", "teal": "#42949E",
     "violet": "#9A4D8E", "gold": "#C9A227", "neutral_dark": "#4D4D4D",
     "neutral_mid": "#767676", "neutral_light": "#CFCECE",
     "delta_up": "#2E9E44", "delta_down": "#E53935"}

NAME = {"X_pca": "Unintegrated", "X_combat": "ComBat", "X_pca_harmony": "Harmony",
        "X_scanorama": "Scanorama", "X_scVI": "scVI", "X_scANVI": "scANVI"}
COLOR = {"Unintegrated": P["neutral_dark"], "ComBat": P["blue_secondary"],
         "Harmony": P["blue_main"], "Scanorama": P["teal"], "scVI": P["violet"],
         "scANVI": P["gold"]}
ORDER = ["Unintegrated", "ComBat", "Harmony", "scVI", "scANVI"]

# ── Load ──────────────────────────────────────────────────────────────────────
m = pd.read_csv(HERE / "integration_benchmark_metrics.csv", index_col=0)
mtype = m.loc["Metric Type"]
m = m.drop(index="Metric Type").astype(float)
m.index = [NAME[i] for i in m.index]
m = m.loc[ORDER]
rec = pd.read_csv(HERE / "integration_cluster_recovery.csv")
rec["disp"] = rec["method"].map({"unintegrated": "Unintegrated", "combat": "ComBat",
    "harmony": "Harmony", "scanorama": "Scanorama", "scvi": "scVI", "scanvi": "scANVI"})
rec = rec.set_index("disp").loc[ORDER]

bio_base = m.loc["Unintegrated", "Bio conservation"]

# ── Figure: 2x2 quantitative grid ─────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
fig.subplots_adjust(left=0.085, right=0.985, top=0.93, bottom=0.085,
                    wspace=0.32, hspace=0.42)

# Panel A — batch vs bio tradeoff scatter
ax = axes[0, 0]
ax.axhline(bio_base, color=P["neutral_light"], lw=0.8, ls="--", zorder=0)
ax.text(0.985, bio_base + 0.004, "un-integrated bio level", ha="right", va="bottom",
        fontsize=5.6, color=P["neutral_mid"], transform=ax.get_yaxis_transform())
for met in ORDER:
    x, y = m.loc[met, "Batch correction"], m.loc[met, "Bio conservation"]
    sup = met == "scANVI"
    ax.scatter(x, y, s=46, color=COLOR[met], edgecolor="white", lw=0.6,
               marker="D" if sup else "o", zorder=3)
    dy = 0.012 if met not in ("scVI", "Harmony") else -0.02
    ax.annotate(met + ("*" if sup else ""), (x, y), xytext=(0, 9 if dy > 0 else -11),
                textcoords="offset points", ha="center", fontsize=6,
                color=COLOR[met] if met != "scANVI" else "#8a6d00")
ax.set_xlabel("Batch correction  (mixing →)")
ax.set_ylabel("Bio conservation  (biology →)")
ax.set_title("A  Batch–bio trade-off", loc="left", fontweight="bold", fontsize=8)
ax.set_xlim(0.55, 0.88); ax.set_ylim(0.37, 0.70)
ax.annotate("supervised on C0–C5\n(not unsupervised-comparable)", (m.loc["scANVI", "Batch correction"], m.loc["scANVI", "Bio conservation"]),
            xytext=(-4, -26), textcoords="offset points", ha="center", fontsize=5.2,
            color="#8a6d00")

# Panel B — aggregate scores grouped bars
ax = axes[0, 1]
sub = ["Batch correction", "Bio conservation", "Total"]
xb = np.arange(len(ORDER)); w = 0.26
hatch = {"Batch correction": "", "Bio conservation": "////", "Total": ".."}
shade = {"Batch correction": 1.0, "Bio conservation": 0.55, "Total": 0.8}
for j, s in enumerate(sub):
    for i, met in enumerate(ORDER):
        c = COLOR[met]
        ax.bar(xb[i] + (j - 1) * w, m.loc[met, s], w, color=c, alpha=shade[s],
               edgecolor="white", lw=0.4, hatch=hatch[s],
               label=s if i == 0 else None)
ax.set_xticks(xb); ax.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=6)
ax.set_ylabel("scIB score"); ax.set_ylim(0, 1)
ax.set_title("B  Aggregate scores", loc="left", fontweight="bold", fontsize=8)
from matplotlib.patches import Patch
leg = [Patch(facecolor=P["neutral_mid"], hatch=hatch[s], label=s) for s in sub]
ax.legend(handles=leg, fontsize=5.4, loc="upper right", ncol=1, handlelength=1.3)

# Panel C — sub-metric heatmap
ax = axes[1, 0]
cols = list(m.columns[:-3])  # drop aggregate cols
H = m[cols].T
im = ax.imshow(H.values, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=6)
ax.set_yticks(range(len(cols)))
ax.set_yticklabels([c.replace(" comparison", "") for c in cols], fontsize=5.6)
for yi in range(len(cols)):
    for xi in range(len(ORDER)):
        v = H.values[yi, xi]
        ax.text(xi, yi, f"{v:.2f}", ha="center", va="center", fontsize=4.6,
                color="white" if v > 0.6 else "#222")
nbio = (mtype[cols] == "Bio conservation").sum()
ax.axhline(nbio - 0.5, color="white", lw=1.6)
ax.text(-1.65, (nbio - 1) / 2, "bio", rotation=90, va="center", fontsize=5.6, color=P["neutral_dark"])
ax.text(-1.65, nbio + (len(cols) - nbio - 1) / 2, "batch", rotation=90, va="center", fontsize=5.6, color=P["neutral_dark"])
ax.set_title("C  scIB sub-metrics", loc="left", fontweight="bold", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=5)

# Panel D — cluster recovery
ax = axes[1, 1]
xb = np.arange(len(ORDER)); w = 0.38
ax.bar(xb - w / 2, rec["leiden_ARI_vs_C0C5"], w, color=[COLOR[m_] for m_ in ORDER],
       edgecolor="white", lw=0.4, label="ARI")
ax.bar(xb + w / 2, rec["leiden_NMI_vs_C0C5"], w, color=[COLOR[m_] for m_ in ORDER],
       alpha=0.5, edgecolor="white", lw=0.4, hatch="////", label="NMI")
for i, met in enumerate(ORDER):
    ax.text(xb[i], max(rec.loc[met, "leiden_ARI_vs_C0C5"], rec.loc[met, "leiden_NMI_vs_C0C5"]) + 0.015,
            f"{int(rec.loc[met,'n_leiden_clusters'])}c", ha="center", fontsize=4.8, color=P["neutral_mid"])
ax.set_xticks(xb); ax.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=6)
ax.set_ylabel("Leiden agreement vs C0–C5"); ax.set_ylim(0, 0.72)
ax.set_title("D  Cluster recovery", loc="left", fontweight="bold", fontsize=8)
leg = [Patch(facecolor=P["neutral_mid"], label="ARI"),
       Patch(facecolor=P["neutral_mid"], alpha=0.5, hatch="////", label="NMI")]
ax.legend(handles=leg, fontsize=5.4, loc="upper left", handlelength=1.3)
ax.text(0.5, 0.72, "*scANVI supervised on C0–C5; '#c' = clusters found (res 0.4)",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=4.8, color=P["neutral_mid"])

for ext in ("svg", "pdf", "png"):
    fig.savefig(HERE / f"fig_integration_benchmark.{ext}",
                dpi=600 if ext == "png" else None, bbox_inches="tight")
print("wrote fig_integration_benchmark.{svg,pdf,png}")
