#!/usr/bin/env python3
"""
16c_per_cluster_heatmap.py
Per-cluster pySCENIC combined heatmap (adapted from 1208_Combined_Eos_Mac_Heatmap.py)

Reads per-cluster AUCell matrices from scenic_per_cluster/C0..C5/auc_matrix.csv
and mac_aucell_scored.csv, finds Eos-unique regulons, plots combined heatmap.

SLURM: --partition=batch --mem=32G --cpus-per-task=4 --time=1:00:00
Conda: /bigdata/godziklab/shared/Xinru/cellrank_env
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "savefig.dpi": 300})

ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"])

BASE          = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR    = os.path.join(BASE, "scenic_per_cluster")
MAC_SCORED    = os.path.join(BASE, "scenic", "mac_aucell_scored.csv")
MAC_META      = os.path.join(BASE, "mac_subtyped_metadata.csv")
FIG_DIR       = os.path.join(BASE, "scenic", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

EOS_CLUSTERS = ["C0", "C1", "C2", "C3", "C4", "C5"]
MAC_SUBTYPES = ["M1", "M2_Functional", "M2_Surface"]
MAC_RENAME   = {"Mac_M1": "M1", "Mac_M2_Func": "M2_Functional",
                "Mac_M2_Surf": "M2_Surface"}

CELLS_PER_CLUSTER = 200
RANDOM_SEED       = 42
ZSCORE_CAP        = 2.5
TOP_PER_CLUSTER   = 5

EXCLUDE_REGULONS = {
    "Erf", "Pura", "Zbtb14", "Rfx5", "Smad3", "Jund", "Tef", "Bclaf1",
    "Crem", "Atf4", "Etv5", "Zbtb7b", "Sp3", "Elk1", "Cebpb", "Max",
    "Atf3", "Ikzf1", "Atf1", "Cux1", "Stat6", "Elk3"
}

print("=" * 60)
print("16c_per_cluster_heatmap.py")
print("=" * 60)

# =============================================================================
# Step 1: Load per-cluster AUCell matrices
# =============================================================================
print("\nStep 1: Loading per-cluster AUCell matrices...")
np.random.seed(RANDOM_SEED)

eos_data = {}
gene_counts = {}   # tf -> n_targets (from regulon names like "Mybl2(10)")

for cl in EOS_CLUSTERS:
    auc_path = os.path.join(SCENIC_DIR, cl, "auc_matrix.csv")
    if os.path.exists(auc_path):
        df = pd.read_csv(auc_path, index_col=0)
        eos_data[cl] = df
        # Parse target gene counts from column names e.g. "Mybl2(10)"
        for col in df.columns:
            if "(" in col:
                tf   = col.split("(")[0].strip()
                try:
                    n = int(col.split("(")[1].rstrip(")").strip())
                    if tf not in gene_counts or n > gene_counts[tf]:
                        gene_counts[tf] = n
                except Exception:
                    pass
        print(f"  {cl}: {df.shape[0]} cells x {df.shape[1]} regulons")
    else:
        print(f"  WARNING: {auc_path} not found")

if not eos_data:
    raise FileNotFoundError(
        f"No per-cluster AUCell files found in {SCENIC_DIR}. "
        "Run 14c_export_per_cluster_loom.R + 15c_scenic_per_cluster.sh first.")

# =============================================================================
# Step 2: Load Mac AUCell scores
# =============================================================================
print("\nStep 2: Loading Mac scores...")
mac_mean_df = None

if os.path.exists(MAC_SCORED):
    mac_scored = pd.read_csv(MAC_SCORED, index_col=0)
    if os.path.exists(MAC_META):
        mac_meta   = pd.read_csv(MAC_META, index_col=0)
        sub_col    = "mac_subtype" if "mac_subtype" in mac_meta.columns \
                     else mac_meta.columns[0]
        common     = mac_meta.index.intersection(mac_scored.index)
        mac_scored = mac_scored.loc[common]
        mac_scored["subtype"] = mac_meta.loc[common, sub_col].map(
            lambda x: MAC_RENAME.get(x, x)).values
    elif "mac_subtype" in mac_scored.columns:
        mac_scored["subtype"] = mac_scored["mac_subtype"].map(
            lambda x: MAC_RENAME.get(x, x))
    else:
        print("  WARNING: mac subtype info not found")
        mac_scored["subtype"] = "Mac"

    reg_cols = [c for c in mac_scored.columns if c != "subtype"]
    mac_mean_df = mac_scored.groupby("subtype")[reg_cols].mean()
    print(f"  Mac subtypes: {list(mac_mean_df.index)}")
else:
    print(f"  Mac scores not found at {MAC_SCORED}")

# =============================================================================
# Step 3: Compute per-cluster mean AUC (subsampled, TF-level)
# =============================================================================
print("\nStep 3: Computing per-cluster mean AUC...")

def extract_tf(col_name):
    return col_name.split("(")[0].strip()

# Build unified TF set
all_tfs = set()
eos_tf_map = {}   # cl -> {tf: col_name}
for cl, df in eos_data.items():
    eos_tf_map[cl] = {}
    for col in df.columns:
        tf = extract_tf(col)
        eos_tf_map[cl][tf] = col
        all_tfs.add(tf)

all_tfs = sorted(all_tfs)

# Subsample and compute means
eos_means = {}
for cl in EOS_CLUSTERS:
    if cl not in eos_data:
        continue
    df  = eos_data[cl]
    n   = len(df)
    idx = np.random.choice(n, min(n, CELLS_PER_CLUSTER), replace=(n < CELLS_PER_CLUSTER))
    sub = df.iloc[idx]
    row = {}
    for tf in all_tfs:
        col = eos_tf_map[cl].get(tf)
        row[tf] = sub[col].mean() if col else 0.0
    eos_means[cl] = row
    print(f"  {cl}: {n} cells, {sum(1 for v in row.values() if v > 0)} active regulons")

eos_mean_df = pd.DataFrame(eos_means)  # tfs x clusters

# Mac means aligned to same TF names
mac_tf_means = {}
if mac_mean_df is not None:
    for subtype in mac_mean_df.index:
        row = {}
        for tf in all_tfs:
            # match tf to mac column (may have (n) suffix or (+) suffix)
            matched = [c for c in mac_mean_df.columns
                       if extract_tf(c.replace("(+)","").replace("(-)", "")) == tf]
            row[tf] = mac_mean_df.loc[subtype, matched[0]] if matched else 0.0
        mac_tf_means[subtype] = row
    mac_tf_df = pd.DataFrame(mac_tf_means)   # tfs x mac subtypes
else:
    mac_tf_df = None

# =============================================================================
# Step 4: Find Eos-unique regulons
# =============================================================================
print("\nStep 4: Finding Eos-unique regulons...")

results = []
for tf in all_tfs:
    if tf in EXCLUDE_REGULONS:
        continue

    eos_vals = {cl: eos_mean_df.loc[tf, cl]
                for cl in EOS_CLUSTERS if cl in eos_mean_df.columns}
    if not eos_vals:
        continue

    max_cl  = max(eos_vals, key=eos_vals.get)
    max_auc = eos_vals[max_cl]
    if max_auc < 0.05:
        continue

    mac_mean = (mac_tf_df.loc[tf].mean()
                if mac_tf_df is not None and tf in mac_tf_df.index else 0.0)

    fold_change = max_auc / mac_mean if mac_mean > 0 else float("inf")
    sorted_eos  = sorted(eos_vals.values(), reverse=True)
    uniqueness  = sorted_eos[0] - sorted_eos[1] if len(sorted_eos) >= 2 else sorted_eos[0]

    results.append({
        "tf": tf, "max_cluster": max_cl,
        "max_auc": max_auc, "fold_change": fold_change, "uniqueness": uniqueness
    })

res_df = pd.DataFrame(results)
print(f"  Evaluated: {len(res_df)} TFs")

# Top N per cluster by uniqueness
top_rows = []
for cl in EOS_CLUSTERS:
    sub = res_df[(res_df["max_cluster"] == cl) &
                 (res_df["fold_change"] >= 1.2) &
                 (res_df["uniqueness"] >= 0.001)].nlargest(TOP_PER_CLUSTER, "uniqueness")
    top_rows.append(sub)

top_df = pd.concat(top_rows, ignore_index=True).drop_duplicates("tf")
print(f"  Selected {len(top_df)} Eos-unique regulons")
for cl in EOS_CLUSTERS:
    sub = top_df[top_df["max_cluster"] == cl]
    if len(sub):
        print(f"    {cl}: {', '.join(sub['tf'].tolist())}")

# =============================================================================
# Step 5: Build heatmap matrix
# =============================================================================
print("\nStep 5: Building heatmap...")

selected_tfs = top_df["tf"].tolist()

# Eos means for selected TFs
eos_heat = eos_mean_df.loc[selected_tfs,
                            [c for c in EOS_CLUSTERS if c in eos_mean_df.columns]].T

# Mac means for selected TFs
if mac_tf_df is not None:
    mac_heat = mac_tf_df.loc[
        selected_tfs, [s for s in MAC_SUBTYPES if s in mac_tf_df.columns]
    ].T
    mac_heat.index = [{"M1": "M1", "M2_Functional": "M2\nFunct.",
                       "M2_Surface": "M2\nSurf."}.get(s, s) for s in mac_heat.index]
    combined_raw = pd.concat([eos_heat, mac_heat])
else:
    combined_raw = eos_heat

# Row-wise z-score (per regulon, across all cell types)
mat      = combined_raw.T   # tfs x cell_types
row_mean = mat.mean(axis=1)
row_std  = mat.std(axis=1).replace(0, 1)
combined = ((mat.T - row_mean) / row_std).T.clip(-ZSCORE_CAP, ZSCORE_CAP).T
# combined: cell_types x tfs

plot_data = combined.T  # tfs x cell_types

# Y-axis labels with target gene count
y_labels = []
for tf in selected_tfs:
    n = gene_counts.get(tf, "")
    y_labels.append(f"{tf}({n})" if n != "" else tf)

n_rows = len(y_labels)
n_cols = len(combined.index)
n_eos  = len([c for c in EOS_CLUSTERS if c in eos_heat.index])

cell_h   = 0.38
cell_w   = 0.55
left_pad = 1.8
right_pad= 1.2
top_pad  = 0.7
bot_pad  = 1.0

fig_h = n_rows * cell_h + top_pad + bot_pad
fig_w = n_cols * cell_w + left_pad + right_pad

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

ax_left   = left_pad  / fig_w
ax_bottom = bot_pad   / fig_h
ax_width  = (n_cols * cell_w) / fig_w
ax_height = (n_rows * cell_h) / fig_h
ax.set_position([ax_left, ax_bottom, ax_width, ax_height])

im = ax.imshow(plot_data.values, cmap=ZSCORE_CMAP, aspect="auto",
               vmin=-ZSCORE_CAP, vmax=ZSCORE_CAP)

# X-axis labels
x_tick_labels = []
for lbl in list(combined.index):
    if "Funct" in str(lbl):   x_tick_labels.append("M2")
    elif "Surf" in str(lbl):  x_tick_labels.append("M2")
    else:                      x_tick_labels.append(str(lbl))

ax.set_xticks(range(n_cols))
ax.set_xticklabels(x_tick_labels, fontsize=10, fontweight="bold")

# Y-axis on LEFT
ax.set_yticks(range(n_rows))
ax.set_yticklabels(y_labels, fontsize=9, ha="right")
ax.yaxis.set_ticks_position("left")
ax.tick_params(axis="both", length=0)

# Vertical separator between Eos and Mac
if mac_tf_df is not None and n_eos < n_cols:
    ax.axvline(x=n_eos - 0.5, color="black", linewidth=2.0)

# Group labels at bottom
def data_x_to_fig(x):
    return ax_left + (x + 0.5) / n_cols * ax_width

line_y   = ax_bottom - 0.32 / fig_h
label_y  = ax_bottom - 0.60 / fig_h
sublbl_y = ax_bottom - 0.20 / fig_h

# Eosinophils
fig.add_artist(plt.Line2D(
    [data_x_to_fig(-0.5), data_x_to_fig(n_eos - 0.5)],
    [line_y, line_y], transform=fig.transFigure,
    color="black", linewidth=1.2, clip_on=False))
fig.text(data_x_to_fig((n_eos - 1) / 2), label_y,
         "Eosinophils", ha="center", va="top", fontsize=10,
         transform=fig.transFigure)

if mac_tf_df is not None and n_eos < n_cols:
    n_mac = n_cols - n_eos
    fig.add_artist(plt.Line2D(
        [data_x_to_fig(n_eos - 0.5), data_x_to_fig(n_cols - 0.5)],
        [line_y, line_y], transform=fig.transFigure,
        color="black", linewidth=1.2, clip_on=False))
    fig.text(data_x_to_fig(n_eos + (n_mac - 1) / 2), label_y,
             "Mac", ha="center", va="top", fontsize=10,
             transform=fig.transFigure)
    if n_mac >= 3:
        fig.text(data_x_to_fig(n_eos + 1), sublbl_y,
                 "Funct.", ha="center", va="top", fontsize=8,
                 transform=fig.transFigure)
        fig.text(data_x_to_fig(n_eos + 2), sublbl_y,
                 "Surf.", ha="center", va="top", fontsize=8,
                 transform=fig.transFigure)

ax.set_title("Eos Subcluster Unique Regulons", fontsize=13,
             fontweight="bold", pad=10)

# Colorbar
cbar_ax = fig.add_axes([ax_left + ax_width + 0.03 / fig_w * fig_h,
                         ax_bottom + ax_height * 0.2,
                         0.03 / fig_w * fig_h,
                         ax_height * 0.6])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("Z-score", fontsize=9)
cbar.set_ticks([-2, -1, 0, 1, 2])
cbar.ax.tick_params(labelsize=8)

out_pdf = os.path.join(FIG_DIR, "Fig5_regulon_heatmap_per_cluster.pdf")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {out_pdf}")

top_df.to_csv(os.path.join(BASE, "scenic", "eos_unique_regulons_per_cluster.csv"),
              index=False)
print("\n16c complete.")
