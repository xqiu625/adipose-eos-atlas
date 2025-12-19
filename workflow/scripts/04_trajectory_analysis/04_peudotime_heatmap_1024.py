#!/usr/bin/env python3
"""
Generate scVelo heatmap and cluster legend in multiple formats
Saves in PDF, SVG, and EPS formats
"""

import matplotlib
matplotlib.use('Agg')  # Set backend before importing pyplot
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
import os
from matplotlib.patches import Patch

print("\n" + "="*70)
print("SCVELO HEATMAP AND LEGEND GENERATOR")
print("="*70)

# Define cluster colors
CLUSTER_COLORS = {
    "C0": "#6CA6CD",
    "C1": "#E84A5F",
    "C2": "#FF6F1F",
    "C3": "#CBAACB",
    "C4": "#B03A48",
    "C5": "#5BA4A4",
    "C6": "#2F3C7E",
    "Gurtner_Normal_Blood": "#FFD700",
    "Gurtner_Normal_Bone Marrow": "#32CD32"
}

# Set output directory
output_dir = "/bigdata/godziklab/shared/xinruqiu/302025/R16_20251002/CellRank_analysis/violin_plots/"
scv.settings.figdir = output_dir
os.makedirs(output_dir, exist_ok=True)

# =============================================================================
# PART 1: Generate Heatmap with Legend
# =============================================================================
print("\n[1/2] Generating heatmap with legend...")

# Load data
print("   Loading data...")
adata = sc.read('/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_this_study_gurtner_harmonized_by_organ_cellrank_pseudotime_analyzed.h5ad')

# Filter adata to keep only clusters with assigned colors
wanted_clusters = list(CLUSTER_COLORS.keys())
adata = adata[adata.obs['cluster'].isin(wanted_clusters)].copy()

print(f"   Filtered to {adata.n_obs} cells across {len(wanted_clusters)} clusters")
print(f"   Clusters kept: {sorted(adata.obs['cluster'].unique())}")

# Get top 1000 highly variable genes
print("   Selecting top 1000 HVGs...")
top_hvg = adata.var[adata.var['highly_variable']].sort_values('dispersions_norm', ascending=False).index[:1000]

# Generate heatmap
print("   Generating heatmap...")
fig = scv.pl.heatmap(adata, 
                     var_names=top_hvg, 
                     sortby='dpt_pseudotime', 
                     col_color='cluster', 
                     n_convolve=100,
                     show=False)

# Create legend patches for each cluster
unique_clusters = adata.obs['cluster'].unique()
legend_patches = [Patch(facecolor=CLUSTER_COLORS[cluster], 
                       label=cluster) 
                 for cluster in sorted(unique_clusters)]

# Add legend to the figure
plt.legend(handles=legend_patches, 
          loc='center left', 
          bbox_to_anchor=(1.05, 0.5),
          frameon=False,
          title='Cluster')

plt.tight_layout()

# Save in multiple formats
print("   Saving heatmap in multiple formats...")
heatmap_base = os.path.join(output_dir, 'heatmap_hvg_1000_dpt_with_legend')
for ext in ['pdf', 'svg', 'eps']:
    filepath = f'{heatmap_base}.{ext}'
    plt.savefig(filepath, bbox_inches='tight', dpi=300)
    if os.path.exists(filepath):
        file_size = os.path.getsize(filepath) / (1024 * 1024)  # MB
        print(f"      ✓ Saved: heatmap_hvg_1000_dpt_with_legend.{ext} ({file_size:.1f} MB)")
    else:
        print(f"      ✗ Failed: heatmap_hvg_1000_dpt_with_legend.{ext}")

plt.close()

# =============================================================================
# PART 2: Generate Standalone Legend
# =============================================================================
print("\n[2/2] Generating standalone cluster legend...")

# Create a figure with just the legend
fig, ax = plt.subplots(figsize=(4, 6))
ax.axis('off')  # Hide axes

# Create legend patches
legend_patches = [Patch(facecolor=color, label=cluster) 
                  for cluster, color in CLUSTER_COLORS.items()]

# Add legend
legend = ax.legend(handles=legend_patches, 
                   loc='center',
                   frameon=True,
                   title='Cluster',
                   fontsize=12,
                   title_fontsize=14)

plt.tight_layout()

# Save in multiple formats
print("   Saving legend in multiple formats...")
legend_base = os.path.join(output_dir, 'cluster_legend')
for ext in ['pdf', 'svg', 'eps']:
    filepath = f'{legend_base}.{ext}'
    plt.savefig(filepath, bbox_inches='tight', dpi=300)
    if os.path.exists(filepath):
        file_size = os.path.getsize(filepath) / 1024  # KB
        print(f"      ✓ Saved: cluster_legend.{ext} ({file_size:.1f} KB)")
    else:
        print(f"      ✗ Failed: cluster_legend.{ext}")

plt.close()

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "="*70)
print("GENERATION COMPLETE!")
print("="*70)

print(f"\nAll files saved to:")
print(f"   {output_dir}")

print("\nGenerated files:")
print("   Heatmap with legend:")
print("      • heatmap_hvg_1000_dpt_with_legend.pdf")
print("      • heatmap_hvg_1000_dpt_with_legend.svg")
print("      • heatmap_hvg_1000_dpt_with_legend.eps")
print("\n   Standalone legend:")
print("      • cluster_legend.pdf")
print("      • cluster_legend.svg")
print("      • cluster_legend.eps")

print("\n" + "="*70 + "\n")