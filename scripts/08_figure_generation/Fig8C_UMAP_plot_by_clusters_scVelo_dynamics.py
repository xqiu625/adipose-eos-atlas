! source activate /bigdata/godziklab/shared/Xinru/velocyto_env

import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

# Set global font to Arial
plt.rcParams['font.family'] = 'Arial'

# Load the AnnData object
adata = scv.read("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_dynamics.h5ad")

# Perform velocity analysis
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# Create and save UMAP plot with velocity streams
fig, ax = plt.subplots(figsize=(14, 12))  # Increased figure size

# Define a custom colormap for clusters
cluster_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

scv.pl.velocity_embedding_stream(adata, basis='umap', color='Clusters', 
                                 title='', 
                                 ax=ax, show=False,
                                 arrow_size=2,
                                 arrow_length=2,
                                 density=1.5,
                                 color_map=cluster_colors,
                                 )

# Add labels to the clusters
for cluster in adata.obs['Clusters'].cat.categories:
    cells = adata[adata.obs['Clusters'] == cluster]
    centroid = np.median(cells.obsm['X_umap'], axis=0)
    ax.text(centroid[0], centroid[1], cluster, 
            fontsize=20, fontweight='bold', 
            ha='center', va='center',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.5))

# Customize the plot further
ax.set_xlabel('UMAP1', fontsize=24)
ax.set_ylabel('UMAP2', fontsize=24)
ax.tick_params(axis='both', which='major', labelsize=18)

# Remove the legend as cluster labels are now on the plot
ax.legend().remove()

# Adjust layout to prevent cutting off labels
plt.tight_layout()

# Save the plot as PDF
umap_output_file = "/bigdata/nairlab/shared/XQiu/EOS_project/figures/grant_figures/UMAP_plot_by_clusters_scVelo_dynamics.pdf"
plt.savefig(umap_output_file, format='pdf', dpi=300, bbox_inches='tight')
plt.close()
