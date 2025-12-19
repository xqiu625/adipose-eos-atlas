import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

# Load the merged AnnData object
adata = scv.read("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata.h5ad")

# Investigate AnnData metadata
print("AnnData object loaded.")
print(f"Shape: {adata.shape}")
print(adata.obs.head())

print(adata.obs['condition'].value_counts())

print("\nSex composition:")
print(adata.obs['sex'].value_counts())

print("\nGenotype composition:")
print(adata.obs['genotype'].value_counts())

print("\nInfection status composition:")
print(adata.obs['infection'].value_counts())



# Add metadata columns to adata.obs based on CellID

def parse_cellid(cellid):
    """
    Parse the CellID to extract condition, sex, genotype, and infection status.
    Assumes CellID format: <Condition>_<Sex>_<HFD or HFDInfected>_<SampleNum>_<Barcode>
    Example: WT_Male_HFD_1_AAAGGCAAGACTAACAACTTTAGG
    """
    # Split by underscore, but only first 4 are metadata, rest is barcode
    parts = cellid.split('_')
    if len(parts) < 5:
        # fallback: return Nones
        return pd.Series({
            "condition": None,
            "sex": None,
            "genotype": None,
            "infection": None
        })
    genotype = parts[0]
    sex = parts[1]
    hfd_part = parts[2]
    # Determine infection status
    if "Infected" in hfd_part:
        infection = "Infected"
        condition = f"{genotype}_{sex}_HFDInfected"
    else:
        infection = "Uninfected"
        condition = f"{genotype}_{sex}_HFD"
    return pd.Series({
        "condition": condition,
        "sex": sex,
        "genotype": genotype,
        "infection": infection
    })

# Apply the parser to the index (CellID)
meta_df = adata.obs.index.to_series().apply(parse_cellid)

# Add new columns to adata.obs
for col in meta_df.columns:
    adata.obs[col] = meta_df[col].values

adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_250923.h5ad")



# Recover dynamics for the dynamical model
scv.tl.recover_dynamics(adata)

# Save the AnnData object to a file
adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_dynamics.h5ad")

scv.tl.velocity(adata, mode = "dynamical")
scv.tl.velocity_graph(adata)

# Create and save UMAP plot with velocity streams
fig, ax = plt.subplots(figsize=(10, 8))
scv.pl.velocity_embedding_stream(adata, basis='umap', color='Clusters', 
                                 title='UMAP colored by Clusters', 
                                 ax=ax, show=False)

# Add labels to the clusters (if needed)
scv.pl.scatter(adata, basis='umap', color='Clusters', 
               legend_loc='on data', size=10, ax=ax, show=False)

# Save the plot
umap_output_file = "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_clusters_scVelo_dynamics.pdf"
plt.savefig(umap_output_file, dpi=300, bbox_inches='tight', format='pdf')
plt.close()



adata = scv.read("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_dynamics.h5ad")
# Apply EOS color palette
eos_palette = {
    # Direct cluster names (C0, C1, etc.)
    "C0": "#6CA6CD",
    "C1": "#E84A5F", 
    "C2": "#FF6F1F",
    "C3": "#CBAACB",
    "C4": "#B03A48",
    "C5": "#5BA4A4",
    "C6": "#2F3C7E",
    # Also keep Eos_ prefixed versions as backup
    "Eos_C0": "#6CA6CD",
    "Eos_C1": "#E84A5F", 
    "Eos_C2": "#FF6F1F",
    "Eos_C3": "#CBAACB",
    "Eos_C4": "#B03A48",
    "Eos_C5": "#5BA4A4",
    "Eos_C6": "#2F3C7E"
}

# Ensure Clusters column is categorical
if not pd.api.types.is_categorical_dtype(adata.obs['Clusters']):
    adata.obs['Clusters'] = adata.obs['Clusters'].astype('category')

# Set colors for clusters
categories = adata.obs['Clusters'].cat.categories
colors = []
for cat in categories:
    if str(cat) in eos_palette:
        colors.append(eos_palette[str(cat)])
    elif f"Eos_{cat}" in eos_palette:
        colors.append(eos_palette[f"Eos_{cat}"])
    else:
        # Fallback color for undefined clusters
        colors.append("#808080")

# Apply colors to adata
adata.uns['Clusters_colors'] = colors

# Run velocity analysis
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# Create and save UMAP plot with velocity streams
fig, ax = plt.subplots(figsize=(10, 8))
scv.pl.velocity_embedding_stream(adata, basis='umap', color='Clusters',
                                  title='UMAP colored by Clusters',
                                  ax=ax, show=False)

# Add labels to the clusters (if needed)
scv.pl.scatter(adata, basis='umap', color='Clusters',
                legend_loc='on data', size=10, ax=ax, show=False)

# Save the plot
umap_output_file = "/bigdata/nairlab/shared/XQiu/EOS_project/R08_20250716/01_EOS_Trajectory/EOS_scVelo_dynamics.png"
plt.savefig(umap_output_file, dpi=300, bbox_inches='tight', format='png')
plt.close()

print(f"✅ UMAP plot saved to: {umap_output_file}")
print(f"🎨 Colors applied: {len(colors)} clusters")
print(f"📊 Cluster categories: {list(categories)}")