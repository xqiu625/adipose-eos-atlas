
! source activate /bigdata/godziklab/shared/Xinru/velocyto_env

# Import necessary libraries

import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import scanpy as sc

# Define paths to .loom files for each sample

loom_files = {
    "RR0104_1": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0104_1/sample_alignments_3AOPN.loom",
    "RR0104_2": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0104_2/sample_alignments_8E9WA.loom",
    "RR0104_3": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0104_3/sample_alignments_8ITFG.loom",
    "RR0104_4": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0104_4/sample_alignments_B68J7.loom",
    "RR0508_5": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0508_5/sample_alignments_PHBRL.loom",
    "RR0508_6": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0508_6/sample_alignments_WXYY1.loom",
    "RR0508_7": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0508_7/sample_alignments_HCBDA.loom",
    "RR0508_8": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0508_8/sample_alignments_T9Y06.loom",
    "RR0912_9": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0912_9/sample_alignments_RDLTG.loom",
    "RR0912_10": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0912_10/sample_alignments_1BHVN.loom",
    "RR0912_11": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0912_11/sample_alignments_4KB06.loom",
    "RR0912_12": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR0912_12/sample_alignments_QMY1F.loom",
    "RR1316_13": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR1316_13/sample_alignments_4TBD1.loom",
    "RR1316_14": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR1316_14/sample_alignments_5XJ7Z.loom",
    "RR1316_15": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR1316_15/sample_alignments_KY4NG.loom",
    "RR1316_16": "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/RR1316_16/sample_alignments_W87X0.loom",
    
}

# Function to load and modify the CellID format of an AnnData object
def load_and_modif_loom(file_path):
    adata = anndata.read_loom(file_path)
    adata.obs.index = adata.obs.index.str.replace(":", "_").str.replace("x", "")
    return adata
    
# Load and modify loom file
anndata_objects = {name: load_and_modif_loom(path) for name, path in loom_files.items()}

# Define the mapping of AnnData object keys to their new prefixes
prefix_map = {
    "RR0104_1": "WT_Male_HFD_1",
    "RR0104_2": "WT_Female_HFD_1",
    "RR0104_3": "KO_Male_HFD_1",
    "RR0104_4": "KO_Female_HFD_1",
    "RR0508_5": "WT_Female_HFDInfected_1",
    "RR0508_6": "WT_Female_HFDInfected_2",
    "RR0508_7": "KO_Female_HFDInfected_1",
    "RR0508_8": "KO_Female_HFDInfected_2",
    "RR0912_9": "WT_Male_HFD_2",
    "RR0912_10": "WT_Male_HFD_3",
    "RR0912_11": "WT_Female_HFD_2",
    "RR0912_12": "WT_Female_HFD_3",
    "RR1316_13": "KO_Male_HFD_2",
    "RR1316_14": "KO_Male_HFD_3",
    "RR1316_15": "KO_Female_HFD_2",
    "RR1316_16": "KO_Female_HFD_3"
}

# Replace prefixes in cell IDs
for key, adata in anndata_objects.items():
    if key in prefix_map:
        new_prefix = prefix_map[key]
        adata.obs_names = adata.obs_names.map(lambda x: replace_prefix(x, new_prefix))
    print(f"Updated {key}:")
    print(adata.obs_names[:5])  # Print first 5 obs_names to verify
    print()

# Make var_names unique for each AnnData object
for adata in anndata_objects.values():
    adata.var_names_make_unique()
  
# Concatenate all AnnData objects
try:
    combined_adata = anndata.concat(
        anndata_objects.values(),
        join='outer',
        label='batch',
        keys=anndata_objects.keys()
    )
    print("Concatenation successful!")
    print(f"Combined AnnData shape: {combined_adata.shape}")
except Exception as e:
    print(f"Concatenation failed with error: {e}")


# Load additional data (UMAP coordinates, cell observations, clusters)
umap_cord = pd.read_csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/EOS_embeddings.csv")
sample_obs = pd.read_csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/EOS_cellID_obs.csv")
cell_clusters = pd.read_csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/EOS_clusters.csv")

print(umap_cord.head())

# Preprocess cell IDs in the loaded data
sample_obs['x'] = sample_obs['x'].str[:-2]
print(sample_obs['x'].head())

umap_cord['Unnamed: 0'] = umap_cord['Unnamed: 0'].str[:-2]
print(umap_cord['Unnamed: 0'].head())

cell_clusters['CellID'] = cell_clusters['CellID'].str[:-2]
print(cell_clusters['CellID'].head())

# filter for EOS cells
adata = combined_adata[np.isin(combined_adata.obs.index,sample_obs["x"])]

# Prepare and merge UMAP coordinates
adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {0:'CellID'})
print(adata_index.head()) 

umap_cord = umap_cord.rename(columns = {'Unnamed: 0':'CellID'})
print(umap_cord.head()) 

umap_ordered = adata_index.merge(umap_cord, on = "CellID")
print(umap_ordered.head()) 

umap_ordered = umap_ordered.iloc[:,1:]
adata.obsm['X_umap'] = umap_ordered.values


cell_clusters = cell_clusters.rename(columns = {'CellID':'CellID'})
clusters_ordered = adata_index.merge(cell_clusters, on = "CellID")
print(clusters_ordered.head()) 

# Prepare and merge cluster information
clusters_ordered = clusters_ordered.iloc[:,1:]
adata.obs['Clusters'] = clusters_ordered['SeuratClusters'].values
print(adata.obs[['Clusters']].head())

adata.obsm['X_cluster'] = clusters_ordered.values
print(adata.obsm['X_cluster'][:5])  # Print first 5 rows

# Save the processed AnnData object
adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata.h5ad")

# Perform scVelo preprocessing and analysis
scv.pp.filter_and_normalize(adata)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, basis = 'umap')

# Create and save UMAP plot with velocity streams
fig, ax = plt.subplots(figsize=(10, 8))
scv.pl.velocity_embedding_stream(adata, basis='umap', color='Clusters', 
                                 title='UMAP colored by Clusters', 
                                 ax=ax, show=False)

# Add labels to the clusters (if needed)
scv.pl.scatter(adata, basis='umap', color='Clusters', 
               legend_loc='on data', size=10, ax=ax, show=False)

# Save the plot
umap_output_file = "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_clusters_scVelo_stochastic.png"
plt.savefig(umap_output_file, dpi=300, bbox_inches='tight')
plt.close()



