import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors

# Load the mouse_ref sample dataset
mouse_ref = sc.read_h5ad('/home/xinruqiu/UCE/data/mouse_kidney_norn.h5ad')
print(mouse_ref)
print(mouse_ref.obs_keys())  # Check available metadata
print(mouse_ref.var_names[:10])  # Check some gene names
print(mouse_ref.obs['metacell'].value_counts())
print(mouse_ref.obs['group'].value_counts())
eosinophil_mask = mouse_ref.obs['group'] == 'eosinophil'
eosinophil_count = eosinophil_mask.sum()
print(f"Number of eosinophil cells: {eosinophil_count}")
print(mouse_ref.obs['group'].unique())
unique_groups = mouse_ref.obs['group'].unique()
print("Unique groups in mouse_ref dataset:")
for group in unique_groups:
    print(group)

print(f"\nTotal number of unique groups: {len(unique_groups)}")

def find_nearest_neighbors(query_embeddings, reference_embeddings, k=1):
    nn = NearestNeighbors(n_neighbors=k, metric='euclidean')
    nn.fit(reference_embeddings)
    distances, indices = nn.kneighbors(query_embeddings)
    return distances, indices

def annotate_cells(query_adata, reference_adata, embedding_key='X_uce', label_key='cell_type', k=1):
    distances, indices = find_nearest_neighbors(query_adata.obsm[embedding_key], reference_adata.obsm[embedding_key], k=k)
    if k == 1:
        query_adata.obs[f'predicted_{label_key}'] = reference_adata.obs[label_key].iloc[indices.flatten()].values
    else:
        # For k > 1, we'll use the most common label among the k nearest neighbors
        predicted_labels = []
        for idx_set in indices:
            labels = reference_adata.obs[label_key].iloc[idx_set].values
            predicted_labels.append(max(set(labels), key=list(labels).count))
        query_adata.obs[f'predicted_{label_key}'] = predicted_labels
    return query_adata

# load the query data
adata = sc.read_h5ad('/home/xinruqiu/302025/data/processed_data/302025_counts_uce_adata.h5ad')
adata = annotate_cells(adata, mouse_ref, embedding_key='X_uce', label_key='group', k=5) # Using 5 nearest neighbors
# Check results
print(adata.obs['predicted_group'].value_counts())
print(adata.obs['predicted_group'].unique())

# Create a dataframe with cell barcodes, predicted cell types, and coarse cell types
metadata_df = pd.DataFrame({
    'cell_barcode': adata.obs_names,
    'predicted_group': adata.obs['predicted_group']
    })

# Save the dataframe as a CSV file
output_file = '/home/xinruqiu/302025/data/processed_data/predicted_cell_types_metadata_mouseKidney.csv'
metadata_df.to_csv(output_file, index=False)
print(metadata_df.head())
