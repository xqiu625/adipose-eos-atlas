import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors

# Load the human_ref sample dataset
human_ref = sc.read_h5ad('/home/xinruqiu/UCE/data/human_lung_disease.h5ad')
print(human_ref)
print(human_ref.obs_keys())  # Check available metadata
print(human_ref.var_names[:10])  # Check some gene names

eosinophil_mask = human_ref.obs['Subclass_Cell_Identity'] == 'Eosinophil'
eosinophil_count = eosinophil_mask.sum()
print(f"Number of eosinophil cells: {eosinophil_count}")
print(human_ref.obs['CellType_Category'].unique())
unique_groups = human_ref.obs['Subclass_Cell_Identity'].unique()
print("Unique groups in human_ref dataset:")
for group in unique_groups:
    print(group)

print(f"\nTotal number of unique groups: {len(unique_groups)}")

