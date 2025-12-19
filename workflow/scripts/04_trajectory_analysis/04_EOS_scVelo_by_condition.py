! source activate /bigdata/godziklab/shared/Xinru/velocyto_env


import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc


adata = scv.read("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata.h5ad")
# Create mapping dictionary from batch IDs to sample descriptions
batch_to_sample = {
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

# Map the batch column to create new sample_description column
adata.obs['sample_description'] = adata.obs['batch'].map(batch_to_sample)

# Verify the mapping worked correctly
print("Sample description value counts:")
print(adata.obs['sample_description'].value_counts())

# Check for any unmapped values (should be empty if all batches were mapped)
unmapped = adata.obs[adata.obs['sample_description'].isna()]
if len(unmapped) > 0:
    print(f"\nWarning: {len(unmapped)} cells have unmapped batch values:")
    print(unmapped['batch'].value_counts())
else:
    print("\nAll batch values successfully mapped to sample descriptions!")

# Optional: Extract additional categorical information from sample descriptions
adata.obs['genotype'] = adata.obs['sample_description'].str.split('_').str[0]  # WT or KO
adata.obs['sex'] = adata.obs['sample_description'].str.split('_').str[1]       # Male or Female
adata.obs['condition'] = adata.obs['sample_description'].str.extract(r'_(HFD(?:Infected)?)_')  # HFD or HFDInfected
adata.obs['infection_status'] = adata.obs['sample_description'].apply(
    lambda x: 'infection' if 'HFDInfected' in x else 'naive'
)
adata.obs['replicate'] = adata.obs['sample_description'].str.split('_').str[-1]  # Replicate number


# Subset the data by condition
ko_adata = adata[adata.obs['genotype'] == 'KO']
wt_adata = adata[adata.obs['genotype'] == 'WT']
female_adata = adata[adata.obs['genotype'] == 'KO']
male_adata = adata[adata.obs['genotype'] == 'WT']
naive_adata = adata[adata.obs['infection_status'] == 'naive']
infected_adata = adata[adata.obs['infection_status'] == 'Infected']


# Corrected boolean indexing with proper parentheses and correct values

# Use 'infection' and 'naive' (not 'Infected' and 'Not_Infected')
female_wt_infection_adata = adata[
    (adata.obs['sex'] == 'Female') & 
    (adata.obs['infection_status'] == 'infection') & 
    (adata.obs['genotype'] == 'WT')
].copy()

female_wt_naive_adata = adata[
    (adata.obs['sex'] == 'Female') & 
    (adata.obs['infection_status'] == 'naive') & 
    (adata.obs['genotype'] == 'WT')
].copy()

male_wt_naive_adata = adata[
    (adata.obs['sex'] == 'Male') & 
    (adata.obs['infection_status'] == 'naive') & 
    (adata.obs['genotype'] == 'WT')
].copy()

female_ko_infection_adata = adata[
    (adata.obs['sex'] == 'Female') & 
    (adata.obs['infection_status'] == 'infection') & 
    (adata.obs['genotype'] == 'KO')
].copy()

female_ko_naive_adata = adata[
    (adata.obs['sex'] == 'Female') & 
    (adata.obs['infection_status'] == 'naive') & 
    (adata.obs['genotype'] == 'KO')
].copy()

male_ko_naive_adata = adata[
    (adata.obs['sex'] == 'Male') & 
    (adata.obs['infection_status'] == 'naive') & 
    (adata.obs['genotype'] == 'KO')
].copy()


# Print the sizes of each subset to verify
print("Subset sizes:")
print(f"Female WT Infection: {female_wt_infection_adata.n_obs} cells")
print(f"Female WT Naive: {female_wt_naive_adata.n_obs} cells")
print(f"Male WT Naive: {male_wt_naive_adata.n_obs} cells")
print(f"Female KO Infection: {female_ko_infection_adata.n_obs} cells")
print(f"Female KO Naive: {female_ko_naive_adata.n_obs} cells")
print(f"Male KO Naive: {male_ko_naive_adata.n_obs} cells")

# Verify the total adds up to the original dataset
total_cells = (female_wt_infection_adata.n_obs + female_wt_naive_adata.n_obs + 
               male_wt_naive_adata.n_obs +  
               female_ko_infection_adata.n_obs + female_ko_naive_adata.n_obs +
               male_ko_naive_adata.n_obs )
print(f"\nTotal cells in subsets: {total_cells}")
print(f"Original dataset: {adata.n_obs} cells")
print(f"Match: {total_cells == adata.n_obs}")

# Show the breakdown by group
print("\nBreakdown by experimental groups:")
breakdown = adata.obs.groupby(['sex', 'genotype', 'infection_status']).size()
print(breakdown)

# Save the subsetted AnnData objects
ko_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_ko.h5ad")
wt_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_wt.h5ad")
female_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_female.h5ad")
male_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_male.h5ad")
naive_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_naive.h5ad")
infected_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_infected.h5ad")

female_wt_infection_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_female_wt_infection.h5ad")
female_wt_naive_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_female_wt_naive.h5ad")
male_wt_naive_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_male_wt_naive.h5ad")

female_ko_infection_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_female_ko_infection.h5ad")
female_ko_naive_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_female_ko_naive.h5ad")
male_ko_naive_adata.write("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_male_ko_naive.h5ad")




print(adata.obs['infection_status'].value_counts())

print(adata.obs['condition'].value_counts())
print(adata.obs['genotype'].value_counts())
print(adata.obs['sex'].value_counts())
print(adata.obs['infection_status'].value_counts())
print(adata.obs['replicate'].value_counts())




print(adata.obs['condition'].value_counts())