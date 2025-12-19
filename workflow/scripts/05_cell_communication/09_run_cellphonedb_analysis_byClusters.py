import pandas as pd
import os
import anndata as ad
import pandas as pd

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# Path to your counts file
counts_file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_mac_cells_eos_subclusters.h5ad"  
cpdb_file_path = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/v5.0.0/cellphonedb.zip"  

# Read the metadata
metadata = pd.read_csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_mac_metadata_with_eos_subclusters.csv", index_col=0)

# Load the existing h5ad file
adata = ad.read_h5ad(counts_file)

# Get unique eosinophil clusters
eos_clusters = metadata[metadata['cell_type'] == 'Eosinophils']['Clusters'].unique()


for eos_cluster in eos_clusters:
    print(f"Analyzing Eosinophil cluster {eos_cluster} vs Macrophages")
    
    # Create a subset of metadata for this analysis
    subset_metadata = metadata[(metadata['Clusters'] == eos_cluster) | (metadata['Clusters'] == 'Mac')].copy()
    subset_metadata.loc[:, 'cell_type'] = subset_metadata['Clusters'].map({'Mac': 'Macrophage'}).fillna(f'Eos_{eos_cluster}')
    subset_metadata_path = f"metadata_Eos{eos_cluster}_vs_Mac.csv"
    subset_metadata.to_csv(subset_metadata_path)
    
    output_path = f"./cellphonedb_results_Eos{eos_cluster}_vs_Mac"
    
    # Ensure output directory exists
    os.makedirs(output_path, exist_ok=True)
    
    try:
        # Run CellPhoneDB
        cpdb_statistical_analysis_method.call(
            cpdb_file_path=cpdb_file_path,
            meta_file_path=subset_metadata_path,
            counts_file_path=counts_file,
            counts_data='hgnc_symbol',
            score_interactions=True,
            output_path=output_path
        )
        print(f"Analysis completed for Eosinophil cluster {eos_cluster}")
    except Exception as e:
        print(f"Error analyzing Eosinophil cluster {eos_cluster}: {str(e)}")

print("All analyses completed.")
