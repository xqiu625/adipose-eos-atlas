# 1. Create separate metadata files for each condition

import pandas as pd

# Read the original metadata file
metadata = pd.read_csv("eos_mac_metadata.csv")

# Assuming you have a 'condition' column
conditions = metadata['condition'].unique()

for condition in conditions:
    condition_metadata = metadata[metadata['condition'] == condition]
    condition_metadata.to_csv(f"metadata_{condition}.csv", index=False)

# 2. CellPhoneDB analysis script to run the analysis for each condition
import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

conditions = ['WT_Male_HFD', 'WT_Female_HFD', 'KO_Male_HFD', 'KO_Female_HFD','WT_Female_HFDInfected', 'KO_Female_HFDInfected'] 

for condition in conditions:
    print(f"Running analysis for {condition}")
    
    cpdb_statistical_analysis_method.call(
        cpdb_file_path = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/v5.0.0/cellphonedb.zip",
        meta_file_path = f"metadata_{condition}.csv",
        counts_file_path = "eos_mac_cells.h5ad",
        counts_data = 'hgnc_symbol',
        score_interactions = True,
        output_path = f"./cellphonedb_results_{condition}"
    )

print("All analyses completed.")
