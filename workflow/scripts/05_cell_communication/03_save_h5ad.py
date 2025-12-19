# Note: Activate conda environment before running: conda activate velocyto_env
# Working directory: /bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB

import pandas as pd
import anndata as ad
import numpy as np

# Read the CSV files
counts = pd.read_csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/hgnc_eos_mac_counts.csv", index_col=0)
metadata = pd.read_csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_mac_metadata.csv", index_col=0)

# Create AnnData object
adata = ad.AnnData(X=counts.T.values, obs=metadata, var=pd.DataFrame(index=counts.index))

# Save as h5ad
adata.write_h5ad("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_mac_cells.h5ad")


import pandas as pd
import anndata as ad
import numpy as np

# Read the CSV files
counts = pd.read_csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/hgnc_eos_mac_counts.csv", index_col=0)
metadata = pd.read_csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_mac_metadata_with_eos_subclusters.csv", index_col=0)

# Create AnnData object
adata = ad.AnnData(X=counts.T.values, obs=metadata, var=pd.DataFrame(index=counts.index))

# Save as h5ad
adata.write_h5ad("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_mac_cells_eos_subclusters.h5ad")
