! source activate /bigdata/godziklab/shared/Xinru/velocyto_env


import pandas as pd
import sys
import os

# -- Version of the databse
cpdb_version = 'v5.0.0'

# -- Path where the input files to generate the database are located
cpdb_target_dir = os.path.join('/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB', cpdb_version)
from cellphonedb.utils import db_utils
db_utils.download_database(cpdb_target_dir, cpdb_version)
pd.set_option('display.max_columns', 100)
os.chdir(cpdb_target_dir)


import sys
import os
import numpy as np
import pandas as pd
import cellphonedb
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

print(f"Python version: {sys.version}")
print(f"NumPy version: {np.__version__}")
print(f"Pandas version: {pd.__version__}")
print(f"CellPhoneDB version: {cellphonedb.__version__}")

# Print current working directory and list files
print(f"Current working directory: {os.getcwd()}")
print("Files in current directory:")
print(os.listdir('.'))


# Run CellPhoneDB analysis
cpdb_statistical_analysis_method.call(
    cpdb_file_path = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/v5.0.0/cellphonedb.zip",  # This is the default database file
    meta_file_path = "eos_mac_metadata.csv",
    counts_file_path = "eos_mac_cells.h5ad",
    counts_data = 'hgnc_symbol',
    score_interactions = True,
    output_path = "./cellphonedb_results"
)

