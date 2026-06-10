#!/usr/bin/bash

# 04_cellranger_multi.sh — Upstream raw-data processing, step 4 of 4
# SLURM `cellranger multi` (Cell Ranger 8.0.1) runs that demultiplex each of the
# 4 GEM-reaction pools (RR0104, RR0508, RR0912, RR1316) into per-sample count
# matrices using the config CSVs from 03_create_multiconfig_csv.sh.
# Output count matrices feed the analysis pipeline starting at 00_decontx_per_sample.R.
# Both standard and highmem SLURM variants are recorded; paths are HPCC-specific.

# RR1316

#SBATCH --ntasks=10
#SBATCH --mem=256gb
#SBATCH --time=4-00:00:00
#SBATCH -p batch
#SBATCH --job-name="RR1316"


cellranger multi --id=RR1316 --csv=/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR1316_config.csv

# RR1316 highmem

#SBATCH --ntasks=16
#SBATCH --mem=512gb
#SBATCH --time=4-00:00:00
#SBATCH -p highmem
#SBATCH --job-name="hRR1316"


cellranger multi --id=RR1316 --csv=/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR1316_config.csv


## RR0104


#SBATCH --ntasks=10
#SBATCH --mem=256gb
#SBATCH --time=4-00:00:00
#SBATCH -p intel
#SBATCH --job-name="RR0104"


cellranger multi --id=RR0104 --csv=/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR0104_config.csv

## RR0104 highmem

#SBATCH --ntasks=16
#SBATCH --mem=512gb
#SBATCH --time=4-00:00:00
#SBATCH -p highmem
#SBATCH --job-name="hRR0104"


cellranger multi --id=RR0104 --csv=/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR0104_config.csv



## RR0508


#SBATCH --ntasks=10
#SBATCH --mem=256gb
#SBATCH --time=4-00:00:00
#SBATCH -p intel
#SBATCH --job-name="RR0508"


cellranger multi --id=RR0508 --csv=/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR0508_config.csv


## RR0912


#SBATCH --ntasks=16
#SBATCH --mem=256gb
#SBATCH --time=10-00:00:00
#SBATCH -p intel
#SBATCH --job-name="RR0912"


cellranger multi --id=RR0912 --csv=/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR0912_config.csv
