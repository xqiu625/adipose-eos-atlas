#!/bin/bash
#SBATCH --job-name=14b_export_mac
#SBATCH --output=/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/logs/14b_export_mac.out
#SBATCH --error=/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/logs/14b_export_mac.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=epyc

export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/cellrank_env/lib:/opt/linux/rocky/8.x/x86_64/pkgs/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
module load R/4.5.0
Rscript /bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts/14b_export_mac_h5ad.R
