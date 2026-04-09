#!/bin/bash
#SBATCH --job-name=16b_scenic_heatmap
#SBATCH --output=/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/logs/16b_scenic_heatmap.out
#SBATCH --error=/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/logs/16b_scenic_heatmap.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=epyc

source /bigdata/godziklab/shared/Xinru/cellrank_env/bin/activate
python /bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts/16b_scenic_eos_mac_heatmap.py
