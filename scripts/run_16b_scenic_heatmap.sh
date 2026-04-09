#!/bin/bash
#SBATCH --job-name=16b_scenic_heatmap
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/16b_scenic_heatmap.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/16b_scenic_heatmap.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=epyc

/bigdata/godziklab/shared/Xinru/cellrank_env/bin/python /bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts/16b_scenic_eos_mac_heatmap.py
