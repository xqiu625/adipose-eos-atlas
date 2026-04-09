#!/bin/bash
#SBATCH --job-name=10c_cellchat_fig3d
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/10c_cellchat_fig3d.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302025/06_LOGS/10c_cellchat_fig3d.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=epyc

export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/cellrank_env/lib:/opt/linux/rocky/8.x/x86_64/pkgs/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
module load R/4.5.0
Rscript /bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts/10c_cellchat_fig3d.R
