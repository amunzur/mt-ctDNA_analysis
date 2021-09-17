#!/bin/bash
#SBATCH --job-name=WPS_mito
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=2
#SBATCH --mem 100000 # memory pool for all cores
#SBATCH -t 05:30:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

# Script to call compute_WPS.sh
bed_file="$1"
window_size="$2"

conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${conda_profile_path};
conda activate r_env_v1
printf "Rscript /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/scripts/compute_WPS.R  --path-to-bed ${bed_file} --window-size ${window_size}\n";
Rscript /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/scripts/compute_WPS.R  --path-to-bed ${bed_file} --window-size ${window_size};
