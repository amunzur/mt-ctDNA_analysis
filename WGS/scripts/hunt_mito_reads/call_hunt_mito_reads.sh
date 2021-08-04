#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/scripts/hunt_mito_reads/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};
conda activate alignment 

cd /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams

for bam_file in *.bam 
do
	printf "bash ${script_dir}/hunt_mito_reads.sh ${bam_file}"
	sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/scripts/hunt_mito_reads/hunt_mito_reads.sh ${bam_file};
done

# bam_file="AE-156-WBC.bam"
# printf "bash ${script_dir}/hunt_mito_reads.sh ${bam_file} "
# sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/scripts/hunt_mito_reads/hunt_mito_reads.sh ${bam_file};

