#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/scripts/hunt_mito_reads/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};
conda activate alignment 

cd /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams

for bam_file in *.bam 
do
	mito_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito/${bam_file//.bam/_sorted.bam}"
	HG_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_HG/${bam_file//.bam/_sorted.bam}"

	printf "bash ${script_dir}/extract_compare_reads.sh ${bam_file} ${mito_bam} ${HG_bam}"
	sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/scripts/hunt_mito_reads/extract_compare_reads.sh ${bam_file} ${mito_bam} ${HG_bam};
done

# to test on one sample
# bam_file="CTRL-AE-201-WBC.bam"
# mito_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito/${bam_file//.bam/_sorted.bam}"
# HG_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_HG/${bam_file//.bam/_sorted.bam}"

# printf "bash ${script_dir}/extract_compare_reads.sh ${bam_file} ${mito_bam} ${HG_bam}"
# sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/scripts/hunt_mito_reads/extract_compare_reads.sh ${bam_file} ${mito_bam} ${HG_bam};

