#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};
conda activate mito 

if [[ $type == "WGS" ]]
then 
	main_dir="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/"
	cd "${main_dir}bams/original_bams"	
fi

if [[ $test == "YES" ]]
then 
	cd "${main_dir}bams/test_bams"	
fi

for bam_file in *.bam 
do
	printf "bash ${script_dir}/process_samples.sh ${bam_file} ${main_dir} "
	sbatch --exclude=cn[01-05] ${path_to_script} ${bam_file} ${main_dir};
done
