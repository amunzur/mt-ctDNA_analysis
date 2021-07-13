#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/config.txt" # no need to modify
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh" # no need to modify
source ${configfile};
source ${conda_profile_path};
conda activate variant_calling 

if [[ $type == "WGS" ]] # comes from the the config file
then 
	main_dir="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/"
	cd "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams"	
fi

if [[ $test == "YES" ]] #comes from the config file
then 
	cd "${main_dir}bams/test_bams"	
fi

for bam_file in *.bam # for all bam files in the current dir, so we need the "cd" command above
do
	printf "bash ${script_dir}/process_samples.sh ${bam_file} ${main_dir} "
	sbatch --exclude=cn[01-05] ${path_to_script} ${bam_file} ${main_dir}; # path_to_script is specified in the config file
done
