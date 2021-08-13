#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/scripts/VarScan2/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};
conda activate variant_calling 

# cd /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/mpileup/mt_bams_filtered
while read name_tumor name_normal
do
    echo $name_tumor
    echo $name_normal

	printf "bash ${script_dir}run_VarScan2.sh  ${name_tumor} ${name_normal}"
	sbatch --exclude=cn[01-05] ${script_dir}/run_VarScan2.sh  ${name_tumor} ${name_normal};

done < /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/scripts/VarScan2/bamsList.txt