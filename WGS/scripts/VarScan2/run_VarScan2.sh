#!/bin/bash 

#SBATCH --job-name=VarScan2
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=31
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -t 05:30:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

tumor=$1
normal=$2
dir_to_mpileups="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/mpileup/mt_bams_filtered/"
path_vcf_snv="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarScan2/snv/${tumor/.mpileup/.vcf}"
path_vcf_indel="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarScan2/indel/${tumor/.mpileup/.vcf}"

java -jar /home/amunzur/VarScan.v2.3.9.jar somatic "${dir_to_mpileups}${normal}" "${dir_to_mpileups}${tumor}" ${tumor/.mpileup/.vcf} \
--output-snp $path_vcf_snv \
--output-indel $path_vcf_indel \
--min-coverage-normal 8 \
--min-coverage-tumor 4 \
--min_var_freq 0.20 \
--p-value 0.05 \
--somatic-p-value 0.05