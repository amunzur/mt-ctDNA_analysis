tumor_bam=$1 # path to tumor sample 
normal_bam=$2 # path to normal sample
main_dir=$3 # main dir where the magic happens
ref=$4 # path to human ref genome
name_tumor=$5 # name of the tumor sample
vardict_java_path=$6 # location of all vardict scripts
vardict_perl_path=$7 # location of all vardict scripts
outputdir="${main_dir}variant_calling/VarDict_results/" # output path where the vcf files from variant calling will be saved

# make sure end all paths with /, makes file concatenation easier
MT_bam_dir="${main_dir}bams/mt_bams/" # dir for saving bams with chrM only. should be done prior to running this script - INPUT
MT_filtered_dir="${main_dir}bams/mt_bams_filtered/" # dir where we save the bams with read group information, final output - OUTPUT

markdups_metric_dir="${main_dir}metrics/markDuplication_metrics/" # files with duplication metrics from picard are saved here
mutect2_output_dir="${main_dir}variant_calling/mutect2_results/" # vcf files from mutec
mutect2_filtered_output_dir="${main_dir}variant_calling/mutect2_results_filtered/" # filtered vcf files after filtering mutect results
read_numbers_dir="${main_dir}metrics/read_numbers/"

# cd $MT_filtered_dir # switch to the dir where we will be saving stuff, makes the naming below more legible

printf "\n"
printf "*******************************\n"
printf "RUNNING VARDICT - $bam_file\n"
printf "*******************************\n"

VarDict -G $ref -N $name_tumor -b "$tumor_bam|$normal_bam" -f 0.01 -R chrM:1-16569 | $vardict_perl_path/testsomatic.R | 
$vardict_perl_path/var2vcf_somatic.pl -N "$tumor_bam|$normal_bam" 0.01 > ${outputdir}/$name_tumor.vcf;


/home/amunzur/VarDictJava/VarDict/vardict.pl -G /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa \
-N RG_AE-132-Baseline-cfDNA \
-b "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/RG_AE-132-Baseline-cfDNA.bam|" \
-f 0.01 \
-R chrM:1-16569 | /home/amunzur/VarDictJava/VarDict/teststrandbias.R | /home/amunzur/VarDictJava/VarDict/testsomatic.R | /home/amunzur/VarDictJava/VarDict/var2vcf_paired.pl \
-N "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/RG_AE-132-Baseline-cfDNA.bam" \
-f 0.01 -x 200 \
> /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarDict_results/RG_AE-132-Baseline-cfDNA.vcf

VarDict -G /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa -N RG_AE-132-Baseline-cfDNA -b /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/RG_AE-132-Baseline-cfDNA.bam|/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/RG_AE-132-WBC.bam -f 0.01 -R chrM:1-16569 | $vardict_java_path/testsomatic.R | $vardict_java_path/var2vcf_paired.pl -N /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/RG_AE-132-Baseline-cfDNA.bam|/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/RG_AE-132-WBC.bam -f 0.01 > /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarDict_results/RG_AE-132-Baseline-cfDNA.vcf

