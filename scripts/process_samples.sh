#!/bin/bash

#SBATCH --job-name=Mutect2
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=31
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -t 08:30:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

bam_file=$1
main_dir=$2

markdups_metric_dir="${main_dir}metrics/markDuplication_metrics/"
RG_bam_dir="${main_dir}bams/rg_bams/"
MT_bam_dir="${main_dir}bams/mt_bams/" # dir for saving bams with chrM only
mutect2_output_dir="${main_dir}variant_calling/mutect2_results/" # vcf files from mutec
mutect_filtered_output_dir="${main_dir}variant_calling/mutect2_results_filtered/" # filtered vcf files after filtering mutect results
read_numbers_dir="${main_dir}metrics/read_numbers/"

cd RG_bam_dir # switch to the dir where we will be saving stuff, makes the naming below more legible

printf "\n"
printf "*******************************\n"
printf "MARK DUPLICATES - $bam_file\n"
printf "*******************************\n"

picard MarkDuplicates I=$bam_file O="rmdup_${bam_file}" M="${markdups_metric_dir}rmdup_${bam_file}" REMOVE_DUPLICATES=true

printf "\n"
printf "*******************************\n"
printf "FIXMATE - $bam_file\n"
printf "*******************************\n"

picard FixMateInformation I="rmdup_${bam_file}" O="fixed_mate_${bam_file}" ADD_MATE_CIGAR=true

picard AddOrReplaceReadGroups \
    I="fixed_mate_${bam_file}" \
    O="RG_${bam_file}" \
    RGID=1 \
    RGLB=library \
    RGPL=ILLUMINA \
    RGPU=unit \
    RGSM=sample

samtools index "RG_${bam_file}" # index the final product

samtools view -b "RG_${bam_file}" chrM > "$MT_bam_dir${bam_file}" # filter to keep chrM only
samtools index "$MT_bam_dir${bam_file}" # index after filtering

rm rmdup* fixed* # remove intermediate files

printf "\n"
printf "******************************\n*"
printf "VARIANT CALLING\n"
printf "*******************************\n"

/home/amunzur/gatk-4.2.0.0/gatk Mutect2 \
-R /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa \
--mitochondria-mode \
-I "$MT_bam_dir${bam_file}" \
-O "$mutect2_output_dir${bam_file}_vcf"

printf "\n"
printf "******************************\n*"
printf "FILTERING CALLED VARIANTS\n"
printf "*******************************\n"

/home/amunzur/gatk-4.2.0.0/gatk FilterMutectCalls \
-R /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa \
-V "$mutect2_output_dir${bam_file}_vcf" \
-O "$mutect_filtered_output_dir${bam_file}_FILTERED_vcf"

printf "\n"
printf "******************************\n*"
printf "COUNTING READS\n"
printf "*******************************\n"

samtools view -c "${bam_file}" >> "${read_numbers_dir}original_bams/${bam_file}" # original bam
samtools view -c "RG_${bam_file}" >> "${read_numbers_dir}rg_bams/${bam_file}" # after removing duplicates
samtools view -c "$MT_bam_dir${bam_file}" >> "${read_numbers_dir}mt_bams/${bam_file}" # after removing duplicates

