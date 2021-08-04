#!/bin/bash

#SBATCH --job-name=extract
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=31
#SBATCH --mem 737000 # memory pool for all cores
#SBATCH -t 12:00:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

# This script extracts read ids from mito and HG alignments, and then finds the common ones. 

initial_bam=$1 # initial bam from finland filtered to keep the mito reads only
mito_bam=$2 # $1 aligned to mito
HG_bam=$3 # $1 aligned to nuclear genome

# save text files here
path_to_HG_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG60//${initial_bam//.bam/.txt}"
path_to_mito_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_mito/${initial_bam//.bam/.txt}"
path_to_common_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG_mito_common/${initial_bam//.bam/.txt}" # one text file per sample, containing the common reads
dir_to_common_reads_SPLIT="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG_mito_common_SPLIT/${initial_bam}/" # splitted files containing reasd names, make sure to end with /

path_to_common_reads_awk="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG_mito_common/${initial_bam//.bam/awk.txt}"

path_to_mito_filtered_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito_filtered/${initial_bam}" # and this is the filtered mito bam files - without reads that appeared in both alignments - MERGED BAM 
dir_to_mito_filtered_bam_SPLIT="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito_filtered_SPLIT_BAMS/${initial_bam}/" # after grepping multiple text files with read names 

path_to_figures="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/insert_size_mito_filtered/${initial_bam//.bam/.pdf}" # picard fragment size plots
path_to_metrics="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/CollectInsertSizeMetrics_mito_filtered/${initial_bam//.bam/.txt}" # picard fragment size metrics

path_to_metrics_rmdup="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/markDuplication_metrics/realigned_mito_filtered/${initial_bam//.bam/.txt}"

# echo "Extracting reads from mito bams."
# samtools view $2 |  cut -f1 | LC_ALL=C sort > $path_to_mito_reads

# echo "Extracting reads from nuclear bams."
# samtools view $3 |  cut -f1 | LC_ALL=C sort > $path_to_HG_reads
# samtools view -q 60 $3 |  cut -f1 | LC_ALL=C sort > $path_to_HG_reads

# echo "Finding common reads."
# LC_ALL=C comm -12 $path_to_mito_reads $path_to_HG_reads > $path_to_common_reads
# awk 'NR==FNR{seen[$0]=1; next} seen[$0]' $path_to_mito_reads $path_to_HG_reads > $path_to_common_reads_awk

# make the dirs we need 
mkdir -p $dir_to_mito_filtered_bam_SPLIT # split bams are saved here 
mkdir -p $dir_to_common_reads_SPLIT # split read names txt files are saved here

# divide the commond reads file into smaller chunks 
echo "Splitting the common reads file. New files are here: $dir_to_common_reads_SPLIT"
read_nums=`wc -l $path_to_common_reads | awk '{print $1}'`
split -l $(($read_nums/5)) $path_to_common_reads $dir_to_common_reads_SPLIT

echo "Started removing the common reads from initial bams." "The file we are working on is" $initial_bam
samtools view -h $initial_bam | grep -vf "${dir_to_common_reads_SPLIT}aa" | samtools view -bS -o "${dir_to_mito_filtered_bam_SPLIT}${initial_bam}aa" 
echo "Finished file 1. Saved the bam to ${dir_to_mito_filtered_bam_SPLIT}${initial_bam}aa"
samtools view -h $initial_bam | grep -vf "${dir_to_common_reads_SPLIT}ab" | samtools view -bS -o "${dir_to_mito_filtered_bam_SPLIT}${initial_bam}ab"
echo "Finished file 2. Saved the bam to ${dir_to_mito_filtered_bam_SPLIT}${initial_bam}aa"
samtools view -h $initial_bam | grep -vf "${dir_to_common_reads_SPLIT}ac" | samtools view -bS -o "${dir_to_mito_filtered_bam_SPLIT}${initial_bam}ac"
echo "Finished file 3. Saved the bam to ${dir_to_mito_filtered_bam_SPLIT}${initial_bam}aa"
samtools view -h $initial_bam | grep -vf "${dir_to_common_reads_SPLIT}ad" | samtools view -bS -o "${dir_to_mito_filtered_bam_SPLIT}${initial_bam}ad"
echo "Finished file 4. Saved the bam to ${dir_to_mito_filtered_bam_SPLIT}${initial_bam}aa"
samtools view -h $initial_bam | grep -vf "${dir_to_common_reads_SPLIT}ae" | samtools view -bS -o "${dir_to_mito_filtered_bam_SPLIT}${initial_bam}ae"
echo "Finished file 5. Saved the bam to ${dir_to_mito_filtered_bam_SPLIT}${initial_bam}aa"
samtools view -h $initial_bam | grep -vf "${dir_to_common_reads_SPLIT}af" | samtools view -bS -o "${dir_to_mito_filtered_bam_SPLIT}${initial_bam}af"
echo "Finished file 6. Saved the bam to ${dir_to_mito_filtered_bam_SPLIT}${initial_bam}aa"

echo "Started merging split bam files."
cd ${dir_to_mito_filtered_bam_SPLIT}${initial_bam} # go into sample dir with all the split bam files
samtools merge -O $path_to_mito_filtered_bam ${dir_to_mito_filtered_bam_SPLIT}* # merge all files here
echo "Finished merging. Merged bam is here: $path_to_mito_filtered_bam"

# duplicate removal
echo "Started processing the merged bam with picard."
picard MarkDuplicates I=$path_to_mito_filtered_bam​ O="${path_to_mito_filtered_bam}.rmdup" M="${path_to_metrics_rmdup}" REMOVE_DUPLICATES=true
picard FixMateInformation I="${path_to_mito_filtered_bam}.rmdup" O="${path_to_mito_filtered_bam}.fixmate" ADD_MATE_CIGAR=true
picard AddOrReplaceReadGroups \
    I="${path_to_mito_filtered_bam}.fixmate" \
    O="${path_to_mito_filtered_bam}.RG" \
    RGID=1 \
    RGLB=library \
    RGPL=ILLUMINA \
    RGPU=unit \
    RGSM=sample

# run picard to get the fragment sizes
echo "Running picard to generate figures."
picard CollectInsertSizeMetrics \
      I=${path_to_mito_filtered_bam} \
      O=${path_to_metrics} \
      H=${path_to_figures} \
      M=0.5
