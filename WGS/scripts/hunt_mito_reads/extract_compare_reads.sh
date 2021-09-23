#!/bin/bash

#SBATCH --job-name=extract
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=6
#SBATCH --mem 40000 # memory pool for all cores
#SBATCH -t 5:00:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

# This script extracts read ids from mito and HG alignments, and then finds the common ones. 

initial_bam=$1 # initial bam from finland filtered to keep the mito reads only
mito_bam=$2 # $1 aligned to mito
HG_bam=$3 # $1 aligned to nuclear genome

# save text files here
path_to_HG_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG60/${initial_bam//.bam/.txt}"
path_to_mito_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_mito/${initial_bam//.bam/.txt}"
path_to_common_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG_mito_common/${initial_bam//.bam/.txt}" # one text file per sample, containing the common reads

path_to_mito_filtered_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito_filtered2/${initial_bam}" # and this is the filtered mito bam files - without reads that appeared in both alignments 

# OUTPUT FILES FOR PICARD
path_to_figures="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/insert_size_mito_filtered2/${initial_bam//.bam/.pdf}" # picard fragment size plots - PDF
path_to_figures_png="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/insert_size_mito_filtered_PNG2/${initial_bam//.bam/}" # picard fragment size plots - PNG

path_to_metrics="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/CollectInsertSizeMetrics_mito_filtered2/${initial_bam//.bam/.txt}" # picard fragment size metrics

path_to_metrics_rmdup="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/markDuplication_metrics/realigned_mito_filtered2/${initial_bam//.bam/.txt}"

echo "Extracting reads from mito bams."
samtools view -q60 $2 |  cut -f1 | LC_ALL=C sort > $path_to_mito_reads

echo "Extracting reads from nuclear bams."
samtools view -q50 $3 |  cut -f1 | LC_ALL=C sort > $path_to_HG_reads
# samtools view -q 60 $3 |  cut -f1 | LC_ALL=C sort > $path_to_HG_reads

echo "Finding common reads."
LC_ALL=C comm -12 $path_to_mito_reads $path_to_HG_reads > $path_to_common_reads
# awk 'NR==FNR{seen[$0]=1; next} seen[$0]' $path_to_mito_reads $path_to_HG_reads > $path_to_common_reads_awk

echo "Started filtering reads from mito bams with picard from " $initial_bam
picard -Xmx40g FilterSamReads \
    I=$initial_bam \
    O=$path_to_mito_filtered_bam \
    READ_LIST_FILE=$path_to_common_reads\
    FILTER=excludeReadList 

# duplicate removal
echo "Started processing the merged bam with picard."

picard MarkDuplicates I=$path_to_mito_filtered_bam O="${path_to_mito_filtered_bam}.rmdup" M="${path_to_metrics_rmdup}" REMOVE_DUPLICATES=true

picard FixMateInformation I="${path_to_mito_filtered_bam}.rmdup" O="${path_to_mito_filtered_bam}.fixmate" ADD_MATE_CIGAR=true

picard AddOrReplaceReadGroups \
    I="${path_to_mito_filtered_bam}.fixmate" \
    O="${path_to_mito_filtered_bam}.RG" \
    RGID=1 \
    RGLB=library \
    RGPL=ILLUMINA \
    RGPU=unit \
    RGSM=sample

# samtools sort "${path_to_mito_filtered_bam}.sorted" "${path_to_mito_filtered_bam}.RG"

# some renaming right here 
# for file in *.fixmate; do mv -- "$file" "${file%.log}.txt"; done

# for i in $( ls *.bam ); do mv $i ${i%.*}; done

# run picard to get the fragment sizes
echo "Running picard to generate figures."
picard CollectInsertSizeMetrics \
      I="${path_to_mito_filtered_bam}.RG" \
      O=${path_to_metrics} \
      H=${path_to_figures} \
      M=0.5