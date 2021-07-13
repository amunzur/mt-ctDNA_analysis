#!/bin/bash
# this script is for converting the bam files containing mito reads to fastq, and realigning to HG using bwa.
bam_file=$1 # input, with "bam" extension
fq_path=$2 # output dir where all fq will be saved
ref=$3 # human ref genome
R1_name="${bam_file/.bam/_R1.fastq}"
R2_name="${bam_file/.bam/_R2.fastq}"
# sort before converting
samtools sort -n $bam_file -o sorted_RG_AE-132-WBC.bam
# convert to fasq
samtools fastq -@ 8 -1 "${fq_path}${R1_name}" -2 "${fq_path}${R2_name}" -0 /dev/null -s /dev/null -n sorted_RG_AE-132-WBC.bam

