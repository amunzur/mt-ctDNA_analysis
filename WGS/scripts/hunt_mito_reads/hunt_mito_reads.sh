#!/bin/bash

#SBATCH --job-name=Alignment-mito
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=31
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -t 12:00:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

bam_file=$1
R1_fasta="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/processed_fasta/${bam_file//.bam/.R1.fasta}" # string replacement
R2_fasta="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/processed_fasta/${bam_file//.bam/.R2.fasta}"

# NUCLEAR GENOME PARAMS
HG_sam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/sams/realigned_HG/${bam_file//.bam/.sam}" # sam file generated after aligning to human genome 
HG_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_HG/${bam_file}" 
HG_bam_sorted="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_HG/${bam_file//.bam/_sorted.bam}" 
path_to_ref_HG="/groups/wyattgrp/users/amunzur/chip_project/references/hg38_no_mito.fa" # chrM removed from hg38
path_to_flagstat_HG="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/flagstat/nuclear_genome_alignment/${bam_file//.bam/.txt}"

# MITO GENOME PARAMS
mito_sam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/sams/realigned_mito/${bam_file//.bam/.sam}" # sam file generated after aligning to human genome 
mito_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito/${bam_file}"  
mito_bam_sorted="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito/${bam_file//.bam/_sorted.bam}" 
path_to_ref_mito="/groups/wyattgrp/users/amunzur/mt-ctDNA/references/hg38_mito.fa"
path_to_flagstat_mito="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/flagstat/mito_genome_alignment/${bam_file//.bam/.txt}"

# align and index
align_and_samtools () {

	echo $5
	echo $6
	echo $7
	echo "Started aligning."
	bwa mem -t 6 -R "@RG\tID:sample" $1 $2 $3 > $4
	echo "Finished aligning."

	samtools view -bS $4 > $5
	samtools sort $5 -o $6
	samtools index $6
	samtools flagstat $6 > $7

}

echo "ALIGNING THE NUCLEAR GENOME."
align_and_samtools $path_to_ref_HG $R1_fasta $R2_fasta $HG_sam $HG_bam $HG_bam_sorted $path_to_flagstat_HG # nuclear

# echo "ALIGNING THE MITO GENOME."
# align_and_samtools $path_to_ref_mito $R1_fasta $R2_fasta $mito_sam $mito_bam $mito_bam_sorted $path_to_flagstat_mito # mito

# echo "PROCESSING THE NUCLEAR ALIGNMENT."
# filter_with_GATK $HG_bam filtered_${HG_bam}

# echo "PROCESSING THE MITO ALIGNMENT."

