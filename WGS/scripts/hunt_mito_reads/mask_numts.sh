#!/bin/bash
# This script removes the regions from bam files already subsetted and filtered to contain mitochondrial reads only. 
# The bam files we modify come from the script named "extract_compare_reads.sh"

path_to_input_dir="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito_filtered"
path_to_output_dir="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/numts_blacklisted/" # bams will be saved here after numt regions are removed
path_to_bed="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/NUMTS/data/numts_hg38.bed"
path_to_metrics="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/CollectInsertSizeMetrics_metrics_numts_blacklisted/"
path_to_figures="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/insert_size_numts_blacklisted/"

cd ${path_to_input_dir}
for bam_file in RG_filtered_RG_*.bam
do 
	file="${bam_file//RG_filtered_RG_/}" # remove the unnecessarily long prefix

	echo $file
	bedtools intersect -abam ${file} -b ${path_to_bed} -v > "${path_to_output_dir}${file}"

	picard CollectInsertSizeMetrics \
      I="${path_to_output_dir}${file}" \
      O="${path_to_metrics}${file//.bam/.txt}" \
      H="${path_to_figures}${file//.bam/.pdf}" \
      M=0.5
done

# index all the files here
ls *.bam | parallel samtools index '{}'