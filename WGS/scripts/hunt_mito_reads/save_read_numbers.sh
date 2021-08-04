#!/bin/bash 
path_to_readnums="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_numbers/realignment_bams/realignment_bams_awk.txt" 
rm $path_to_readnums
touch $path_to_readnums

cd /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams

for initial_bam in *bam 
do
	path_to_mito_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG/${initial_bam//.bam/.txt}"
	path_to_HG_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_mito/${initial_bam//.bam/.txt}"
	path_to_common_reads="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG_mito_common/${initial_bam//.bam/.txt}"
	path_to_common_reads_awk="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_names/realigned_HG_mito_common/${initial_bam//.bam/awk.txt}"

	echo $initial_bam
	echo $initial_bam $(wc -l $path_to_mito_reads | awk '{print $1}') $(wc -l $path_to_HG_reads | awk '{print $1}') $(wc -l $path_to_common_reads_awk | awk '{print $1}') >> $path_to_readnums
done;




