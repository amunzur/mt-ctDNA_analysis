#!/bin/bash
mt_path="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams/"

for sample in /groups/wyattgrp/data/bam/WGS_PCa/alignments/*.bam
do 
	printf "\n"
	printf "Started sample ${sample}\n"

	sample_name=`basename $sample` # this helps strip the sample name from the full path

	samtools view -b ${sample} chrM > "${mt_path}${sample_name}" # we need the -b flag here otherwise it outputs an empty file
	samtools index "${mt_path}${sample_name}"
done