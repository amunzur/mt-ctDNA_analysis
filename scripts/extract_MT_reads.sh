#!/bin/bash
input_path="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/original_bams"
mt_path="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams"

for sample in /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/original_bams/*.bam
do 
	printf "\n"
	printf "Started sample ${sample}\n"

	samtools view ${sample} chrM > ${mt_path}${sample}
	samtools index ${mt_path}${sample}
done