#!/bin/bash

# script to prepare reads for bedtools
# cd /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/realigned_mito_filtered
# ls RG_filtered_RG_*.bam | parallel -j6 -k bash /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/scripts/pairize_reads.sh

# add paired read information, just adds /1 and /2 to the end of the read name
f=$1 
output_bam="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/paired_realigned_mito_filtered/" # where the bams with paired read info will be saved
output_bed="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/data/beds/"

echo $f 

# samtools sort -n $f | samtools view | awk 'BEGIN{FS=OFS="\t"} NR % 2 == 1 {$1=$1"\/1"} NR % 2 == 0 {$1=$1"\/2"} {print}' > "${output_path}${f}_temp1.sam" # headerless sam file
# cat <(samtools view -H $f) <(cat "${output_path}${f}_temp1.sam" ) > "${output_path}${f}_temp2.sam" # add header
# samtools view -S -b "${output_path}${f}_temp2.sam" > ${output_path}${f} # convert to bam and save

samtools sort -n ${f} | samtools view -hf 3 > ${output_bam}${f} # sort by name and remove singletons
bedtools bamtobed -i ${output_bam}${f} > ${output_bed}${f/.bam/.bed} # convert to bed
