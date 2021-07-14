#!/bin/bash
# This script calculates depth for each sample in each position, and saves them to text files. 
# It also computes the mean coverage per sample and saves this information to ONE file  

dir_to_bams="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/" # INPUT
dir_to_read_depth="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth/" # the dir where read depth information will be saved, OUTPUT

# note that the read depth was calculated on samples after removing reads with low MQ and too short reads
#####################################################
# Calculate read depth - INDIVIDUAL samples
#####################################################

conda activate variant_calling # we need to use this samtools with the right version
cd $dir_to_bams # dir with all bams

for bam in *.bam
do 
    echo $bam 
    echo 
    path_to_metrics="${dir_to_read_depth}${bam}" # we keep the "bam" extension
    samtools depth $bam > $path_to_metrics

done

#####################################################
# Calculate read depth - COMBINED
#####################################################
dir_to_read_depth="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth/" # the dir where read depth information will be saved, OUTPUT
temp1="${dir_to_read_depth}temp1" # the file with sample names
temp2="${dir_to_read_depth}temp2" # the file with depth information, will be concatted with temp1 later on
combined="${dir_to_read_depth}AVERAGE_DEPTH_PER_SAMPLE"

cd $dir_to_read_depth # dir where the text files with depth information are saved

# clean up
rm $temp1 $temp2 # remove any earlier versions
touch $temp1 $temp2 # initialize an empty text file in the read depth directory

cd $dir_to_bams # dir to filtered bams
for bam in *.bam; do 
    read_num=$(samtools view $bam | wc -l)
    if (( $read_num > 0 )); then
        echo $bam # print the sample name to terminal 
        echo $bam >> $temp1; # sample name
        samtools depth $bam |  awk '{sum+=$3} END { print sum/NR}' >> $temp2; # average read depth from one sample
    fi
done


if [ `wc -l temp1` -eq `wc -l temp2` ]; then # if both have the same number of lines
    paste $temp1 $temp2 > $combined; # remove intermediate files only if the paste is successful
    rm $temp1 $temp2
fi

