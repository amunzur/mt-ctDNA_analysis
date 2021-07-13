
# This script calculates depth for each sample in each position, and saves them to text files. 
# It also computes the mean coverage per sample and saves this information to ONE file  

dir_to_bams="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/" # INPUT
dir_to_read_depth="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth/" # the dir where read depth information will be saved, OUTPUT

# note that the read depth was calculated on samples after removing reads with low MQ and too short reads
#####################################################

#####################################################

conda activate variant_calling # we need to use this samtools with the right version
cd $dir_to_bams # dir with all bams

for bam in *.bam
do 
	echo $bam 
	path_to_metrics="${dir_to_read_depth}${bam/.bam/""}" # remove the bam extension, just write a normal text file
	samtools depth $bam > $path_to_metrics

done

# and this one piece of code calculates the average read depth per sample 
cd $dir_to_read_depth # dir where the text files with depth information are saved
touch "${dir_to_read_depth}ALL_SAMPLES_READ_DEPTH" # initialize an empty text file in the read depth directory

cd $dir_to_bams

for bam in *.bam
do 
	echo $bam 
	path_to_metrics="${dir_to_read_depth}${bam/.bam/""}" # remove the bam extension, just write a normal text file
	samtools depth $bam > $path_to_metrics

done


# remove intermediate files 
rm filtered*
