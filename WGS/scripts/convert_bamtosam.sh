dir_to_bams="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/" 
dir_to_sams="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/sams/mt_bams_filtered/"

cd $dir_to_bams

for bam in *.bam
do 
	echo $bam 
	path_to_sam="${dir_to_sams}${bam/.bam/.sam}" # replace the bam extension with sam while writing to file
	
	samtools view -bSq 20 $bam > "filtered_${bam}" # remove reads with MQ less than 20
	samtools view -h "filtered_${bam}" | awk 'length($10) > 40 || $1 ~ /^@/' | samtools view -b - > "filtered2_${bam}" # remove reads with insert size less than 40

	samtools view -h "filtered2_${bam}" > $path_to_sam # save the final bam file as sam after the filtering above
	mv "filtered2_${bam}" $bam # some renaming here

done

# remove intermediate files 
rm filtered_*






for bam in RG*.bam
do 
	echo $bam 
	mv "filtered2_${bam}" $bam # some renaming here

done
