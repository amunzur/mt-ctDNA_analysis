#!/bin/bash
# This script combined the vcf files in a given directory to one large file, great for using to take snapshots in igv later on 
vcf_dir="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/mutect2_results_filtered/"
cd $vcf_dir

for vcf in *vcf
do 
	egrep -v "^#" $vcf > "${vcf}_temp" # remove the header using reverse grep, assign to new file
	# sed -i "s/$/\t$f/" $vcf
	awk '{print FILENAME (NF?",":"") $0}' "${vcf}_temp" > "${vcf}_no_header"
	perl -pi -e "s/bam_FILTERED_vcf_temp/bam/g" "${vcf}_no_header"
done 

rm *_temp
cat *no_header > COMBINED.vcf # combine all headerless files into one large file