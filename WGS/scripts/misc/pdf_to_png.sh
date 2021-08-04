#!/bin/bash 
# This script converts all pdf in a directory to png.

pdf_dir="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/insert_size_mito_filtered" # INPUT
png_dir="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/insert_size_mito_filtered_PNG" # OUTPUT

cd $pdf_dir

# convert 
for file in *pdf
do
	echo $file
	pdftoppm  $file ${file//.pdf/} -png
	mv "${file//.pdf/}-1.png" "${file//.pdf/.png}"
done

# move the png to a new dir
find . -name "*png" -type f | xargs mv -t $png_dir/