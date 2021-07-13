#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 13:45:38 2021
@author: amunzur
"""
import os
import pandas as pd

# in igv, write the range you want to visualize as an integer. if 100 for example, there will be 50 bp on both sides of the snv. 
# if no range is given, write an empty string: ""
given_range = 60
rewrite_file = True # rewrite batch script everytime we rerun this
show_WBC = True # do you also want to show the WBC samples next to tumor?
add_prefix = ""
add_suffix = ".bam"
input_type="csv" # options are vcf or csv, no capital letters 

# WINDOWS VERSION TO RUN IGV ON THE LOCAL MACHINE 
file_path = "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/results/somatic_cleaned.csv" # the location of the vcf of interest
snap_dir_windows = r"Y:\users\amunzur\mt-ctDNA\WGS\variant_calling\matti_pipeline_results\snapshots" # where snapshots will be saved
bam_dir = r"Y:\users\amunzur\mt-ctDNA\WGS/bams\mt_bams_filtered" # both WBC and tumor bams are here
path_batch = "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/snapshots/WINDOWS_batch_file_finland_bams" # where the batch file will be saved, same place as the png snapshots
path_bamList = "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/tumor_normal_pairs.txt"

if input_type == "vcf": 
	maf = pd.read_table(file_path, header = None) # read the vcf into a delim file, we call it a maf but honestly it really isnt a maf but rather a vcf 
	maf[['sample_t','chr']] = maf.iloc[:, 0].str.split(",", expand = True) # split the col with sample info into two cols 
	maf.columns[1] = "POS" # some renaming so that the col naming in maf and vcf files is consistent
else: 
	maf = pd.read_csv(file_path) # read the maf file as csv 

if rewrite_file and os.path.isfile(path_batch): 
	os.remove(path_batch)
	
for index, row in maf.iterrows(): # iterate through each snv in the maf
	tumor_id = row["sample_names"] # just the id
	bamtumor = os.path.join(bam_dir, str(add_prefix + tumor_id + add_suffix)).replace("/", "\\") # exact path to the tumor

	if show_WBC: # only consider if the user wants to see them, we need to refer to the bams list to find the corresponding wbc id
		bamList = pd.read_table(path_bamList, sep = "\t") # read the text file that corresponds the tumor bam to the wbc bams
		WBC_id = bamList.loc[bamList["TEST"] == tumor_id, ["REF"]].iloc[0, 0] # subset the df based on tumor id to get the wbc id
		bamWBC = os.path.join(bam_dir, str(add_prefix + WBC_id + add_suffix)).replace("/", "\\") # exact path to the wbc id

	position = int(row["POSITION"])
	if given_range: # if the user wants to see a given range, consider the end position of the range as well
		start_position = str(int(position - given_range/2))
		end_position = str(int(position + given_range/2))
		
	output_file_name = row["sample_names"] + "_" + str(row["CHROM"]) + "_" + str(position) + ".png" # one snapshot for each snv
	
	with open(path_batch, 'a') as the_file:
		the_file.write('new\n')
		if index == 0: the_file.write('genome hg38\n') # only load the genome if it is the first time we are starting IGV
		the_file.write(str("load " + bamtumor + '\n'))
		if show_WBC: the_file.write(str("load " + bamWBC + '\n'))
		the_file.write(str("snapshotDirectory " + snap_dir_windows + '\n'))
			
		chrom = str(row["CHROM"])
		if given_range: the_file.write(str('goto ' + chrom + ":" + start_position + "-" + end_position + "\n"))
		else: the_file.write(str('goto ' + chrom + ":" + str(position) + "\n"))
		the_file.write('sort base\n')
		the_file.write(str('snapshot ' + output_file_name + '\n'))
		the_file.write("\n")
			

