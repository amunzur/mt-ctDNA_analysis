library(tidyverse)
library(erer)

################################################################
# MATTI PIPELINE
################################################################
path_to_combined_matti <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/bamsList" # this is the path that will contain the combined sample information with cfdna and wbc combinations 
path_to_combined_tnvstats <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/bamsList_tnvstats.txt"

bam_dir <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered"
cfdnaList_all <- grep(as.list(list.files(bam_dir, glob2rx("RG*.bam"))), pattern = "WBC|CTRL", invert=TRUE, ignore.case=TRUE, value=TRUE) # starting with RG, ending with .bam

# clean up the strings to retain the sample id only, this will help identify the corresponding wbc id
remove <- c(".bam", "-Baseline-cfDNA", "-Baseline-Tumor2", "-Baseline-Tumor", "-Progression-cfDNA", "-cfDNA", 
	"-1st-Progression-cfDNA", "-2017Sep28", "-2019Mar04", "-2018Nov26", "-Progression-Tumor2", "-Progression-Tumor", "-Progression2")

cfdnaList <- lapply(cfdnaList_all, function(some_name) str_remove_all(some_name, paste(remove, collapse = "|")))
wbcList_all <- grep(as.list(list.files(bam_dir, glob2rx("RG*.bam"))), pattern = "WBC", ignore.case=TRUE, value=TRUE) # starting with RG, ending with .bam, all WBC samples

# find the wbc matches 
wbcList <- list()
for (sample in cfdnaList) {

	wbc_sample <- grep(sample, wbcList_all, value=TRUE) # find the wbc match 
	wbcList <- append(wbc_sample, wbcList) # append the wbc samples to a list

	print(c(sample, wbc_sample))

}

# clean up the wbc matches 
wbcList <- unlist(lapply(as.list(wbcList), function(some_name) str_remove_all(some_name, ".bam")))
cfdnaList <- unlist(lapply(as.list(cfdnaList_all), function(some_name) str_remove_all(some_name, ".bam")))

combined <- data.frame("TEST" = sort(cfdnaList), "REF" = sort(wbcList))
write_delim(combined, path_to_combined_matti)

################################################################
# TNVSTATS
################################################################
combined$TEST <- paste(combined$TEST, "bam", sep=".")
combined$REF <- paste(combined$REF, "bam", sep=".")

write_delim(combined, path_to_combined_tnvstats)

# make sure to run the following code in the terminal to remove the first line of the file, which contains the colnames. 
setwd("/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/bams/mt_bams_filtered/")
x <- "tail -n +2 bamsList_tnvstats.txt > bamsList_tnvstats.txt.tmp && mv bamsList_tnvstats.txt.tmp bamsList_tnvstats.txt"

system(x) # run within the R console
