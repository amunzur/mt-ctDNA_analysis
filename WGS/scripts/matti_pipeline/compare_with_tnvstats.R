library(tidyverse)
library(knitr)
library(ggrepel)

# this script is for combining the positions from vcf and tnvstats files to make one large file

# helper functions
hh <- function(d) 
   if(class(d)=="matrix"|class(d)=="data.frame") d[1:20, 1:10]

get_vaf <- function(tnvstat){

	# Given a tnvstat, find the vaf from the allele frequencies and add this information as a new col to make it more legible.
	# Assuming the tnvstat already has a col named alt that contains the information about the identity of the alt allele. 
	
	tnvstat$tumor_vaf <- "dummy_data"
	tnvstat$normal_vaf <- "dummy_data"

	i <- 1
	while (i <= dim(tnvstat)[1]){

		i_row <- tnvstat[i, ]
		alt <- as.character(i_row$alt)
		# print(alt) # print to terminal what the mutant allele is 
		
		tumor_vaf <- paste0(alt, "AF_t") # go to the correct col to get the vaf based on what the alt is
		normal_vaf <- paste0(alt, "AF_n")# same as above but this one is for the wbc sample

		# subsetting here and adding to the df 
		idx <- which(names(i_row) == tumor_vaf) # get the idx of pos we are interested in 
		tumor_vaf <- as.character(i_row[idx]) # get the vaf based on that

		idx <- which(names(i_row) == normal_vaf)
		normal_vaf <- as.character(i_row[idx])

		# add to the df 
		i_row$tumor_vaf <- tumor_vaf
		i_row$normal_vaf <- normal_vaf

		# and modify the original df as well with the modified row as well 
		tnvstat[i, ] <- i_row

		i <- i + 1 # next row 

		} # end of while loop

	return(tnvstat) 

	} # end of function

path_to_cleaned_vcf <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results_realigned/mutations/somatic_cleaned.csv"
dir_to_tnvstats <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/tnvstats_realigned"
path_to_combined_tnvstats <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/tnvstats_realigned/combined.csv"

vcf <- as.data.frame(read_csv(path_to_cleaned_vcf))
vcf$sample_names <- paste0(vcf$sample_names, ".bam") # to make sure they match across tnvstats and the vcf file, add a prefix

pos_list <- as.list(vcf$POSITION) # extract a list of positions where the snvs occurred
alt_list <- as.list(vcf$ALT) # extract a list of mutant alleles

tumor_names <- lapply(as.list(vcf$sample_names), function(some_name) paste0(some_name)) # with the correct suffix
tnvstats_paths <- lapply(tumor_names, function(some_name) file.path(dir_to_tnvstats, some_name, paste0(some_name, ".tnvstat"))) # add the paths to the tumor ids 
tnvstats <- lapply(tnvstats_paths, read.delim) # save each tnvstat to a list as a data frame

some_func <- function(tnvstat, snv_pos){
	tnvstat <- tnvstat[!duplicated(tnvstat), ] # drop the duplicated rows 
	tnvstat <- as.data.frame(tnvstat)
	tnvstat <- list(filter(tnvstat, pos == snv_pos))

	# check if the pos exists in the tnvstats
	return(tnvstat)}

tnvstats_filtered <- mapply(some_func, tnvstats, pos_list) # for each tnvstat, only extract the pos where the snvs happened
combined <- do.call("rbind", tnvstats_filtered) # concatenate them all

# sanity check 
if (identical(as.character(combined$sample_t), as.character(vcf$sample_names)) == TRUE) {noquote("Fantastic job.")} else {"The order of samples is not the same across files."} # make sure they are identical
combined$alt <- vcf$ALT # get the alt allele from the matti's vcf
combined <- get_vaf(combined) # based on the allele freqs in the tnvstats, pull the vaf from normal and tumor samples

# further modify the combined df for plotting
combined <- as.data.frame(pivot_longer(data = combined, cols = c(tumor_vaf, normal_vaf), names_to = c("mut_type"))) # in the long format, now with a separate row for wbc and tumor for each sample
names(combined)[ncol(combined)] <- "mut_vafs" # need to rename the last col after pivot longer 
combined$mut_type <- ifelse(combined$mut_type == "tumor_vaf", "tumor", "normal") # a small ifelse to rename the entries in the mut_type column

# add a patient id and mut_id (2 extra columns) for easy plotting and legibility
x <- lapply(combined$sample_t, function(sample_id) head(str_split(sample_id, "-")[[1]], -2)) # just str split by - and disregard the last two elements
combined$patient_id <- unlist(lapply(x, function(sample_id) paste(sample_id[1], sample_id[2], sep = "-")))
combined$patient_id <- unlist(lapply(combined$patient_id, function(some_name) strsplit(some_name, "filtered_RG_")[[1]][[2]])) # remove the unnecessary prefix in front of the patient id
combined$mut_id <- paste(combined$patient_id, combined$pos, sep = "_")

# add the ctDNA percentage
samples_path <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/sample_information.csv" # a copy of the google sheets with sample information
samples <- as.data.frame(read_csv(samples_path)[, c(1, 5, 6)])
samples <- samples[-1, ] # the first row is the colname
names(samples) <- c("sample_id", "percentage_genePanel", "percentage_WGS") # rename the cols
samples$sample_id <- paste0("filtered_RG_", samples$sample_id, ".bam")# some string modifications to make sure the names match across data frames

samples <- samples %>%filter(sample_id %in% combined$sample_t)

combined <- inner_join(combined, samples, by = c("sample_t" = "sample_id")) # only keep the samples that are found in both dfs to have t
combined$percentage_WGS <- as.numeric(combined$percentage_WGS)
combined <- combined[order(combined$percentage_WGS),]

write_csv(combined, path_to_combined_tnvstats) # it's ok to overwrite!