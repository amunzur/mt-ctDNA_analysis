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

path_to_cleaned_vcf <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/results/somatic_cleaned.csv"
dir_to_tnvstats <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/tnvstats"
path_to_combined_tnvstats <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/tnvstats/combined.csv"

vcf <- as.data.frame(read_csv(path_to_cleaned_vcf))
vcf$sample_names <- paste0(vcf$sample_names, ".bam") # to make sure they match across tnvstats and the vcf file

pos_list <- as.list(vcf$POSITION) # extract a list of positions where the snvs occurred
alt_list <- as.list(vcf$ALT) # extract a list of mutant alleles

tumor_names <- lapply(as.list(vcf$sample_names), function(some_name) paste0(some_name)) # with the correct suffix
tnvstats_paths <- lapply(tumor_names, function(some_name) file.path(dir_to_tnvstats, some_name, paste0(some_name, ".tnvstat"))) # add the paths to the tumor ids 
tnvstats <- lapply(tnvstats_paths, read.delim) # save each tnvstat to a list as a data frame

some_func <- function(tnvstat, snv_pos){
	tnvstat <- tnvstat[!duplicated(tnvstat), ] # drop the duplicated rows 
	tnvstat <- as.data.frame(tnvstat)
	tnvstat <- list(filter(tnvstat, pos == snv_pos))
	return(tnvstat)}

tnvstats_filtered <- mapply(some_func, tnvstats, pos_list) # for each tnvstat, only extract the pos where the snvs happened
combined <- do.call("rbind", tnvstats_filtered) # concatenate them all

write_delim(combined, path_to_combined_tnvstats) # save the cleaned up tnvstats

##########################################################################
if (identical(as.character(combined$sample_t), as.character(vcf$sample_names)) == TRUE) {noquote("Fantastic job.")} else {"The order of samples is not the same across files."} # make sure they are identical
combined$alt <- vcf$ALT
combined <- get_vaf(combined) # based on the allele freqs in the tnvstats, pull the vaf from normal and tumor samples

# some minor modifications to fix the naming
combined <- combined %>% mutate_all(funs(str_replace(., "RG_DTB-149-Baseline-Tumor2.bam", "RG_DTB-149-Baseline-Tumor.bam")))

# save the tnvstat & vcf file
write_csv(combined, path_to_combined_tnvstats) # it's ok to overwrite!

########################################################################## not done just yet! 
# we also need the allele frequencies of all samples at the positions we identified. The following 
# loads all tnvstats from all samples, extracts the position and the allele frequency for the alt allele. 
# This happends regardless of whether they have the mutation at that location or not. 
some_func_all <- function(tnvstat, snv_pos_list){

	# This is slighly different from the some_func above. This one will extract ALL mutations
	# from all samples, regardless of whether they have a mut at that position or not
	
	tnvstat <- as.data.frame(tnvstat[!duplicated(tnvstat), ]) # drop the duplicated rows from a given tnvstat file
	tnvstat <- list(tnvstat[tnvstat$pos %in% snv_pos_list ,]) # only this row is different from the initial function

	return(tnvstat)}

add_variant_information <- function(tnvstat, vcf) {

	# For a tnvstat file, go through the vcf, identify the correct sample, find the related mutations from that sample and 
	# add them to the tnvstat. Also make a new col to indicate which snvs are muts and which are not.
	# This function also adds what the alt is, and finds the corresponding freq information from other samples without the mutation.

	tumor <- as.character(unique(tnvstat$sample_t)) # tumor name taken from the tnvstat
	tnvstat$mut_status <- "normal"

	vcf_filtered <- vcf %>% filter(sample_names == tumor) # filter to keep the tumor id only from the whole vcf
	
	if (dim(vcf_filtered)[1] > 0) { # if we have at least one mutation from a given sample

		tnvstat$mut_status[which(tnvstat$pos %in% vcf_filtered$POSITION)] <- "variant" # find the idx in the tnvstat where the variant exists, turn it from normal to variant

		# identify the alt from the vcf, unduplicated form
		idx <- which(duplicated(vcf$POSITION))
		tnvstat$alt <- vcf$ALT[-idx] # add the alt alleles to the tnvstats as an extra col 

		# now identify the vafs of the variant, and the general freq of the allele across all samples
		i <- 1
		tumor_vaf_list <- list()
		normal_vaf_list <- list()
		while (i <= dim(tnvstat)[1]) { # go through the rows (variants) in the tnvstat to identify the vaf based on the alt

			colname <- paste0(tnvstat$alt[i], "AF_t")
			tumor_vaf <- tnvstat[i, colname] # identify the vaf 
			tumor_vaf_list <- append(tumor_vaf_list, tumor_vaf) # add to list

			colname <- paste0(tnvstat$alt[i], "AF_n")
			normal_vaf <- tnvstat[i, colname]
			normal_vaf_list <- append(normal_vaf_list, normal_vaf)

			i <- i + 1

			} # end of while loop 

		# add the previously calculated tumor and normal vafs to the tnvstat 
		tnvstat$tumor_vaf <- unlist(tumor_vaf_list)
		tnvstat$normal_vaf <- unlist(normal_vaf_list)

	return(list(tnvstat))

	} # end of if loop

} # end of function

tnvstats_filtered <- mapply(some_func_all, tnvstats, rep(list(unlist(pos_list)), length(pos_list))) # for each tnvstat, extract all snv positions, regardless of whether they occurred in that sample or not
tnvstats_filtered2 <- mapply(add_variant_information, tnvstats_filtered, rep(list(vcf), length(tnvstats_filtered))) # for each tnvstat, identify the variants from the rest of the snvs
# At this point we made tnvstats for each individual sample, computed vaf for tumor and normal for each snv position we identified. 

all_snvs <- do.call("rbind", tnvstats_filtered2)
all_snvs$mut_status <- as.factor(all_snvs$mut_status)

##########################################################################
# PLOTTING - utilities

cool_theme <- 
  
  theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black", size = 1), 
      axis.ticks = element_line(colour = "black", size = 2),
      axis.text = element_text(size=10),
      axis.text.x = element_text(vjust=0.5, colour = "black", size=8),
      axis.text.y = element_text(vjust=0.5, colour = "black", size=8),
      axis.title = element_text(size=10,face="bold"), 
      legend.title = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 12),
      axis.ticks.length=unit(0.15, "cm"), 
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 20)), 
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

##########################################################################
# PLOTTING ALL SNVS IN TWO SEPARATE PLOTS - not super informative 
##########################################################################
tumor_vaf_plot <- ggplot(data = all_snvs, aes(x = sample_t, y = tumor_vaf, color = mut_status)) + 
	geom_point(position=position_jitterdodge())	+ 
	xlab("Sample name") + 
	ylab("Allele frequency") +
	ggtitle("TUMOR") +
	theme_bw() + 
	cool_theme +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

normal_vaf_plot <- ggplot(data = all_snvs, aes(x = sample_t, y = normal_vaf, color = mut_status)) + 
	geom_point(position=position_jitterdodge())	+ 
	xlab("Sample name") + 
	ylab("Allele frequency") +
	ggtitle("NORMAL") +
	theme_bw() + 
	cool_theme +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

path_to_tumor_plot <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/figures/tumor_allele_freq.pdf"
path_to_normal_plot <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/figures/normal_allele_freq.pdf"

ggsave(path_to_tumor_plot, tumor_vaf_plot, width = 10, height = 7, units = "in")
ggsave(path_to_normal_plot, normal_vaf_plot, width = 10, height = 7, units = "in")

##########################################################################
# MAKING A CLASSIC VARIANT PLOT TO COMPARE TUMOR AND NORMAL VAFS
##########################################################################
# For this figure we go back to the df made earlier: combined. This one contains 
# all the identified variants, with vaf and sample information. 

# some minor modifications so that the df is in the right format for plotting 
mut_type <- c(rep("tumor", dim(combined)[1]), rep("normal", dim(combined)[1]))
mut_vafs <- as.numeric(c(combined$tumor_vaf, combined$normal_vaf))

combined <- rbind(combined, combined)
combined <- combined %>%
	mutate(mut_type = mut_type, # three new cols 
		   mut_vafs = mut_vafs, 
		   variant_id = paste(combined$sample_t, pos, sep = "_")) %>% # this one is so that the plot we have is more informative about the variant ids
		   select(-tumor_vaf, -normal_vaf) # delete

# i am also doing some string modification here to have the patient ids as a separate column. This will help coloring the axis labels later on while plotting.
# i actually ended up not using colors. 
x <- lapply(combined$sample_t, function(sample_id) head(str_split(sample_id, "-")[[1]], -2)) # just str split by - and disregard the last two elements
combined$patient_id <- unlist(lapply(x, function(sample_id) paste(sample_id[1], sample_id[2], sep = "-")))

p <- ggplot(data = combined, aes(x = mut_vafs*100, y = variant_id, color = mut_type)) + 
	geom_point(size = 2) +
	# geom_text_repel(label = round(mut_vafs, 2)) + # can uncomment later to label the dots with the allele freq 
	scale_x_continuous(breaks = seq(0, 100, 10)) + 
	geom_hline(yintercept = seq(0, length(combined$variant_id), 1), color = "gray", linetype = "dotted") +
	geom_vline(xintercept = seq(0, 100, 10), color = "black", linetype = "dotted") +
	xlab("VAF percentage") + 
	ylab("Variants") +
	theme_bw() + 
	theme(legend.position = c(0.95, 0.1)) + 
	cool_theme 

path_to_classic_plot <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/figures/classic_plot.pdf"
ggsave(path_to_classic_plot, p, width = 12, height = 7, units = "in")

##########################################################################
# A PLOT SHOWING ctDNA PERCENTAGE IN SAMPLES WHERE WE IDENTIFIED VARIANTS
##########################################################################
samples_path <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/sample_information.csv" # a copy of the google sheets with sample information
samples <- as.data.frame(read_csv(samples_path)[, c(1, 5, 6)])
samples <- samples[-1, ] # the first row is the colname
names(samples) <- c("sample_id", "percentage_genePanel", "percentage_WGS") # rename the cols
samples$sample_id <- paste0("RG_", samples$sample_id, ".bam")# some string modifications to make sure the names match across data frames

# filter the big samples df to keep the samples we have present in the combined df. Note that not all the samples have a corresponding tumor percentage in the samples df
samples <- samples %>%
	filter(sample_id %in% combined$sample_t)

combined <- inner_join(combined, samples, by = c("sample_t" = "sample_id")) # only keep the samples that are found in both dfs
combined$percentage_WGS <- as.numeric(combined$percentage_WGS)
combined <- combined[order(combined$percentage_WGS),]

p2 <- ggplot(data = combined, aes(x = percentage_WGS, y = variant_id)) + 
	geom_point(size = 2) +
	scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) + 
	geom_hline(yintercept = seq(0, length(combined$variant_id), 1), color = "gray", linetype = "dotted") +
	geom_vline(xintercept = seq(0, 100, 10), color = "black", linetype = "dotted") +
	xlab("ctDNA Percentage (WGS)") + 
	ylab("Variants") +
	theme_bw() + 
	cool_theme

path_to_WGS_percentage <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/figures/tumor_percentage.pdf"
ggsave(path_to_WGS_percentage, p2, width = 12, height = 7, units = "in")

gridExtra::grid.arrange(p, p2, nrow = 1)