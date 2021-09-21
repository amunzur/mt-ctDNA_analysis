library(tidyverse)
library(R.utils)
######################################################################
# SUBSET AFTER MANUAL IGV CURATION
######################################################################

DIR_to_curated_snapshots <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/snapshot/VarScan2/snv_curated"
PATH_to_anno_df <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/annovar/annovar_results/VarScan2/snv.hg38_multianno.txt"
output_dir <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarScan2/realigned_mito_filtered/curated_dfs"

subset_to_snapshots <- function(DIR_to_curated_snapshots, PATH_to_anno_df, output_dir){

	bams <- unlist(sapply(str_split(list.files(DIR_to_curated_snapshots, pattern = "filtered"), "_chrM_"), function(x) x[1]))
	positions <- unlist(sapply(str_split(list.files(DIR_to_curated_snapshots, pattern = "filtered"), "_chrM_"), function(x) strsplit(x[2], ".png")[1]))

	df <- data.frame(Tumor = bams, position = positions) # df with variants we want 
	anno_df <- as.data.frame(read.delim(PATH_to_anno_df))
	curated <- inner_join(anno_df, df)

	write_delim(curated, file.path(output_dir, "snv.txt"), delim = "\t")
	write_csv(curated, file.path(output_dir, "snv.csv"))
	write_delim(curated_distinct, file.path(output_dir, "snv_distinct.txt"), delim = "\t")
	write_csv(curated_distinct, file.path(output_dir, "snv_distinct.csv"))

	return(curated)

}

# VarScan2 calls
curated <- subset_to_snapshots(DIR_to_curated_snapshots, PATH_to_anno_df, output_dir)
curated_distinct <- distinct(curated, position, .keep_all = TRUE) # only unique variants

######################################################################
# ANNOTATIONS
#####################################################################
# CLINVAR
clinvar <- as.data.frame(read.delim("/groups/wyattgrp/software/annovar/annovar/humandb/hg38_clinvar_20170130_MT.txt"), header = FALSE))
names(clinvar)[c(1, 2, 3, 4, 5, 6)] <- c("chrom", "position", "stop", "ref", "var", "pathogenecity")
curated_distinct[c("position", "ref", "var")] <- lapply(curated_distinct[c("position", "ref", "var")], as.character)
clinvar[c("position", "ref", "var")] <- lapply(clinvar[c("position", "ref", "var")], as.character)
curated_distinct <- left_join(curated_distinct, clinvar, by = c("position" = "position", "ref" = "ref", "var" = "var"))

# YUAN YUAN 2020 STUDY
yuan2020 <- as.data.frame(read.delim("/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/database/TCMA-MutationSNV.tsv"))
yuan2020[c("position", "ref", "var")] <- lapply(yuan2020[c("position", "ref", "var")], as.character)
curated2 <- inner_join(curated_distinct, yuan2020, by = c("position" = "position", "ref" = "ref", "var" = "var"))

# HELIX DB 
helix <- as.data.frame(read.delim("/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/database/HelixMTdb_20200327.tsv", header = TRUE))
helix$alleles <- gsub("\\[|\\]", "", helix$alleles) # remove square brackets 

helix$ref <- unlist(lapply(helix$alleles, function(x) strsplit(x, ",")[[1]][[1]])) # extract the ref alleles 
helix$var <- as.vector(lapply(helix$alleles, function(x) strsplit(x, ",")[[1]][-1])) # extract the ref alleles 

# disregard indels for now
indel_locations <- which(sapply(var, function(x) length(x)) == 2)
helix <- helix[-indel_locations, ]
helix <- helix %>% separate(locus, c("chrom", "position")) 
helix$var <- unlist(helix$var)

# refs <- insert(helix$ref, indel_locations + 1, helix$ref[indel_locations]) # repeated refs inserted into the right position
curated_helix <- inner_join(curated_distinct, helix, by = c("chrom" = "chrom", "position" = "position", "ref" = "ref", "var" = "var")) %>%
						select(position, ref, var, samplename, AF_hom, AF_het, mean_ARF, max_ARF, counts_hom, counts_het) %>%
						arrange(as.numeric(position)) 

helix$position	<- factor(as.numeric(helix$position), levels = unique(c(helix$position, curated_distinct$position)))

# log transformation
curated_helix$counts_hom[which(curated_helix$counts_hom == 0)] <- NA
curated_helix$counts_het[which(curated_helix$counts_het == 0)] <- NA

p <- curated_helix %>% 
		gather("variant_type", "counts", counts_hom, counts_het) %>%
		ggplot(aes(x = factor(position, level = unique(sort(as.numeric(as.character(position))))), y = log(counts), fill = variant_type)) + 
		geom_bar(position="dodge", stat="identity") + 
		scale_fill_manual(values = c("#F2BC94", "#1F8AC0")) + 
		xlab("Genomic location of the SNV") + 
		ylab("Number of individuals, log transformed")

ggsave("/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/VarScan2/helixdb_fig.pdf", p, width = 15, height = 5)


make_plasmicity_type_plot <- function(df, output_dir) {

	df_plot <- data.frame(both = length(which(as.numeric(as.logical(df$counts_hom)) + as.numeric(as.logical(df$counts_het)) == 2)), 
						hom_only = length(which(as.numeric(as.logical(df$counts_hom)) - as.numeric(as.logical(df$counts_het)) == -1)), 
						het_only = length(which(as.numeric(as.logical(df$counts_het)) - as.numeric(as.logical(df$counts_hom)) == -1)))

	counts_hom <- df_plot$counts_hom[is.na(df_plot$counts_hom)] <- 0





	df_plot$counts_hom[is.na(df_plot$counts_het)] <- 0

	p <- data.frame(both = length(which(as.numeric(as.logical(df$counts_hom)) + as.numeric(as.logical(df$counts_het)) == 2)), 
					hom_only = length(which(as.numeric(as.logical(df$counts_hom)) - as.numeric(as.logical(df$counts_het)) == -1)), 
					het_only = length(which(as.numeric(as.logical(df$counts_het)) - as.numeric(as.logical(df$counts_hom)) == -1))) %>%
					gather("plasmicity", "freq", both, hom_only, het_only) %>% 
					ggplot(aes(x = plasmicity, y = freq)) + geom_bar(stat="identity")

	ggsave(file.path(output_dir, "plasmicity_fig.pdf"), p, width = 5, height = 5)}

make_lollipop_plot <- function(helix, curated_distinct, output_dir) {

	# two extra cols needed for plotting 
	helix <- helix %>% 
			mutate(point_size = 10,
					point_shape = 0,
					point_color = ifelse(helix$feature == "non_coding", "chartreuse3", 
									ifelse(helix$feature == "protein_coding_gene", "brown2", 
									ifelse(helix$feature == "rRNA_gene", "dodgerblue2", 
									ifelse(helix$feature == "tRNA_gene", "goldenrod2", "black")))),
					tumor_var_freq = 0)

	helix_slice <- helix[1:dim(curated_distinct)[1], ] %>%
					mutate(point_size = 3, 
							point_shape = 16, 
							point_color = "gray",
							tumor_var_freq = curated_distinct$tumor_var_freq, 
							position = curated_distinct$position)
	
	helix <- rbind(helix, helix_slice)

	p <- ggplot(data = helix, aes(x = position, y = tumor_var_freq)) + 
		geom_segment(aes(x = position, xend = position, y = 0, yend = tumor_var_freq), color = "gray") +
		geom_point(shape = helix$point_shape, size = helix$point_size, color = helix$point_color) + 
		theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

	ggsave(file.path(output_dir, "lollipop.pdf"), p, width = 10, height = 5)}

make_coverage_plot <- function(curated, DIR_to_tnvstats, output_dir) {

	tnvstats_paths <- file.path(rep(DIR_to_tnvstats, 
		dim(curated)[1]), 
		gsub(".vcf", ".bam", curated$samplename), 
		paste0(gsub(".vcf", ".bam", curated$samplename), ".tnvstat"))

	tnvstats_all <- mapply(function(x, idx) {list(x[which(x$pos == idx), ])}, x = lapply(tnvstats_paths, read.delim), idx = curated$position) # x is just a list of all tnvstats read into a list
	tnvstats <- as.data.frame(do.call(rbind, tnvstats_all)) %>% 
				select(pos, reads_all_n, reads_all_t, sample_n, sample_t) %>%
				mutate(sample_name = gsub("filtered_RG_|.bam", "", sample_t)) %>%
				unite("pos_id", pos, sample_name, sep = "_", remove = FALSE) %>%
				gather("sample_type", "depth", reads_all_n, reads_all_t) %>%
				arrange(pos)

	tnvstats$pos_id <- factor(tnvstats$pos_id, level = unique(tnvstats$pos_id))
	tnvstats$sample_type <- ifelse(tnvstats$sample_type == "reads_all_n", "WBC", "Tumor")
	rownames(tnvstats) <- NULL

	p <- ggplot(data = tnvstats, aes(x = pos_id, y = depth, fill = sample_type)) + 
		geom_bar(position="dodge", stat="identity") + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
		coord_flip()

	ggsave(file.path(output_dir, "coverage.pdf"), p, width = 10, height = 20)}
}

############################################
# VARSCAN2
############################################
DIR_to_curated_snapshots <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/snapshot/VarScan2/snv_curated"
PATH_to_anno_df <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/annovar/annovar_results/VarScan2/snv.hg38_multianno.txt"
output_dir <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarScan2/realigned_mito_filtered/curated_dfs"

fig_dir <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/VarScan2"
DIR_to_tnvstats <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/tnvstats_realigned_NEW"

curated <- subset_to_snapshots(DIR_to_curated_snapshots, PATH_to_anno_df, output_dir)
curated_distinct <- distinct(curated, position, .keep_all = TRUE) # only unique variants

make_plasmicity_type_plot(curated_helix, fig_dir)
make_lollipop_plot(helix, curated_distinct, fig_dir)
make_coverage_plot(curated, DIR_to_tnvstats, fig_dir)


############################################
# MUTATO
############################################
DIR_to_curated_snapshots <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results_realigned/snapshots_curated"
PATH_to_anno_df <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/annovar/annovar_results/mutato/ALL_MERGED_FILTERED0.01.hg38_multianno.txt"
output_dir <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results_realigned/curated_dfs"

fig_dir <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/VarScan2"
DIR_to_tnvstats <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/tnvstats_realigned_NEW"

curated <- subset_to_snapshots(DIR_to_curated_snapshots, PATH_to_anno_df, output_dir)
curated_distinct <- distinct(curated, position, .keep_all = TRUE) # only unique variants

make_plasmicity_type_plot(curated_helix, fig_dir)
make_lollipop_plot(helix, curated_distinct, fig_dir)
make_coverage_plot(curated, DIR_to_tnvstats, fig_dir)



