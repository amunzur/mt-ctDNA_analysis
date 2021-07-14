library(tidyverse)

# start with the average read depth plots
path_to_combined <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth/AVERAGE_DEPTH_PER_SAMPLE"
path_to_variants <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/results/somatic_cleaned.csv"
dir_to_readdepths <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth"

df <- read.delim(path_to_combined, header = FALSE) # main file we'll be working on 
names(df) <- c("bam_name", "depth")
df$depth <- round(df$depth, 0) # who likes decimals anyways
df <- df[order(df$depth),] # order by depth

var_df <- read_csv(path_to_variants) # matti variant results 
var_df$sample_names <- paste0(var_df$sample_names, ".bam") # just so that the names match across files 

df$variant_identified <- as.integer(df$bam_name %in% var_df$sample_names) # 1 means the sample has a variant, 0 means it doesn't.
df$colors <- ifelse(df$variant_identified == 1, "red", "black") # helps coloring the plots later on

# Now look into the depth at exact position where the snvs occurred. For this, we need to do a bit of a file hunting: identify the 
# sample with the variant, go to the correct file, and pull the read depth at that position. 
i <- 1 
snv_depth_list <- list() # depths at snv
sample_depth_list <- list() # average sample depths, taken from df
while (i <= dim(var_df)[1]){ # iterate through the variants
	i_row <- var_df[i, ]
	
	# INDIVIDUAL DEPTHS
	df_individual_depth <- as.data.frame(read.delim(file.path(dir_to_readdepths, i_row$sample_names), header = FALSE)) # read the depth file from the sample of interest
	names(df_individual_depth) <- c("CHROM", "POSITION", "DEPTH") # we are interested in the 3rd col 
	
	df_individual_depth <- filter(df_individual_depth, POSITION == i_row$POSITION) # filter to keep the snv pos only
	depth <- df_individual_depth$DEPTH # extract the depth at the snv
	snv_depth_list <- append(depth, snv_depth_list)

	# AVERAGE DEPTHS 
	d <- as.numeric(filter(df, bam_name == i_row$sample_names)[2])
	sample_depth_list <- append(sample_depth_list, d)

	i <- i + 1
}

var_df$SNV_DEPTH <- unlist(snv_depth_list)
var_df$VARIANT_ID <- paste(var_df$sample_names, var_df$POSITION, sep = "_") # modify the names slightly so that they match across files
var_df$AVERAGE_DEPTH <- unlist(sample_depth_list)

# we might as well save this now that we did some useful modifications 
write_csv(var_df, "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/matti_pipeline_results/results/somatic_cleaned_depth.csv")

##################################
# PLOTTING
##################################
cool_theme <- 
  
  theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black", size = 1), 
      axis.ticks = element_line(colour = "black", size = 2),
      axis.text = element_text(size=10),
      axis.text.x = element_text(vjust=0.5, colour = "black", size=8),
      axis.text.y = element_text(vjust=0.5, colour = "black", size=6),
      axis.title = element_text(size=10,face="bold"), 
      legend.title = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 12),
      axis.ticks.length=unit(0.15, "cm"), 
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 20)), 
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

# A summary of plots below: 
# p1: all samples, showing average sample depth as a dot plot. Samples with variants are indicated in red.
# p2: all samples, showing average sample depth as a histogram
# p3: samples with variants, average sample depth as a dot plot
# p4: samples with variants, average sample depth as a histogram
# p5: exact read depth at the positions snvs are identified. 

# plot showing depth in all individual samples
p1 <- ggplot(data = df, aes(x = depth, y = reorder(bam_name, -depth))) + 
	geom_point(color = df$colors) + 
	geom_hline(yintercept = seq(1, length(df$bam_name), 2), color = "gray", linetype="dotted") +
	geom_hline(yintercept = seq(0, length(df$bam_name), 2), color = "black", linetype="dotted") +
	xlab("Average read depth") + 
	ylab(" ") +
	ggtitle("Average read depth in WGS samples", subtitle = "Samples with variant are shown in red.") +
	scale_x_continuous(breaks = seq(0, 8000, 1000), limits = c(0, 8000)) + 
	theme_bw() + 
	cool_theme + 
	theme(axis.text.y = element_text(colour = df$colors), 
		axis.ticks = element_line(colour = "black", size = 1))

# a histogram, all samples
bin_number <- 20
p2 <- ggplot(data = df, aes(x = depth)) + 
	geom_hline(yintercept = seq(1, 11, 1), color = "black", linetype="dotted") + # horizontal grid lines
	geom_histogram(bins=bin_number, fill="lightblue") + 
	geom_vline(aes(xintercept=mean(depth)), color="blue", linetype="dashed", size=1) +
	# geom_vline(aes(xintercept=median(depth)), color="red", linetype="dashed", size=1) +
	xlab("Average read depth") + 
	ylab("Number of samples") +
	ggtitle("Histogram of average read depth in WGS samples", subtitle = "Dashed line indicates mean.") +
	scale_x_continuous(breaks = seq(0, 8000, 800)) + 
	scale_y_continuous(breaks = seq(0, 12, 1)) +
	theme_bw() + 
	cool_theme + 
	theme(axis.text.y = element_text(vjust=0.5, colour = "black", size=8))

# This time a plot for the samples with variants only. First from wide to long format. 
# depth_values: the name of the new col with the values 
# depth_types: variable names stored here 
# 3rd one is the cols to use
# var_df <- as.data.frame(gather(var_df, depth_types, depth_values, c("SNV_DEPTH", "AVERAGE_DEPTH")))
# var_df$colors <- ifelse(var_df$depth_types == "SNV_DEPTH", "#1B9E77", "#D95F02")

p3 <- df %>% filter(variant_identified == 1) %>%
	ggplot(aes(x = depth, y = bam_name)) + 
	geom_point() + 
	xlab("Average read depth") + 
	ylab(" ") +
	ggtitle("Average read depth in WGS samples with variants") +
	scale_x_continuous(breaks = seq(0, 8000, 1000), limits = c(0, 8000)) + 
	geom_hline(yintercept = seq(1, length(df$bam_name), 2), color = "gray", linetype="dotted") +
	geom_hline(yintercept = seq(0, length(df$bam_name), 2), color = "black", linetype="dotted") +
	theme_bw() + 
	cool_theme + 
	theme(axis.text.x = element_text(vjust=0.5, colour = "black", size=10),
		axis.text.y = element_text(vjust=0.5, colour = "black", size=10))

# same as p3 but now we are making a histogram
p4 <- df %>% filter(variant_identified == 1) %>%
	ggplot(aes(x = depth)) + 
	geom_histogram(bins=16, fill="lightblue") + 
	xlab("Average read depth") + 
	ylab("Number of samples") +
	ggtitle("Histogram of average read depth in WGS samples with variants", subtitle = "Dashed line indicates mean.") +
	scale_x_continuous(breaks = seq(0, 8000, 1000), limits = c(0, 8000)) + 
	scale_y_continuous(breaks = seq(0, 7, 1), limits = c(0, 7)) + 
	geom_hline(yintercept = seq(0, 7, 1), color = "black", linetype="dotted") + # horizontal grid lines
	geom_vline(aes(xintercept=mean(depth)), color="blue", linetype="dashed", size=1) + # the mean
	theme_bw() + 
	cool_theme + 
	theme(axis.text.x = element_text(vjust=0.5, colour = "black", size=10),
		axis.text.y = element_text(vjust=0.5, colour = "black", size=10))

# showing the depth of all the snvs
p5 <- ggplot(data = var_df, aes(x = DEPTH, y = VARIANT_ID)) + 
	geom_point() + 
	geom_hline(yintercept = seq(1, length(df$bam_name), 2), color = "gray", linetype="dotted") +
	geom_hline(yintercept = seq(0, length(df$bam_name), 2), color = "black", linetype="dotted") +
	xlab("Average read depth") + 
	ylab(" ") +
	ggtitle("Read depth in WGS samples at each SNV") +
	scale_x_continuous(breaks = seq(0, 8000, 1000), limits = c(0, 8000)) + 
	theme_bw() + 
	cool_theme + 
	theme(axis.text.x = element_text(vjust=0.5, colour = "black", size=10),
		axis.text.y = element_text(vjust=0.5, colour = "black", size=10))

# all paths are below 
p1_path <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/sample_depth/average_depth_dot_ALL.pdf"
p2_path <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/sample_depth/average_depth_hist_ALL.pdf"
p3_path <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/sample_depth/average_depth_dot_VARIANTS.pdf"
p4_path <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/sample_depth/average_depth_hist_VARIANTS.pdf"
p5_path <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/sample_depth/snv_depth_VARIANTS.pdf"

ggsave(p1_path, p1, width = 12, height = 12, units = "in")
ggsave(p2_path, p2, width = 12, height = 7, units = "in")
ggsave(p3_path, p3, width = 12, height = 7, units = "in")
ggsave(p4_path, p4, width = 12, height = 7, units = "in")
ggsave(p5_path, p5, width = 12, height = 7, units = "in")

