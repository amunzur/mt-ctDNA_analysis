library(tidyverse)

# This script generates a histogram of fragment sizes for all samples in a given window size

cool_theme <-

  theme(panel.border = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank(),
	  axis.line = element_line(colour = "black", size = 1),
	  axis.ticks = element_line(colour = "black", size = 2),
	  axis.text = element_text(size=10),
	  axis.text.x = element_text(vjust=0.5, colour = "black", size=12),
	  axis.text.y = element_text(vjust=0.5, colour = "black", size=12),
	  axis.title = element_text(size=14, face="bold"),
	  legend.title = element_text(color = "black", size = 12),
	  legend.text = element_text(color = "black", size = 12),
	  axis.ticks.length=unit(0.15, "cm"),
	  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 20)),
	  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
	  strip.background = element_blank(),
	  strip.text = element_text(size = 12))

# given one bed file, generate a list of all fragments=
generate_fragments <- function(path_to_bed) {

	bed <- as.data.frame(read.delim(path_to_bed, sep = "\t", header = FALSE))
	names(bed) <- c("chr", "start", "end", "readname", "mapq", "strand")

	pair1_start <- bed$start[seq(1, dim(bed)[1], 2)]
	pair2_start <- bed$start[seq(2, dim(bed)[1], 2)]

	pair1_end <- bed$end[seq(1, dim(bed)[1], 2)]
	pair2_end <- bed$end[seq(2, dim(bed)[1], 2)]

	fragment_start <- pmin(pair1_start, pair2_start) # smaller of the two
	fragment_end <- pmax(pair1_end, pair2_end) # larger of the two

	# put in a df
	fragment_df <- data.frame(fr_start = fragment_start, fr_end = fragment_end, fr_size = fragment_end - fragment_start) %>% arrange(fragment_start)

	return(fragment_df)

}

fragment_size_HISTOGRAM <- function(fragment_df){

	p <- ggplot(data = fragment_df, aes(x = fr_size)) + 
		geom_histogram(bins = 40) +
		scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) + 
		scale_y_continuous(breaks = scales::pretty_breaks(n = 15)) + 
		xlab("Fragment size (bp)") +
		ylab("Number of fragments") +
		cool_theme

	return(p)

}

main <- function(dir_to_beds, dir_to_figs){

	for (path_to_bed in list.files(dir_to_beds, full.names = TRUE)) {
		
		fragment_df <- generate_fragments(path_to_bed)
		message(path_to_bed)
		p <- fragment_size_HISTOGRAM(fragment_df)

		# extract sample name from the bed path, will be needed to save the fig
		sample_name <- gsub("RG_filtered_RG_|.bed", "", basename(path_to_bed))
		message(sample_name)
		figure_path <- file.path(dir_to_figs, paste0(sample_name, ".png"))

		ggsave(figure_path, p, height = 15, width = 15)

		} # end of for loop

	} # end of function

dir_to_beds <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/data/beds"
dir_to_figs <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/figures/fragment_size_histograms"
main(dir_to_beds, dir_to_figs)