library(tidyverse)
# This script takes in a WPS df with the peak information and visualized the peaks

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


return_peak_df <- function(path_to_WPS){
	# Given a WPS df, return a df with the idx of peaks in the format ready to plot a histogram
	df <- read.delim(path_to_WPS, sep = "\t")

	df$peak <- ifelse(df$peak == "True", TRUE, FALSE) # from python syntax to R syntax for boolean
	idx <- which(df$peak) # idx of the peaks

	df_peak <- data.frame(position = idx, 
						dummy_value = rep(0, length(idx)), 
						sample_name = gsub(".tsv", "", basename(path_to_WPS)))
	return(df_peak)

}

visualize_peaks_HISTOGRAM <- function(df_peak, dir_to_figs){

	mat <- as.matrix(dist(df_peak$position))
	peak_diff <- mat[row(mat) == col(mat) + 1] # distance between all peaks

	df <- data.frame(peak_diff = peak_diff, dummy_value = 0)
	window_size <- gsub("_", " ", grep("window_size", strsplit(path_to_WPS, "/")[[1]], value = TRUE))
	message(unique(df_peak$sample_name)) # print to console to track progress

	p <- ggplot(data = df, aes(x = peak_diff)) + 
		geom_histogram(bins = 40) + 
		xlab("Distance between peaks (bp)") +
		ylab("Counts") + 
		scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) + 
		scale_y_continuous(breaks = scales::pretty_breaks(n = 15)) + 
		ggtitle(paste(unique(df_peak$sample_name), window_size, sep = " / ")) +
		cool_theme

	return(p)

}

main <- function(dir_to_WPS, dir_to_figs, window_size_vector){
	# Given proper params, call the functions above on all samples where WPS df is available

	for (window_size in window_size_vector){
		message("WINDOW SIZE ", window_size)

		for (path in list.files(file.path(dir_to_WPS, paste0("window_size_", window_size)), full.names = TRUE)){
		
			df_peak <- return_peak_df(path)
			p <- visualize_peaks_HISTOGRAM(df_peak, dir_to_figs)

			figure_path <- file.path(dir_to_figs, paste0("window_size_", window_size), paste0(unique(df_peak$sample_name), ".png"))

			ggsave(figure_path, p, height = 15, width = 15)

		} # end of for loop for window sizes

	} # end of for loop for paths 

} # end of function

dir_to_WPS <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/data/WPS_files"
dir_to_figs <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/figures/histograms"
main(dir_to_WPS, dir_to_figs, c(15, 30, 60, 120))