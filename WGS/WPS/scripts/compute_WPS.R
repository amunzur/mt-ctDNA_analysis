library(tidyverse)
library(argparse)

parser <- ArgumentParser(description = "Calculate window protection score bed files.")

parser$add_argument('--path-to-bed', metavar='FILE', type='character', help="Path to the bed file from one single sample.")
parser$add_argument('--window-size')

args <- parser$parse_args()

# given a genome size in bp, generate a list of all windows with a given size. start and end pos given.
generate_windows <- function(genome_length, window_size) {

	window_size <- as.numeric(window_size)

	mid_pos <- seq(0, genome_length-1, 1) # because bed is 0 based
	start_pos <- mid_pos - window_size/2
	end_pos <- mid_pos + window_size/2

	# put in a df
	window_df <- data.frame(window_start = start_pos, window_mid = mid_pos, window_end = end_pos)

	return(window_df)
}

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

# given windows and fragments, compute WPS for each window 
compute_WPS <- function(window_df, fragment_df, output_dir){

	compute_WPS_per_window <- function(x, fragment_df){

		if (x$window_start %% 1000 == 0) {message("Position ", x$window_start)}
		
		# right end point in window
		right_end <- dim(fragment_df %>% filter(fr_start <= x$window_start & 
												fr_end <= x$window_end & 
												fr_end >= x$window_start))[1]
		# left end point in window
		left_end <- dim(fragment_df %>% filter(fr_start >= x$window_start & 
												fr_start <= x$window_end & 
												fr_end >= x$window_start))[1]

		number_ends <- right_end + left_end

		# fragments in the window
		number_in_window <- dim(fragment_df %>% filter(fr_start < x$window_start & 
											fr_end > x$window_end))[1]

		df <- data.frame(right_end_pos = right_end, 
							left_end_pos = left_end, 
							number_ends = number_ends, 
							number_in_window = number_in_window)
		return(df)

	}

	# convert window df to a list of one line dfs 
	window_list <- lapply(as.list(1:dim(window_df)[1]), function(x) window_df[x[1],])
	landscape <- do.call(rbind, lapply(window_list, compute_WPS_per_window, fragment_df = fragment_df)) # computing WPS for all windows
	WPS_df <- cbind(window_df, landscape)
	WPS_df$WPS <- WPS_df$number_in_window - WPS_df$number_ends

	# remove the positions outside of the range of mito
	WPS_df <- WPS_df %>% filter(window_start >= 0, 
								window_end <= 19569)

	write_delim(WPS_df, output_dir, delim = "\t")

	return(WPS_df)
}

# given a path to a bed file extract the sample name and determine the txt file with the cov info
add_coverage_info <- function(path_to_bed, dir_to_coverage, WPS_df){

	path_to_coverage <- file.path(dir_to_coverage, gsub(".bed", ".txt", basename(path_to_bed)))
	coverage_df <- read.delim(path_to_coverage, header = FALSE)
	names(coverage_df) <- c("chr", "position", "depth")
	
	# because WPS df is 0 based
	coverage_df$position <- coverage_df$position - 1
	WPS_df <- merge(WPS_df, coverage_df, by.x = "window_start", by.y = "position")

	return(WPS_df)
}

plotting <- function(WPS_df, fig_path){

WPS_df <- gather(WPS_df, key = "K", value = "V", WPS, depth)
WPS_df$K <- ifelse(WPS_df$K == "WPS", "Window protection score", "Coverage")
	
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

  	p <- ggplot(data = WPS_df, aes(x = window_start, y = V)) + 
		facet_wrap(~ K, ncol = 1) +
		geom_bar(stat = 'identity') + 
		scale_x_continuous(seq(0, 16000, 2000), name = "mitochondrial genome") + 
		geom_vline(xintercept = seq(0, 16000, 2000), color = "gray", linetype = "dotted") +
		cool_theme

	ggsave(fig_path, p, width = 30, height = 10)

}

main <- function(window_size, path_to_bed){

	main_dir <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS"

	# input
	dir_to_beds <- file.path(main_dir, "WPS/data/beds")
	dir_to_coverage <- file.path(main_dir, "metrics/read_depth/read_depth_realigned_paired")
	# output
	dir_to_figs <- file.path(main_dir, "WPS/figures", paste0("window_size_", window_size))
	dir_to_WPS_files <- file.path(main_dir, "WPS/data/WPS_files", paste0("window_size_", window_size))

	# create dirs for the window sizes we chose 
	dir.create(dir_to_figs, showWarnings = FALSE)
	dir.create(dir_to_WPS_files, showWarnings = FALSE)

	sample_name <- gsub(".bed|RG_filtered_RG_", "", basename(path_to_bed)) # will be used to name the plot
	message("Started sample ", sample_name, ".")

	window_df <- generate_windows(16569, window_size)
	fragment_df <- generate_fragments(path_to_bed)

	message("Computing WPS.")
	output_dir <- file.path(dir_to_WPS_files, paste0(sample_name, ".tsv")) # WPS df will be saved here
	WPS_df <- compute_WPS(window_df, fragment_df, output_dir)
	WPS_df <- add_coverage_info(path_to_bed, dir_to_coverage, WPS_df)

	message("Plotting.")
	fig_path <- file.path(dir_to_figs, paste0(sample_name, ".png"))
	plotting(WPS_df, fig_path)

	} # end of function

# CALL THE FUNCTION WITH THE CM ARGS
main(window_size = args$window_size,
	path_to_bed = args$path_to_bed)