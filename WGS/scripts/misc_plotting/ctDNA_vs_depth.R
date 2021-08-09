library(tidyverse)

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
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

#################################################################
# WITHOUT NUMTS 
#################################################################

path_to_df <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth_realigned/mito_filtered_realigned.csv"
path_to_sample_information <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/sample_information.csv"
fig_path <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/depth_and_tumor_fraction/main_fig.pdf"

df <- as.data.frame(read.csv(path_to_df, header = FALSE)) # load the sample df data 
names(df) <- c("sample_path", "depth")
tumor_fr <- as.data.frame(read_csv(path_to_sample_information)) %>% select("X1", "X6") # load the sample df data 

# some clean up in the df file so that the file names match across files 
df$sample_names <- unlist(lapply(basename(as.character(df[, 1])), function(x) strsplit(x, "filtered_RG_")[[1]][[2]]))
df$sample_names <- unlist(lapply(df$sample_names, function(x) strsplit(x, ".bam")[[1]][[1]]))

df <- left_join(df, tumor_fr, by = c("sample_names" = "X1"))
names(df) <- c("sample_path", "depth", "sample_name", "tumor_fr")

df <- df[!is.na(df[, 4]), ]# drop na, they are the wbc samples
df$tumor_fr <- as.numeric(df$tumor_fr)
df$depth <- as.numeric(df$depth)

p <- ggplot(data = df, aes(x = depth, y = tumor_fr)) + 
	geom_point() + 
	geom_smooth(method = "lm", se = FALSE, color = "red") + 
	scale_x_continuous(breaks = seq(0, 6000, by = 1000), limits = c(0, 6000)) + 
	scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) + 
	ylab("WGS tumor %") + 
	xlab("mt-ctDNA depth without NUMTs") + 
	geom_hline(yintercept = seq(0, 100, 20), linetype = "dotted", color = "gray") + 
	geom_vline(xintercept = seq(0, 6000, 1000), linetype = "dotted", color = "gray") + 
	theme_bw() + 
	cool_theme

ggsave(filename = fig_path, p, width = 10, height = 10, units = "in", device = "pdf")

#################################################################
# WITH NUMTS 
#################################################################

# path_to_df <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth/AVERAGE_DEPTH_PER_SAMPLE"
# path_to_sample_information <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/sample_information.csv"

# df <- as.data.frame(read.delim(path_to_df, header = FALSE)) # load the sample df data 
# names(df) <- c("sample_names", "depth")

# # some clean up in the df file so that the file names match across files 
# df$sample_names <- unlist(lapply(as.character(df$sample_names), function(x) strsplit(x, "RG_")[[1]][[2]]))
# df$sample_names <- unlist(lapply(as.character(df$sample_names), function(x) strsplit(x, ".bam")[[1]][[1]]))

# df <- left_join(df, tumor_fr, by = c("sample_names" = "X1"))
# names(df) <- c("sample_names", "depth", "tumor_fr")

# df <- df[!is.na(df[, 3]), ]# drop na, they are the wbc samples
# df$tumor_fr <- as.numeric(df$tumor_fr)
# df$depth <- as.numeric(df$depth)

# p2 <- ggplot(data = df, aes(x = depth, y = tumor_fr)) + 
# 	geom_point() + 
# 	geom_smooth(method = "lm", se = FALSE, color = "red") + 
# 	scale_x_continuous(breaks = seq(0, 6000, by = 1000), limits = c(0, 6000)) + 
# 	scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) + 
# 	ylab("WGS tumor %") + 
# 	xlab("mt-ctDNA depth without NUMTs") + 
# 	geom_hline(yintercept = seq(0, 100, 20), linetype = "dotted", color = "gray") + 
# 	geom_vline(xintercept = seq(0, 6000, 1000), linetype = "dotted", color = "gray") + 
# 	theme_bw() + 
# 	cool_theme


