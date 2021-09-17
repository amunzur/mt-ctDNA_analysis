# this script is for calculating CN in samples and plotting. 
library(tidyverse)

path_to_depth_df <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth/read_depth_realigned_median/mito_filtered_realigned.csv"
path_to_depth_df_mito_mean <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/read_depth/read_depth_realigned_mean/mito_filtered_realigned.csv"
path_to_depth_df_WGS <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/coverage_hist_WGS/average_coverage.csv" # average coverage per sample across all positions 
path_to_sample_information <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/sample_information.csv" # google sheet

#####################################################
# COMPUTING CN 
#####################################################
compute_CN <- function(path_to_depth_df, 
                        path_to_depth_df_mito_mean, 
                        path_to_depth_df_WGS, 
                        path_to_sample_information, 
                        WBC) {

    depth_df <- as.data.frame(read.csv(path_to_depth_df, header = FALSE)) # load the sample df data 
    sample_info_df <- as.data.frame(read_csv(path_to_sample_information)) %>% select(X1, X6, X8, X9, X16) %>% mutate(across(everything(), as.character))
    # sample name, cancer % WGS, base ploidy, genome duplication, median coverage 

    # some clean up in the depth_df file so that the file names match across files 
    depth_df$sample_names <- unlist(lapply(basename(as.character(depth_df[, 1])), function(x) strsplit(x, "filtered_RG_")[[1]][[2]]))
    depth_df$sample_names <- unlist(lapply(depth_df$sample_names, function(x) strsplit(x, ".bam")[[1]][[1]]))

    if (WBC == FALSE) {depth_df <- depth_df[-c(grep("WBC", depth_df$sample_names)), ]
    depth_df <- depth_df[-c(grep("DTB-119-Progression-cfDNA|AE-018-Baseline-cfDNA", depth_df$sample_names)), ] # remove some problematic samples

        } else {depth_df <- depth_df[c(grep("WBC", depth_df$sample_names)), ]} # remove all WBC samples 
    
    # add info from sample sheet
    depth_df <- left_join(depth_df, sample_info_df, by = c("sample_names" = "X1"))
    names(depth_df) <- c("sample_path", "mtDNA_median_depth", "sample_name", "tumor_perc", "genome_duplication", "base_ploidy", "gDNA_median_depth")
    depth_df$tumor_frac <- (as.numeric(depth_df$tumor_perc))/100
    depth_df <- depth_df %>% mutate(across(c("mtDNA_median_depth", "tumor_perc", "base_ploidy", "gDNA_median_depth"), as.numeric))

    # add the average WGS sample depth computed by me
    depth_df_WGS <- as.data.frame(read_csv(path_to_depth_df_WGS))
    depth_df <- left_join(depth_df, depth_df_WGS, by = c("sample_name" = "sample_name")) # adds a new column: coverage

    # add the mean mito depth as well, so far we have the median
    depth_mito_mean <- as.data.frame(read.csv(path_to_depth_df_mito_mean, header = FALSE))
    depth_mito_mean$sample_name <- strsplit(basename(as.vector(depth_mito_mean[, 1])), "filtered_RG_")[[1]][[2]]

    # add a color column to indicate WGD 
    depth_df$genome_duplication_color <- ifelse(depth_df$genome_duplication == "Yes", "chartreuse3", "Red")

    # compute copy number, add as a new col to the df
    depth_df <- depth_df %>% mutate(mito_CN = (mtDNA_median_depth / gDNA_median_depth)*(tumor_frac*base_ploidy + (1 - tumor_frac)*2))
    depth_df <- depth_df %>% mutate(mito_CN_simple = log2(mtDNA_median_depth/median(na.omit(depth_df$mtDNA_median_depth)) / 
                                                           (gDNA_median_depth/median(na.omit(depth_df$gDNA_median_depth)))))

    # drop any cols that have NA only 
    depth_df <- depth_df[, colSums(is.na(depth_df)) != nrow(depth_df)]

    return(depth_df)

}

depth_df_tumor <- compute_CN(path_to_depth_df, 
                        path_to_depth_df_mito_mean, 
                        path_to_depth_df_WGS, 
                        path_to_sample_information, 
                        WBC = FALSE)


depth_df_WBC <- compute_CN(path_to_depth_df, 
                        path_to_depth_df_mito_mean, 
                        path_to_depth_df_WGS, 
                        path_to_sample_information, 
                        WBC = TRUE)
# WBC SAMPLES 
# WBC_depth_df <- as.data.frame(read.csv(path_to_depth_df, header = FALSE)) # load the sample df data 
# sample_info_df <- as.data.frame(read_csv(path_to_sample_information)) %>% select(X1, X6, X8, X9, X16) %>% mutate(across(everything(), as.character))
# sample name, cancer % WGS, base ploidy, genome duplication, median coverage 

# some clean up in the WBC_depth_df file so that the file names match across files 
# WBC_depth_df$sample_names <- unlist(lapply(basename(as.character(WBC_depth_df[, 1])), function(x) strsplit(x, "filtered_RG_")[[1]][[2]]))
# WBC_depth_df$sample_names <- unlist(lapply(WBC_depth_df$sample_names, function(x) strsplit(x, ".bam")[[1]][[1]]))

#####################################################
# PLOTTING
#####################################################

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

p1 <- ggplot(data = depth_df, aes(x = mito_CN, y = tumor_perc)) + 
   geom_point(color = depth_df$genome_duplication_color) +
   geom_smooth(method = "lm", se = FALSE, color = "black") + 
   scale_x_continuous(breaks = seq(0, 150, by = 25), limits = c(0, 155)) + 
   scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 110)) + 
   ylab("WGS tumor %") + 
   xlab("Mito copy number") + 
   geom_hline(yintercept = seq(0, 150, 25), linetype = "dotted", color = "gray") + 
   geom_vline(xintercept = seq(0, 150, 25), linetype = "dotted", color = "gray") + 
   theme_bw() + 
   cool_theme

# simple plots to explore CN for both tumor and wbc 
tumor_CN_log <- ggplot(data = depth_df_tumor, aes(mito_CN_simple)) + 
                  geom_histogram() + 
                  cool_theme + 
                  xlab("mito_CN") + 
                  ggtitle("Tumor samples") + 
                  geom_vline(xintercept = )

WBC_CN_log <- ggplot(data = depth_df_WBC, aes(mito_CN_simple)) + 
                  geom_histogram() + 
                  cool_theme + 
                  xlab("mito_CN") + 
                  ggtitle("WBC samples") + 



path_tumor_CN <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/CN/tumor_CN_log.pdf"
path_WBC_CN <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures/CN/WBC_CN_log.pdf"

# ggsave(path_tumor_CN, tumor_CN_log, height = 20, width = 20, units = "cm")
# ggsave(path_WBC_CN, WBC_CN_log, height = 20, width = 20, units = "cm")

# SHOW WBC AND TUMOR TOGETHER
depth_df_WBC$sample_ID <- gsub("-WBC", "", depth_df_WBC$sample_name)
depth_df_tumor$sample_ID <- gsub("-Base.*|-Prog.*|-cfDNA.*", "", depth_df_tumor$sample_name)

combined <- inner_join(depth_df_tumor, depth_df_WBC, by = "sample_ID") %>%
            select(sample_ID, tumor_perc, sample_name.x, mtDNA_median_depth.x, gDNA_median_depth.x, mito_CN_simple.x, 
                                            sample_name.y, mtDNA_median_depth.y, gDNA_median_depth.y, mito_CN_simple.y)
            
#####################################################################################
# COMPARE WBC AND TUMOR CN
#####################################################################################
df_tumor <- combined %>% 
            select(sample_ID, sample_name.x, mito_CN_simple.x) %>%
            mutate(sample_type = "Tumor")

df_WBC <- combined %>% 
            select(sample_ID, sample_name.y, mito_CN_simple.y) %>%
            mutate(sample_type = "WBC")

names(df_tumor) <- names(df_WBC) <- c("sample_ID", "sample_name", "mito_CN", "sample_type")
df_combined <- rbind(df_tumor, df_WBC)

p <- ggplot(data = df_combined, aes(x = sample_ID, y = mito_CN, color = sample_type)) +
      geom_point(size = 2) + 
      geom_segment(aes(x=sample_ID, xend=sample_ID, y=0, yend=mito_CN)) + 
      cool_theme + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      geom_hline(yintercept = 0)








      p <- ggplot(data = df_combined, aes(x = sample_ID, y = mito_CN, color = sample_type)) +
      geom_point() + 
      geom_linerange(aes(color = sample_type, ymin = 0, ymax = mito_CN), position = position_dodge(width = 0.5), size = 1) 