library(tidyverse)

normal_threshold <- 0.10
normal_read_support <- 4 # variants with less than 10 rs are EXCLUDED
tumor_read_support <- 4

dir_unfiltered_results <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarScan2/realigned_mito_filtered/snv"
dir_filtered_results <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarScan2/realigned_mito_filtered/filtered_output/snv"

#INPUTS
process_samples <- function(dir_unfiltered_results, dir_filtered_results, vaf_threshold, normal_read_support, tumor_read_support) {

	setwd(dir_unfiltered_results)
    paths_list <- list.files(dir_unfiltered_results, pattern = "filtered")
    files <- lapply(as.list(paths_list), read.delim)
    # files <- lapply(as.list(files), as.data.frame)

    # names(files) <- paths_list

    i <- 1
    while (i <= length(files)){
        file <- as.data.frame(files[[i]])
       
        if (dim(file)[1] > 0) {
            
            file$samplename <- paths_list[i]
            files[[i]] <- file
        }
       
        i <- i + 1
    }


    files <- files[sapply(files, nrow) > 0] # only keep dfs with at least one variant
    combined <- do.call(rbind, files)
    write_delim(combined, file.path(dir_unfiltered_results, "COMBINED.tsv"), delim = "\t")

    combined$tumor_var_freq <- as.numeric(gsub('%', '', combined$tumor_var_freq))/100
    combined$normal_var_freq <- as.numeric(gsub('%', '', combined$normal_var_freq))/100

    combined <- combined[!is.na(combined$normal_var_freq), ]

    combined[, setdiff(5:11, 8)] <- apply(combined[, setdiff(5:11, 8)], 2, as.character)
    combined[, setdiff(5:11, 8)] <- apply(combined[, setdiff(5:11, 8)], 2, as.numeric)

    filtered <- combined %>%
        mutate(VAF_factor = round(tumor_var_freq / normal_var_freq), 
                normal_reads_total = normal_reads2 + normal_reads1) %>%
        filter(VAF_factor >= 3, 
                normal_reads_total >= normal_read_support, 
                tumor_reads2 >= tumor_read_support, 
                normal_var_freq <= normal_threshold)

    write_delim(filtered, file.path(dir_filtered_results, "FILTERED.tsv"), delim = "\t")










}


