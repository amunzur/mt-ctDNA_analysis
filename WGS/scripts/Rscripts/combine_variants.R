# This script combines variants from mutato caller and VarScan2. 

library(tidyverse)

PATH_mutato <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/annovar/annovar_results/mutato/ALL_MERGED_FILTERED0.01.hg38_multianno.csv"
PATH_varscan <- "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/variant_calling/VarScan2/realigned_mito_filtered/curated_dfs/snv.csv"

mutato <- as.data.frame(read_csv(PATH_mutato)) %>% select(chrom, pos, ref, alt, Func.knownGene, Gene.knownGene, mut_vafs)
varscan <- as.data.frame(read_csv(PATH_varscan)) %>% select(chrom, position, ref, var, Func.knownGene, Gene.knownGene, tumor_var_freq)

mutato$pos <- as.numeric(mutato$pos)
varscan$position <- as.numeric(as.character(varscan$position))

names(mutato) <- names(varscan)

intersected <- inner_join(mutato, varscan, by = c("chrom", "position", "ref", "var"))