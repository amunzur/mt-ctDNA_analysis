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
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
      strip.background = element_blank(), 
      strip.text = element_text(size = 12))
# This script visualizes the variant fractions from one single patient across samples, to show how the VAF changes through disease progression

# 12684 and 12705 
df_180_12684 <- data.frame(sample = c("WBC", "Baseline", "Progression", "Progression2"), VAF = c(0, 4, 31, 35), variant = "12684", patient = "AE-180")
df_180_12705 <- data.frame(sample = c("WBC", "Baseline", "Progression", "Progression2"), VAF = c(1, 4, 30, 34), variant = "12705", patient = "AE-180")

df_18_12684 <- data.frame(sample = c("WBC", "Baseline", "Progression", "Progression2"), VAF = c(4, 6, 28, 40), variant = "12684", patient = "AE-18")
df_18_12705 <- data.frame(sample = c("WBC", "Baseline", "Progression", "Progression2"), VAF = c(4, 6, 29, 39), variant = "12705", patient = "AE-18")

df_151_12684 <- data.frame(sample = c("WBC", "Baseline", "Progression", "Progression2"), VAF = c(1, 6, 6, 27), variant = "12684", patient = "AE-151")
df_151_12705 <- data.frame(sample = c("WBC", "Baseline", "Progression", "Progression2"), VAF = c(2, 6, 7, 27), variant = "12705", patient = "AE-151")

main_df <- rbind( df_180_12684, df_180_12705, df_18_12684, df_18_12705, df_151_12684, df_151_12705)

p <- ggplot(data = main_df, aes(x = factor(sample, levels = c("WBC", "Baseline", "Progression", "Progression2")), y = VAF, fill = patient)) + 
		facet_wrap(~variant) +
		geom_bar(stat = "identity", position = "dodge") + 
		xlab("Sample") + 
		scale_y_continuous(seq(0, 45, 10), name = "Variant Allele Frequency") + 
		geom_hline(yintercept = seq(0, 45, 10), linetype = "dotted") + 
		cool_theme

######################################################################
df_DTB_261_12705 <- data.frame(sample = c("WBC", "Baseline_cfDNA", "Baseline_tumor", "Progression_cfDNA"), VAF = c(2, 30, 0, 13), variant = "12705", patient = "DTB-261")

p <- ggplot(data = df_DTB_261_12705, aes(x = factor(sample, levels = c("WBC", "Baseline_cfDNA", "Baseline_tumor", "Progression_cfDNA")), y = VAF, fill = patient)) + 
		geom_bar(stat = "identity", position = "dodge") + 
		xlab("Sample") + 
		scale_y_continuous(seq(0, 45, 10), name = "Variant Allele Frequency") + 
		geom_hline(yintercept = seq(0, 45, 10), linetype = "dotted") + 
		ggtitle("12705") + 
		cool_theme

######################################################################
df_DTB_216_16184 <- data.frame(sample = c("WBC", "Baseline_cfDNA", "Progression_cfDNA", "Progression_Tumor", "Progression_Tumor2"), VAF = c(6, 69, 66, 3, 3), variant = "16184", patient = "DTB-216")	




path_to_df <- ""