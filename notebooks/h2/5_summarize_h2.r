#!/usr/bin/Rscript
# SCRIPT 6) Generate exome-wide RV heritability estimate for a given bin of RVs
# This script handles pipelines where all chromosomes are processed as ONE RARITY BLOCK;
# or when each chromosome is processed as separate RARITY BLOCK(s), which is the default.

# Output:
#     - Appended to trait_H2_exome_block.txt in '/4_trait_H2_RESULTS'

# Input files:
# 		1. pheno_file: mean-imputed 'resid_final', (before mean = 0, sd = 1)
# 		2. norm_pheno_file: mean-imputed, NORMALIZED 
# 		3. geno_file: one geno 'block' (NA-imputed, MAC > 2, normalized)

# WIN="4096" # 4096
# chr="17"
# subsets=("top" "bottom") # top or bottom?
# variant_subsets=("1000" "5000") # top and bottom XX subset of variants 
# DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2"

# for variant_subset in "${variant_subsets[@]}"; do
# for subset in "${subsets[@]}"; do
#   Rscript "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/h2/5_summarize_h2.r" $DIR_WORK $WIN $chr $subset $variant_subset
# done
# done

####################################################################################
rm(list = ls())  
gc()

args <- commandArgs(trailingOnly = TRUE)
DIR_WORK <- args[1]
WIN <- args[2] # ="4096" # 4096
chr <- args[3] # ="17"
subset <- args[4] # "top" # top or bottom XX variants ?
variant_subset <- args[5] # "1000" # subset rows from the top; "none" means no subsetting needed

# WIN="4096" # 4096
# chr="17"
# subset="top" # top or bottom?
# variant_subset="1000" # subset rows from the top; "none" means no subsetting needed
# DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2"

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
library(openxlsx)
library(ggplot2)

cat("\n")

# loops
scores <- c("rovher", "evo2")
traits <- c("height", "LDL_direct", "Alkaline_phosphatase", "BMI","ApoB" ,"ApoA", "Glycated_haemoglobin" ,"Triglycerides" ,"IGF_1" ,"alanine_aminotransferase" ,"C_reactive_protein", "Phosphate")
# subsets <- c("top", "bottom")

# Directory of TOTAL heritability of Chromosome XX
DIR_TOT <- paste0("/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/chr", chr, "_win",WIN,"_all")
cat(paste0("Tot h2 directory: ", DIR_TOT, "\n\n"))

# ------------------------------------- Begin -----------------------------------

results_df <- data.table(
  SCORE = character(),
  CHR = character(),
  WIN = character(),
  SUBSET = character(),   # top or bottom
  VAR_ROWS = numeric(),   # 1000, 5000
  TRAIT = character(), 
  N = numeric(), 
  N_RVs = numeric(), 
  ADJ_R2 = numeric(), 
  TOT_R2 = numeric(),
  PROP_H2 = numeric(),
  STD2 = numeric(), 
  LCL_adj = numeric(), 
  UCL_adj = numeric()
)

for (score in scores) {
  counter <- 1

  for (trait in traits) {

        # Total h2
        tot_file <- paste0(DIR_TOT, "/4_exome_h2/4_", trait, "_H2_RESULTS/TOTAL_H2_", trait, "_exome.txt")
        if (file.exists(tot_file)) {
          data_tot <- fread(tot_file)
          last_row <- data_tot[nrow(data_tot), ]
          TOT_R2 <- last_row$ADJ_R2
        } else {
          TOT_R2 <- 0.1
        }

        # Data h2 for a subset of variants
        H2_RESULTS <- paste0(DIR_WORK, "/", score, "_chr", chr, "_win", WIN, "_", variant_subset, "RV_", subset, "/4_exome_h2/4_", trait, "_H2_RESULTS/")
        infile <- paste0(H2_RESULTS, "TOTAL_H2_", trait, "_exome.txt")

        if (file.exists(infile)) {
          data_h2 <- fread(infile, header = TRUE, data.table = FALSE)
          last_row <- data_h2[nrow(data_h2), ]

          results_df <- rbind(results_df, data.table(
            SCORE = score, 
            CHR = chr, WIN = WIN, SUBSET = subset, VAR_ROWS = variant_subset,
            TRAIT = trait, N = last_row$N, N_RVs = last_row$N_RVs, 
            ADJ_R2 = last_row$ADJ_R2, TOT_R2 = TOT_R2,
            PROP_H2 = (as.numeric(last_row$ADJ_R2) / TOT_R2) * 100,
            STD2 = last_row$STD2, LCL_adj = last_row$LCL_adj, UCL_adj = last_row$UCL_adj
          ))
          counter <- counter + 1
          cat(paste0(counter, ": ", trait, "    ", round(((as.numeric(last_row$ADJ_R2) / TOT_R2) * 100),2), "% \n"))
        } else {
          cat(paste0("MISSING ", infile, "\n\n"))
        }
  }# --- done looping traits
  
  outfile <- paste0(DIR_WORK, "/", score,"_chr",chr,"_win",WIN,"_", variant_subset, "RV_",subset,"/H2_summary.xlsx")
  write.xlsx(results_df, file = outfile, rowNames = FALSE, colNames = TRUE, append = TRUE)
  cat(paste0("\nResults: ", outfile, "\n\n"))

        results_df <- data.table(
        SCORE = character(),
        CHR = character(),
        WIN = character(),
        SUBSET = character(),   # top or bottom
        VAR_ROWS = numeric(),   # 1000, 5000
        TRAIT = character(), 
        N = numeric(), 
        N_RVs = numeric(), 
        ADJ_R2 = numeric(), 
        TOT_R2 = numeric(),
        PROP_H2 = numeric(),
        STD2 = numeric(), 
        LCL_adj = numeric(), 
        UCL_adj = numeric()
      )
} # --- done looping both scores


########################################################
# Plot bar plot
########################################################

# WIN <- "4096" # 4096
# chr <- "17"
# subset <- "bottom" # top or bottom?
# variant_subset <- "5000" # subset rows from the top; "none" means no subsetting needed
# DIR_WORK <-"/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2"
# score <- "rovher" # or "evo2"


# load results
outfile1 <- paste0(DIR_WORK, "/rovher_chr",chr,"_win",WIN,"_", variant_subset, "RV_",subset,"/H2_summary.xlsx")
df1 <- read.xlsx(outfile1)
outfile2 <- paste0(DIR_WORK, "/evo2_chr",chr,"_win",WIN,"_", variant_subset, "RV_",subset,"/H2_summary.xlsx")
df2 <- read.xlsx(outfile2)

# combine df
df <- rbind(df1, df2)
df$SCORE <- gsub("evo2", "Evo2 7B", df$SCORE)
df$SCORE <- gsub("rovher", "RovHer", df$SCORE)

# Ensure `TRAIT` is a factor for proper ordering in the plot
df$TRAIT <- gsub("LDL_direct", "LDL-C", df$TRAIT)
df$TRAIT <- gsub("Alkaline_phosphatase", "alkaline\nphosphatase", df$TRAIT)
df$TRAIT <- gsub("Glycated_haemoglobin", "HbA1c", df$TRAIT)
df$TRAIT <- gsub("IGF_1", "IGF-1", df$TRAIT)
df$TRAIT <- gsub("alanine_aminotransferase", "alanine\naminotransferase", df$TRAIT)
df$TRAIT <- gsub("C_reactive_protein", "CRP", df$TRAIT)
df$TRAIT <- gsub("Phosphate", "phosphate", df$TRAIT)

df$TRAIT <- factor(df$TRAIT, levels = unique(df$TRAIT))

font <- "Open Sans"  

# get the maximum value for y-axis
ymax <- max(df$PROP_H2, na.rm = TRUE)
ymin <- min(df$PROP_H2, na.rm = TRUE)

if (variant_subset == "1000") {
  ymax <- 25
  ymin <- -10
} else {
  ymax <- 60
  ymin <- -10
}

# take the value under column N_RVs from the first row of df
n_rvs <- df$N_RVs[1]

p <- ggplot(df, aes(x = TRAIT, y = PROP_H2, fill = SCORE)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  labs(
    title = paste0("Chromosome 17 exome-wide RV heritability\n\n", subset," ",variant_subset, " of RVs (clumped to ", n_rvs, " RVs)"),
    x = "Trait",
    y = "Proportion of h2 explained (%)"
  ) +
  scale_y_continuous(limits = c(ymin, ymax)) + # coord_cartesian(ylim = c(0.1, 4))
  theme_minimal() +
  theme(
  panel.background = element_rect(fill = "white", color = NA),  # White background
  plot.background = element_rect(fill = "white", color = NA),   # White plot area
  plot.margin = margin(t = 15, r = 15, b = 20, l = 20),         # Adjust plot margins
  panel.grid.major = element_blank(),  # Removes major gridlines
  panel.grid.minor = element_blank(),   # Removes minor gridlines
  axis.title.x = element_text(margin = margin(t = 15, b = 0.8), family = font),  # Padding for x-axis
  axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0.8), family = font),  # Padding for y-axis
  axis.text = element_text(size = 14, family = font),  # Increase axis text size
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11, family = font),  # Rotate x-axis text
  axis.title = element_text(size = 13, face = "bold", family = font),  # Increase axis title size
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5, family = font),  # Centered, bold title
  legend.text = element_text(size = 10, family = font),
  legend.title = element_blank()
  ) +
  scale_fill_manual(values = c("Evo2 7B" = "steelblue", "RovHer" = "orange")) # Custom colors for bars
outplot <- paste0(DIR_WORK, "/chr",chr,"_win",WIN,"_", variant_subset, "RV_", subset, ".png")
ggsave(outplot, p, width = 10, height = 6, dpi = 300)
cat(paste0("Bar plot: ", outplot, "\n\n"))


