#!/usr/bin/R
# Date: July 29, 2025
####################################################################################
# SCRIPT: GET an excel table summary of the h2_curve for 8 predictors for RRV and RV
# A single excel spreadsheet output with 4 tabs:
    # Tab 1: RV
    # Tab 2: RRV
    # Tab 3: Summary
    # Tab 4: Tallies

# top_list=(1,5,10)
# NN_DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RARity_h2curve/NN"
# MARS_DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RARity_h2curve/MARS"
# rovher_DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RARity_h2curve/rovher"
# Rscript "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RARity_h2curve/get_h2_curve_excel_summary.r" $top_list

####################################################################################

args <- commandArgs(trailingOnly = TRUE)
top_list <- args[1] # top_list=(1,5,10)
NN_DIR <- args[2] # 
MARS_DIR <- args[3] # 
rovher_DIR <- args[4] # 
# top_list <- (1)

# root <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2"
# NN_DIR <- paste0(root, "/NN/refvar/RovHer_chr17_blocks.28.mlp.l3_height_FDR_annoyes")
# MARS_DIR<- paste0(root, "/MARS/delta/RovHer_chr17_MARS_d3_np200_nv0.1_lowerAF0e+00_annoyes_embeddelta_blk28")
# rovher_DIR<- paste0(root, "/RovHer/embedno/RovHer")

suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

########################################################
# Variables
########################################################
teal <- "#0B7C9A"
red <- "#DE433C"
green <- "#76B900"

models <- list(
  "RovHer (anno)" = list(dir = paste0(MARS_DIR, "/4_exome_h2"), color = "steelblue"),
  "RovHer (anno + delta)" = list(dir = paste0(rovher_DIR, "/4_exome_h2"), color = "#76B900"),
  "Neural Network (anno + ref + var)" = list(dir = paste0(NN_DIR, "/4_exome_h2"), color = "orange")
)

# Convert prop_blks to a numeric vector if it is a character
if (is.character(top_list)) {
  top_list <- as.numeric(unlist(strsplit(top_list, ",")))
}
trait_list <- c("height", "BMI", "alanine_aminotransferase","albumin", "Alkaline_phosphatase", "ApoA", 
                "ApoB", "Aspartate_aminotransferase", "Calcium", "C_reactive_protein", 
                "Creatinine", "Cystatin_C", "Gamma_glutamyltransferase", "Glycated_haemoglobin", "IGF_1", 
                "LDL_direct", "Phosphate", "Total_protein", "Triglycerides" ,"Urate", "Urea") 

trait_rename <- c("Height", "Body mass index", "Alanine aminotransferase", "Albumin", "Alkaline phosphatase", "Apolipoprotein A", 
                    "Apolipoprotein B", "Aspartate aminotransferase", "Calcium", "C-reactive_protein", "Creatinine", 
                    "Cystatin-C", "Gamma-glutamyl transferase", "Glycated haemoglobin", "IGF-1", 
                    "LDL cholesterol", "Phosphate", "Total protein", "Triglycerides" ,"Urate", "Urea")

########################################################
# Output directory 
########################################################

DIR_OUT <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/plot_chr17"
if (!dir.exists(DIR_OUT)) {
  dir.create(DIR_OUT, recursive = TRUE)
}

########################################################
# Loop through each model and each trait
########################################################

for (top in top_list) {

df <- data.frame(TRAIT = character(0), top = numeric(0), model = character(0), 
            PROP_H2 = numeric(0), ADJ_R2 = numeric(0), 
            TOT_ADJ_R2 = numeric(0), STD2 = numeric(0), 
            LCL_adj = numeric(0), UCL_adj = numeric(0), 
            N_RVs = numeric(0), N = numeric(0)

)  
  for (model in names(models)) {
  
  # Access directory from the list
  dir <- models[[model]]$dir
  cat("\nProcessing:", model, "\n")
  cat("Directory:", dir, "\n\n")
    
      for (TRAIT in trait_list) {
          cat(TRAIT, " - ", model, "\n")

          # Get total r2 from the top100 file 
          tot_file <- paste0(dir, "/top100/4_", TRAIT, "_H2_RESULTS/TOTAL_H2_", TRAIT, "_exome.txt")
          if (file.exists(tot_file)) {
            tot_data <- fread(tot_file)
            tot_row <- tot_data[nrow(tot_data), ]
            tot_r2 <- tot_row$ADJ_R2
          } else {
            tot_r2 <- 0
          }

            
            in_file <- paste0(dir, "/top", top,"/4_", TRAIT, "_H2_RESULTS/TOTAL_H2_", TRAIT, "_exome.txt")
            out_plot <- paste0(DIR_OUT, "/top",top,".png")
            out_excel <- paste0(DIR_OUT, "/top",top,".xlsx")

            if (file.exists(in_file)) {
              data <- fread(in_file)

              # extract the last line of the file
              row <- data[nrow(data), ]

              prop_h2 <- (row$ADJ_R2/tot_r2)*100
              if (tot_r2 == 0) {
                prop_h2 <- 0
              }
              df <- rbind(df, data.table(TRAIT = TRAIT, model = model, top = top,
                                        PROP_H2 = prop_h2, ADJ_R2 = row$ADJ_R2, 
                                        TOT_ADJ_R2 = row$TOT_ADJ_R2, STD2 = row$STD2, 
                                        LCL_adj = row$LCL_adj, UCL_adj = row$UCL_adj, 
                                        N_RVs = row$N_RVs, N = row$N))
            } else {
                cat("\nMISSING: ", in_file, "\n")
            } 
        } 
    }
}

rename_map <- setNames(trait_rename, trait_list)
df[, TRAIT := rename_map[TRAIT]]
cat("\n")
dim(df)
head(df,2)

# Write the data frame to an Excel file
write.xlsx(df, file = out_excel, sheetName = "top1")
cat("\nSaved: ", out_excel, "\n")

################################################################
# Bar plot
################################################################

for (top in top_list) {
# load results
df <- read.xlsx(out_excel)

# output
outplot <- paste0(DIR_OUT, "/top_", top, ".png")

font <- "Open Sans"  

# Ensure `TRAIT` is a factor for proper ordering in the plot
df$TRAIT <- gsub("LDL_direct", "LDL-C", df$TRAIT)
df$TRAIT <- gsub("Alkaline_phosphatase", "alkaline\nphosphatase", df$TRAIT)
df$TRAIT <- gsub("Glycated_haemoglobin", "HbA1c", df$TRAIT)
df$TRAIT <- gsub("IGF_1", "IGF-1", df$TRAIT)
df$TRAIT <- gsub("alanine_aminotransferase", "alanine\naminotransferase", df$TRAIT)
df$TRAIT <- gsub("C_reactive_protein", "CRP", df$TRAIT)
df$TRAIT <- gsub("Phosphate", "phosphate", df$TRAIT)
df$TRAIT <- factor(df$TRAIT, levels = unique(df$TRAIT))

# get the maximum value for y-axis
ymax <- max(df$PROP_H2, na.rm = TRUE)
ymin <- min(df$PROP_H2, na.rm = TRUE)

# Create legend text
legend_text <- sapply(models, function(x) x$color)
names(legend_text) <- names(models)
# print(legend_text)

# take the value under column N_RVs from the first row of df
n_rvs <- df$N_RVs[1]

p <- ggplot(df, aes(x = TRAIT, y = PROP_H2, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  labs(
    title = paste0("Chromosome 17 exome-wide RV heritability\n\nTop ",top, "% of RVs (m=", n_rvs,")"),
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
  scale_fill_manual(values = legend_text) # Custom colors for bars

ggsave(outplot, p, width = 10, height = 6, dpi = 300)
cat(paste0("Bar plot: ", outplot, "\n\n"))
}
