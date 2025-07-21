#!/usr/bin/Rscript
# To combine Evo2 7B scores for chromosome XX with a given window size and sequence length (i.e. chunk size)
# This script reads in multiple Excel files (SEQ_LENGTH number of chunks) and combines them into a
# single text file by rbind() sequentially.

# Input files:
#     - RovHer_chrX_1_winXXXX_seqYYYY_all.xlsx, ..., RovHer_chrX_101_winXXXX_seqYYYY_all.xlsx

# Output:
#     - Combined Evo2 scores for chromosome XX in a single file: evo2_chrXX_winXXXX.txt 

# DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30"
# WIN="4096" # 4096
# CHR="17"
# SEQ_LENGTH="3000"
# Rscript "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/evo2_pilot_study/rbind_evo2_scores.r" $DIR_WORK $WIN $CHR $SEQ_LENGTH

####################################################################################
rm(list = ls())  
gc()

args <- commandArgs(trailingOnly = TRUE)
DIR_WORK <- args[1]
WIN <- args[2] # ="4096" # 4096
CHR <- args[3] # ="17"
SEQ_LENGTH <- args[4] # # subset rows from the top; "none" means no subsetting needed

# WIN="4096" # 4096
# CHR="17"
# SEQ_LENGTH="3000" # subset rows from the top; "none" means no subsetting needed
# DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30"

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
library(readxl) # For reading Excel files
cat("\n")

# output file
outfile <- paste0(DIR_WORK, "/evo2_chr",CHR,"_win",WIN,".txt")
cat("Output: ", outfile, "\n")

# counter
total_rows <- 0

# ------------------------------------------------------------------------

# Initialize the output file by writing the header of the first file
first_file <- paste0(DIR_WORK, "/RovHer_chr", CHR,"_1_win", WIN, "_seq", SEQ_LENGTH, "_all.xlsx")
# RovHer_chr17_1_win4096_seq3000_all
if (file.exists(first_file)) {
  data <- read_xlsx(first_file)
  total_rows <- total_rows + nrow(data)

  fwrite(data, file = outfile, sep = "\t", append = FALSE) # Create the file and write the header and first file's data
} else {
  stop("Not found: ", first_file, "\n")
}
# ------------------------------------------------------------------------

# Loop through the remaining files (2 to 300) 
for (i in 2:300) {
  file <- paste0(DIR_WORK, "/RovHer_chr", CHR,"_",i, "_win", WIN, "_seq", SEQ_LENGTH, "_all.xlsx")

  if (file.exists(file)) {
    cat("Processing: ", i, "\n")
    data <- read_xlsx(file)
    total_rows <- total_rows + nrow(data)

    fwrite(data, file = outfile, sep = "\t", append = TRUE, col.names = FALSE)
  } else {
      next
  }
}

cat("\nSaved to: ", outfile, "\n\n")
cat("Total rows (excl. headers): ", total_rows, "\n\n")


# intialize a readme file and output the total number of rows
readme_file <- paste0(DIR_WORK, "/README.txt")
cat("This file contains the combined Evo2 7B scores for chromosome ", CHR,
 " with each variant evaluated with a window size of ", WIN, " and each input sequence length is ", SEQ_LENGTH, "nt.\n", 
  "Total number of rows (excluding headers): ", total_rows, "\n", 
  "Each row corresponds to a different variant.\n", 
  file = readme_file)
cat("Readme file created at: ", readme_file, "\n")



