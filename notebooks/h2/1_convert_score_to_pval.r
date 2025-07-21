#!/usr/bin/Rscript
# LD CLUMPING WORKFLOW

# STEP 1) Convert variant scores to p-values (0-1)
# 			Use a known distribution, (e.g.normal distribution), the p-value can be calculated directly using that distribution.

# Example of INPUT score_path:
            #      PLINK_SNP_NAME GENEBASS_AF LOF_DISRUPTIVE      yhat
            # 1: 10:100042482:G:A  2.6588e-06              0 0.9021423

# OUTPUT DIR: $working_dir/1_convert_to_pval
# OUTPUT file: ${score}_01range.txt and ${score}_01range_CHR_1-22.txt 

      # CHR	SNP	BP	A1	A2	pathogenicity_score	rank	P
      # 10	10:100481931:G:C	100481931	G	C	-5.53341715408378	1	3.29765925550726e-07
      # 10	10:100482620:G:A	100482620	G	A	-4.23530024406165	2	3.29765925550726e-07

# score_path="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Training/height_FDR/GENEBASS_RV_train_RV_pred_Oct10/cv_MARS_d3_np15_nv0.1_lowerAF3e-05_upperAF1e-02_PLINK_AF_LOF_refGENE_yhat_200K_missense.txt"
# score_col="yhat"
# score_sort="ascending" # "descending" or "ascending" (lowest to highest)
# working_dir="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RVPRS/0_LD_CLUMP"
# variant_type="LOF" # "missense", "LOF", or "all" (no LOF_DISRUPTIVE column filtering needed)
# Rscript /mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RARity_pipeline/SCRIPTS/LD_CLUMP_SCRIPTS/1_convert_score_to_pval.r $score_path $score_col $score_sort $working_dir $variant_type $variant_subset

###############################################################################################################################################################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
score_path <- args[1] # score file path 
score_col <- args[2] # scores column name 
score_sort <- args[3] # ascending (lowest to highest) or descending (highest to lowest)
working_dir <- args[4] # 
variant_type <- args[5] # "missense", "LOF", or "all" (no LOF_DISRUPTIVE column filtering needed)
variant_subset <- as.numeric(args[6]) # take the top XX subset of variants (rows), after rearranging rows by score_sort

# score_path <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30/RovHer_chr17_win4096_seq3000_all.xlsx"
# score_col <- "evo2_delta_score" # or "yhat"d
# score_sort <- "ascending"
# working_dir <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/1000_RV_win4096"
# variant_type <- "all"
# variant_subset <- 1000

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)
library(readxl)

cat("\n")
cat("---------- Inputs ---------\n\n")
cat("score_path:", score_path, "\n\n")
cat("score_col:", score_col, "\n")
cat("score_sort:", score_sort, "\n\n")
cat("working_dir:", working_dir, "\n\n")
cat("---------------------------\n")

########################################################
# Create directories  
########################################################
dir.create(working_dir, showWarnings = FALSE) # Create root directory 

DIR_OUT <- paste0(working_dir, "/1_convert_to_pval")
dir.create(DIR_OUT) # Create results directory 
cat("DIR_OUT:", DIR_OUT, "\n\n")
setwd(DIR_OUT)

########################################################
# Load input file (scores) 
########################################################

# Read file (given complete path)
if (grepl("\\.xlsx$", score_path)) {
    INPUT_SCORES <- read_excel(score_path)
} else if (grepl("\\.txt$", score_path) || grepl("\\.tsv$", score_path)) {
    INPUT_SCORES <- fread(score_path, header = TRUE, sep = "\t")
} else {
    stop(paste("Unsupported file format for:", score_path))
}
print(head(INPUT_SCORES,2)) 
cat("\n")
cat("INPUT_SCORES: ", dim(INPUT_SCORES)[1], "x", dim(INPUT_SCORES)[2], "\n\n") # INPUT_SCORES:  3279474 x 4 

# Check of column name exists
if (!score_col %in% colnames(INPUT_SCORES)) {
    stop(paste("Column", score_col, "not found in the input file!!\n"))
}

# Remove rows where the column has NA or empty values
INPUT_SCORES <- INPUT_SCORES[!(is.na(INPUT_SCORES[[score_col]]) | INPUT_SCORES[[score_col]] == ""), ]
cat("1. Remove NA: ", dim(INPUT_SCORES)[1], "x", dim(INPUT_SCORES)[2], "\n") # INPUT_SCORES:  3279474 x 4 

####################################################################################
# FILTER [OPTIONAL]
####################################################################################

# STEP 1) Subset variants by MAF
# subsetAF <- 0.0001
# INPUT_SCORES <- subset(INPUT_SCORES, GENEBASS_AF < as.numeric(subsetAF))
# cat("After MAF <", subsetAF, ": ", dim(INPUT_SCORES)[1], "x", dim(INPUT_SCORES)[2], "\n") 

#  STEP 2) Remove LOF variants, keeeping only non-LOF (remove LOF_DISRUPTIVE column has 1)
if ("class" %in% colnames(INPUT_SCORES)) {
    # Rename the column "class" to "LOF_DISRUPTIVE"
    setnames(INPUT_SCORES, old = "class", new = "LOF_DISRUPTIVE")
    
    # Replace "LOF" with 1 and "FUNC/INT" with 0
    INPUT_SCORES[, LOF_DISRUPTIVE := ifelse(LOF_DISRUPTIVE == "LOF", 1, 
                                            ifelse(LOF_DISRUPTIVE == "FUNC/INT", 0, LOF_DISRUPTIVE))]
}
if (variant_type == "missense") {
    INPUT_SCORES <- INPUT_SCORES[INPUT_SCORES$LOF_DISRUPTIVE == 0, ]
    cat("2. After removing LOF: ", dim(INPUT_SCORES)[1], "x", dim(INPUT_SCORES)[2], "\n") #  3171071 (non-LOF) x 4 
} else if (variant_type == "LOF") {
    INPUT_SCORES <- INPUT_SCORES[INPUT_SCORES$LOF_DISRUPTIVE == 1, ]
    cat("2. After removing missense: ", dim(INPUT_SCORES)[1], "x", dim(INPUT_SCORES)[2], "\n") #  3171071 (non-LOF) x 4 
} else if (variant_type == "all") {
    cat("2. No variant class filtering: ", dim(INPUT_SCORES)[1], "x", dim(INPUT_SCORES)[2], "\n") #  3171071 (non-LOF) x 4 
} else {
    stop("Invalid variant_type. Use 'missense', 'LOF', or 'all'.")
}

####################################################################################
# Sort by scores from lowest to highest
####################################################################################

if (score_sort == "descending") {
    INPUT_SCORES <- INPUT_SCORES[order(-INPUT_SCORES[[score_col]]), ]
} else if (score_sort == "ascending") {
    INPUT_SCORES <- INPUT_SCORES[order(INPUT_SCORES[[score_col]]), ]
} else {
    stop("Invalid score_sort value. Use 'ascending' or 'descending'.")
}

# If a variant_subset is specified, take the top XX variants
if (!is.na(variant_subset) && variant_subset > 0 && nrow(INPUT_SCORES) >= variant_subset) {
    INPUT_SCORES <- head(INPUT_SCORES, variant_subset)
    cat("3. After taking top", variant_subset, "variants: ", dim(INPUT_SCORES)[1], "x", dim(INPUT_SCORES)[2], "\n")
} else if (!is.na(variant_subset) && variant_subset > 0) {
    cat("3. Not enough rows in INPUT_SCORES to subset. Using all available rows.\n")
} else {
    cat("3. No variant subset specified. Using all variants.\n\n")
}

####################################################################################
# Assign p-value ranks to the rare variants based on scores (the ranks are normalized).
# Selects columns (chr, PLINK_SNP_NAME, pos, alt, ref, score)
####################################################################################

# Normalization function: to normalize data to a [0, 1] range, which is used later to convert rank values into probabilities (P-values)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
	
# Select 2 columns only 
annodf <- select(INPUT_SCORES, c("PLINK_SNP_NAME", score_col))
cat("Dimension:", dim(annodf)[1], "x", dim(annodf)[2], "\n\n") #  3017195 x 2 

# Extract rows that have no missing values.
annodf<-annodf[complete.cases(annodf), ] 
 
# Rank pathogenicity scores and normalize ranks to get p-values. Computes ranks for the scores with rank(). 
annodf$rank <- rank(annodf[,2], ties.method = "random") # ties.method = "random" handles ties by assigning them a rank randomly.
annodf$P <- range01(annodf$rank)

# Adjustment: Replace zero p-values with the second smallest p-value to avoid having zero p-values which might cause issues in statistical tests.
P2ndlast <- (filter(annodf, annodf$rank == 2))$P
annodf$P <- as.numeric(ifelse(annodf$P == 0, P2ndlast, annodf$P))

# Split the PLINK_SNP_NAME column into ("CHR","SNP","BP", "A1", "A2","P"
split_names <- strsplit(annodf$PLINK_SNP_NAME, ":")

# Create the new columns
annodf$CHR <- sapply(split_names, `[`, 1)
annodf$BP <- sapply(split_names, `[`, 2)
annodf$A1 <- sapply(split_names, `[`, 3)
annodf$A2 <- sapply(split_names, `[`, 4)

# Rename the PLINK_SNP_NAME column in annodf to SNP
colnames(annodf)[colnames(annodf)=="PLINK_SNP_NAME"] <- "SNP"

# Re-order the columns in annodf
annodf <- dplyr::select(annodf, "CHR", "SNP", "BP", "A1", "A2", score_col, "rank", "P")
cat("Final p-value ranks:", dim(annodf)[1], "x", dim(annodf)[2], "\n\n") # RRV:  3017195 x 8; 3169500 x 8 
print(head(annodf,2)) # CHR SNP BP A1 A2 score rank P
cat("\n")

# Split into per-chromosome files
for (ch in 1:22) {
	df<-filter(annodf, CHR==ch)
	df2<-df[order(df$BP),]
	outfile<-as.character(paste0(DIR_OUT,"/",score_col,"_01range_CHR_",ch))
	write.table(df2, outfile,append=F,quote=F,row.names=F,col.names=T,sep="\t")
}

# Write FULL file (Chr 1 -22 included)
write.table(annodf, paste0(DIR_OUT,"/", score_col,"_01range.txt"), append=F,quote=F,row.names=F,col.names=T,sep="\t")
cat("Results:", paste0(DIR_OUT,"/", score_col,"_01range.txt"), "\n\n")
