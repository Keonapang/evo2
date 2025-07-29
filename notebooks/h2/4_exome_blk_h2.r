#!/usr/bin/Rscript
# SCRIPT 4) Generate exome-wide RV heritability estimate for a given bin of RVs
# This script handles pipelines where all chromosomes are processed as ONE RARITY BLOCK;
# or when each chromosome is processed as separate RARITY BLOCK(s), which is the default.

# Output:
#     - Appended to trait_H2_exome_block.txt in '/4_trait_H2_RESULTS'

# Input files:
# 		1. pheno_file: mean-imputed 'resid_final', (before mean = 0, sd = 1)
# 		2. norm_pheno_file: mean-imputed, NORMALIZED 
# 		3. geno_file: one geno 'block' (NA-imputed, MAC > 2, normalized)

# bin_num=1
# anno="yhat"
# trait="height"
# work_dir="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RVPRS"
# DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RVPRS"
# Rscript "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pang/Keona_scripts/RARity/exome_blk/4_exome_blk_h2_cont.r" $anno $trait $threads $DIR_WORK

####################################################################################
rm(list = ls())  
gc()

args <- commandArgs(trailingOnly = TRUE)
anno <- as.character(args[1]) 
trait <- as.character(args[2])
threads <- as.numeric(args[3]) 
DIR_WORK <- args[4]
top <- as.numeric(args[5]) 

# top <- "100"
# DIR_WORK <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/MARS/delta/RovHer_chr17_MARS_d3_np200_nv0.1_lowerAF0e+00_annoyes_embeddelta_blk28"
# trait <- "alanine_aminotransferase"
# threads <- 10 
# anno <- "yhat" # yhat, delta, refvar

if (is.na(threads) || threads <= 0) {
  stop("Error: 'threads' must be a valid positive number.")
}

# install.packages("https://github.com/mrvillage/lmutils.r/archive/refs/heads/master.tar.gz", repos=NULL) # use .zip for Windows
suppressMessages(library("data.table"))
suppressMessages(library("lmutils"))
suppressMessages(library("dplyr"))
suppressMessages(library("MBESS"))
suppressMessages(library("corpcor")) #pseudoinverse

cat("\n=========== SCRIPT 4: Exome-wide h2 RV for", trait, "===========\n")
start_time <- Sys.time()
cat("START:", format(start_time, "%B %d %H:%M:%S"), "\n\n")

################################################################################
# INPUT
################################################################################
root1 <- paste0(DIR_WORK, "/4_exome_h2")
root <- paste0(root1, "/top",top)
cat("1. Output root :", root, "\n\n")

# GENO and PHENO directories
GENO_DIR <- paste0(root, "/3_NORM_MAC2_GENO_0.1LD50kWIN_aligned_",trait)
PHENO_DIR <- paste0(root, "/3_PHENO_aligned_",trait)
if (!file.exists(GENO_DIR) || length(list.files(GENO_DIR))==0) { 
  cat("2. GENO_DIR:", GENO_DIR, "\n\n")
  stop("Geno dir missing\n")
}
if (!file.exists(PHENO_DIR) || length(list.files(PHENO_DIR)) == 0) { 
  cat("3. PHENO_DIR:", PHENO_DIR, "\n\n")
  stop("Pheno dir missing\n")
}
# Input files
norm_df_path <- paste0(PHENO_DIR,"/PHENOS_",trait,"_BRIT_NONBRIT_IMPUTED_aligned_QNORM_RESID_COVAR_normz.RData") 
norm_df_path_no_eid <- paste0(PHENO_DIR,"/PHENOS_",trait,"_BRIT_NONBRIT_IMPUTED_aligned_QNORM_RESID_COVAR_normz_no_eid.RData") # 173599 x 1

################################################################################
# OUTPUT
################################################################################

# directory
H2_RESULTS <-paste0(root, "/4_",trait,"_H2_RESULTS/"); if (!dir.exists(H2_RESULTS)) { dir.create(H2_RESULTS)}
setwd(H2_RESULTS) 

# MODIFY output files if needed!
outfile1 <- paste0(H2_RESULTS, trait, "_H2_exome_raw.txt") 
outfile2 <- paste0(H2_RESULTS, trait, "_H2_exome_raw_processed.txt")
outfile3 <- paste0(H2_RESULTS, "TOTAL_H2_", trait, "_exome.txt")
if (file.exists(outfile3)) {
  stop(paste("\n\nOutput already exists: ", outfile3, "\n\n"))
}

################################################################################
#  1.  Load Phenotype and configure the file
################################################################################

if (!file.exists(norm_df_path_no_eid)) {  # ---------- Start of pheno processing  ----------
base::load(norm_df_path) # norm_df object
norm_df <-as.matrix(norm_df[,-1]) # 173599 x 1 (remove eid column)
base::save(norm_df, file=norm_df_path_no_eid)
} else {
  base::load(norm_df_path_no_eid) # norm_df object
  cat("\nnorm_df :", norm_df_path_no_eid, "\n")
}

################################################################################
# 2. Count the number of genotype data
################################################################################
list_of_geno <- list.files(GENO_DIR, pattern = "\\.rkyv\\.gz$", full.names = TRUE) 
rkyv_count <- length(list_of_geno)
list_of_rdata <- list.files(GENO_DIR, pattern = "\\.RData$", full.names = TRUE)
rdata_count <- length(list_of_rdata)

################################################################################
# 3. Calculate h2 
################################################################################

get_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}
start_time <- Sys.time()

cat("=========== [ SCRIPT 4 :", anno, trait, "h2 RV] ==========","\n\n")                       

if (rdata_count > 0) {
  cat("\nCALCULATING h2 for", rdata_count, " RData files... \n")

  lmutils::set_core_parallelism(threads)
  lmutils::set_num_worker_threads(threads)
  lmutils::set_num_main_threads(threads) #  default 16; number of main threads to use - number of primary operations to run in parallel.
  results <- lmutils::calculate_r2(list_of_rdata, norm_df)
  cat("Results: ", paste(dim(results), collapse = " x "), "\n")
    #  block  n       m          r2           adj_r2
    #  1      50      5000    0.05972631     2.233549e-03

  # Subset the results data frame to include only the specified columns
  results_subset <- results[, c("r2", "adj_r2", "n", "m", "data")]

  # Validate if all entries in the "data" column contain the required format
  valid_data <- grepl("CHR_[0-9]+_[0-9]+", results_subset$data)

  if (all(valid_data)) {
    cat("All file paths are valid. Truncating and rearranging...\n")
    
    # Truncate the filepath in column "data" to keep only chr # and blk #
    results_subset$data <- sapply(strsplit(as.character(results_subset$data), "/"), function(x) tail(x, 1))
    results_subset$data <- sub(".*_(CHR_[0-9]+_[0-9]+).*", "\\1", results_subset$data)

    # Rearrange results_subset by the values in the "data" column so they are in chronological CHR order
    results_subset <- results_subset[order(
      as.numeric(sub("CHR_([0-9]+)_.*", "\\1", results_subset$data)), 
      as.numeric(sub("CHR_[0-9]+_([0-9]+)", "\\1", results_subset$data))
    ), ]
  } else {
    cat("Done processing a single rarity block.\n")
  }
}

write.table(results_subset, file = outfile1, append = F, quote = F, row.names = F, col.names = T, sep = "\t")
cat(paste("Output1: ", outfile1, "\n\n"))

################################################################################
# 4. Calculate CIs, block_LCL_R2, block_UCL_R2 adj_R2_perVar, block_VarR2, block_Var_adj_R2
################################################################################

for (i in seq_len(nrow(results_subset))) {
  r2 <- as.numeric(results_subset$r2[i])
  adj_r2 <- as.numeric(results_subset$adj_r2[i])
  n <- as.numeric(results_subset$n[i])
  m <- as.numeric(results_subset$m[i])

  # Check for NA values or unexpected range (example range check shown, adjust as needed)
  if(is.na(r2) || is.na(adj_r2) || is.na(n) || is.na(m) || n <= m || m <= 0) {
    next 
  }
  ci_output <- ci.R2(r2, df.1 = m, df.2 = n - m, conf.level = 0.95)
  results_subset$block_LCL_R2[i] <- ci_output$Lower.Conf.Limit.R2
  results_subset$block_UCL_R2[i] <- ci_output$Upper.Conf.Limit.R2
  results_subset$adj_R2_perVar[i] <- adj_r2 / m
  rho2 <- r2
  block_VarR2 <- ((4 * rho2) * (1 - rho2)^2 * (n - m - 1)^2) / ((n^2 - 1) * (n + 3))
  results_subset$block_VarR2[i] <- block_VarR2
  results_subset$block_Var_adj_R2[i] <- ((n - 1) / (n - m - 1))^2 * block_VarR2
}
results_subset$trait <- trait
results_subset <- results_subset[, c("trait", setdiff(names(results_subset), "trait"))]

# save 
write.table(results_subset, file = outfile2, row.names = F, quote = F, sep = "\t")
cat(paste("\nOutput blocks: ", outfile2, "\n\n"))

#################################################################################
#   5. Total h2 calculation for single-trait
#################################################################################

results_subset_tot <- fread(outfile2)
setDT(results_subset_tot)

# # First, ensure adj_r2 is numeric if it's stored as a character
results_subset_tot[, adj_r2 := as.numeric(adj_r2)]
results_subset_tot[, duplicate := r2 == adj_r2]
results_subset_tot <- results_subset_tot[!duplicate | !duplicated(results_subset_tot[, .(r2, adj_r2)]), ]
results_subset_tot[, duplicate := NULL] # Remove the temporary column

# Get results/summary for output 
N_RVs <- base::sum(results_subset_tot$m, na.rm = TRUE) # Total number of RVs (m), ignoring NA values
N<-base::mean(results_subset_tot$n, na.rm = TRUE) # Mean number of individuals (n)
ADJ_R2 <- base::sum(as.numeric(results_subset_tot$adj_r2), na.rm = TRUE) # Total adj_R2, skipping NAs
STD2 = base::sqrt(base::sum(results_subset_tot$block_Var_adj_R2)) ## check if this is supposed to be block var r2?
LCL_adj = ADJ_R2 - 1.96 * STD2 
UCL_adj = ADJ_R2 + 1.96 * STD2 

# Define headers
output_headers <- c("TRAIT", "N", "N_RVs", "ADJ_R2", "STD2", "LCL_adj", "UCL_adj")
data_to_write <- cbind(trait, N, N_RVs, ADJ_R2, STD2, LCL_adj, UCL_adj)

if (!file.exists(outfile3)) {
    cat(paste("Creating: ", outfile3, "\n"))
    suppressWarnings(write.table(data_to_write, file = outfile3, append = T, quote = F,
                                row.names = F, col.names = output_headers, sep = "\t"))
	} else {
    cat(paste("Appending: ", outfile3, "\n\n"))
    suppressWarnings(write.table(data_to_write, file = outfile3, append = T, quote = F,
                                row.names = F, col.names = F, sep = "\t"))
}
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat(trait, "h2: ", ADJ_R2, "\n")
cat("RVs (m):", N_RVs, "     Participants (n):", N, "\n\n")
cat("Duration:", round(duration), "mins      End:", format(end_time, "%B %d %H:%M:%S"), "\n")
cat("========= SCRIPT 5:", anno, trait, " =======\n\n")
cat("DIR_WORK:", DIR_WORK, "\n")

#################################################################################
# 6. CLEAN UP: all .RData files from dir
#################################################################################

final_result <- fread(outfile3, header = TRUE, sep = "\t")
if (!is.null(final_result) && ncol(final_result) > 0 && nrow(final_result) > 0 && all(!is.na(final_result))) {
  unlink(GENO_DIR, recursive = TRUE)
  unlink(PHENO_DIR, recursive = TRUE)
  file.remove(outfile1)
} else {
  cat("Number of .RData:", rdata_count, "\n\n")
  cat("GENO_DIR:", GENO_DIR, "\n\n")
}
