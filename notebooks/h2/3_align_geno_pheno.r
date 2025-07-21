#!/usr/bin/Rscript
# SCRIPT 3: Performs genotype and phenotype matrices processing for RARity pipeline
# This script handles pipelines where all chromosomes are processed as ONE RARITY BLOCK;
# or when each chromosome is processed as separate RARITY BLOCK(s), which is the default.
# Inputs:
#    - unaligned/unprocessed clumped genotype and phenotype matrix
#    - covariate ("COVARS_UKB_MERGED_PHENO_BRITISH_NONBRIT_CAUCASIAN_2023_05_22_FINAL_INCLUDE_PCs.RData")
# Steps:
#        1) Align GENOTYPE and PHENOTYPE matrices
#        2) GENOTYPE PROCESSING (remove NA, MAC <2, standardize)
#        3) CONT PHENO PROCESSING (remove NA, quantile norm, adjust covar, standardize)
# Output:
#     - Filtered, processed, aligned GENO in /3_NORM_MAC2_GENO_0.1LD50kWIN_RRV_aligned
#     - Filtered, processed, aligned PHENO in /3_PHENO_aligned

# anno="yhat"
# trait="height"
# DIR_CLUMP="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RVPRS/0c_LD_CLUMP_LOF"
# DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/rovher_chr17_win4096_1000RV"
# Rscript "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pang/Keona_scripts/RARity/exome_blk/3_align_geno_pheno_cont.r" $anno $trait $DIR_WORK $DIR_CLUMP

############################################################################################
args <- commandArgs(trailingOnly = TRUE)
anno <- as.character(args[1]) # score_col
trait <- as.character(args[2])
DIR_WORK <- args[3] # DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/chr${chr}_win${WIN}_${variant_subset}RV"
DIR_CLUMP <- args[4]

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("MBESS"))
suppressMessages(library("corpcor")) #pseudoinverse

total_start_time <- Sys.time()
cat("\n")
cat("=========== SCRIPT 3: Align geno and pheno ===========\n\n")

################################################################################
# STEP 0)    PREPARE PHENOTYPE FILE from UKB data (skip if done)
################################################################################
# /mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RARity_pipeline/SCRIPTS/SCRIPT_0_prepare_pheno.r

################################################################################
# INPUT DIRECTORIES
################################################################################

# Get (un-unprocessed, unaligned) GENOTYPES that were from clumping pipeline  
GENO_DIR <- paste0(DIR_CLUMP,"/8_GENO_0.1LD50kWIN_RDATA")
cat("1. Genotypes :", GENO_DIR, "\n\n")
if (!dir.exists(GENO_DIR) || length(list.files(GENO_DIR)) == 0) { # Check if GENO_DIR doesn't exist or is an empty directory
  stop("Error: GENO_DIR does not exist or is empty")
}

# Get (unaligned) PHENO dir
pheno_dir<-"/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RARity_pipeline/PHENOS_biomarker_2024_02_27/2_pheno_processed" # brit and non-brit imputed values
cat("2. Phenotype :", pheno_dir, "\n\n")

################################################################################
# OUTPUT DIRECTORIES
################################################################################

root <- paste0(DIR_WORK, "/4_exome_h2"); if (!dir.exists(root)) { dir.create(root)}

# Create Output GENO and PHENO directories
GENO_OUTDIR <- paste0(root,"/3_NORM_MAC2_GENO_0.1LD50kWIN_aligned_",trait); if (!dir.exists(GENO_OUTDIR)) { dir.create(GENO_OUTDIR)}
PHENO_OUTDIR <- paste0(root,"/3_PHENO_aligned_",trait); if (!dir.exists(PHENO_OUTDIR)) { dir.create(PHENO_OUTDIR)}

# Output files
norm_df_path<-paste0(PHENO_OUTDIR,"/PHENOS_",trait,"_BRIT_NONBRIT_IMPUTED_aligned_QNORM_RESID_COVAR_normz.RData") 
norm_df_path_no_eid<-paste0(PHENO_OUTDIR,"/PHENOS_",trait,"_BRIT_NONBRIT_IMPUTED_aligned_QNORM_RESID_COVAR_normz_no_eid.RData") # 173599 x 1

##########################################################################################################################################################################################
# Functions
NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)) # needs dataframe, not datatabl
normFunc <- function(x){(x-(mean(x, na.rm = T)))/sd(x, na.rm = T)} #normalsation = x-mean/sd(x) or zVar <- (myVar - mean(myVar)) / sd(myVar)
quantNorm = function(x){qnorm(rank(x,ties.method = "average")/(length(x)+1))}

# ---------- Start of pheno processing  ----------

if (!file.exists(norm_df_path)) {  

  # -----------------------------------------------------------------------------------------------------------
  # STEP 1) Align GENOTYPE and PHENOTYPE matrices
  # -----------------------------------------------------------------------------------------------------------
  
  # Get the reference genotype file (first file in GENO_DIR)
  geno_files <- list.files(GENO_DIR, full.names = TRUE)
  if (length(geno_files) > 0) {
    ref_geno <- geno_files[1]
  } else {
    stop("No files found in GENO_DIR:", GENO_DIR,"\n")
  }
  base::load(ref_geno)
  cat("GENE_DF:", dim(GENE_DF)[1], "x", dim(GENE_DF)[2], "\n")
  colnames(GENE_DF)[1] <- "eid"

  # -----------------------------------------------------------------------------------------------------------
  #      Load phenotype matrix with 'eid' and 'trait' for ALL partipants in BRIT and NON BRIT
  # -----------------------------------------------------------------------------------------------------------
  PHENO <- fread(paste0(pheno_dir,"/", trait,"_BRIT_NONBRIT_IMPUTED.txt"), sep=" ", header=TRUE, data.table=TRUE)
  cat(trait, ":", dim(PHENO)[1], "x", dim(PHENO)[2], "\n")
  print(head(PHENO,2))
  cat("\n")
  
  # -----------------------------------------------
  # Remove 'eid' from pheno that is not in GENE_DF matrix 
  # -----------------------------------------------
  index <- match(GENE_DF$eid, PHENO$eid)
  PHENO <- PHENO[index, ]
  PHENO <- PHENO[order(match(PHENO$eid, GENE_DF$eid)), ]
  PHENO <- PHENO[1:nrow(GENE_DF), ]
  
  # ----------------------------------------------------------------------
  # Remove participants with 'NA' pheno values, and then re-align GENO and PHENO matrix
  # ----------------------------------------------------------------------
  na_count <- sum(is.na(PHENO[,2]))
  cat(trait, "PHENO before NA-removal:", nrow(PHENO), " (Missing values:", na_count, ")\n") #  505 for height, 693 for BMI with NA; 178 for LDL-C

  # PREFERRED: if removing rows with NA, to only contain the eid and adjusted values 
  PHENO <- na.omit(PHENO)
  cat(trait, "PHENO AFTER NA-removal:", nrow(PHENO), "\n\n")

  # Align 'eid' of GENE_DF with PHENO -- takes a long time 
  index <- match(PHENO$eid, GENE_DF$eid)
  GENE_DF <- GENE_DF[index, ] 
  GENE_DF <- GENE_DF[1:nrow(PHENO), ]
  # cat("Any NA in pheno", anyNA(PHENO), "\n")

  # check if the eid columns are identical in G and pheno matrix
  identical_ordered <- identical(GENE_DF$eid, PHENO$eid)
  identical_nrow <- nrow(GENE_DF) == nrow(PHENO)
  
  if(identical_ordered && identical_nrow) {
    cat("--------'eid' in GENO and PHENO are aligned and number of rows match!----------\n")
  } else {
    cat("The GENO and PHENO are not aligned or # of rows don't match.\n")
    if(!identical_ordered) {
      cat("The 'eid' values are not in the same order or is mismatched.\n")
    }
    if(!identical_nrow) {
      cat("Rows in GENE_DF:", nrow(GENE_DF), "   PHENO:", nrow(PHENO), "\n")
    }
  }
  cat("Aligned GENE_DF :", dim(GENE_DF)[1], "x", dim(GENE_DF)[2], " PHENO:", dim(PHENO)[1], "x", dim(PHENO)[2], "\n\n")
  # pheno_file<-paste0(PHENO_OUTDIR, "/PHENOS_", trait, "_top", top, "_BRIT_NONBRIT_IMPUTED_aligned.RData") #  173599 x 2
  # base::save(PHENO, file = pheno_file)

  #####################################################################
  # STEP 2)  CONT PHENOTYPE PROCESSING: 
  #####################################################################

  # STEP 2a: Quantile normalize the pheno column
  # -------------------------------------------------------------------
  cat("PHENO step 1: Quantile normalization....\n")
  PHENO<-as.matrix(PHENO)  # 173183 x 2; from prev step

  # Subset the dataframe to only 1st and 2nd cols
  PHENO <- PHENO[, c(1, 2)]
  PHENO <- PHENO[complete.cases(PHENO), ]  # Remove rows with missing values

  # Quantile normalize the 2nd column
  PHENO[, 2] <- quantNorm(PHENO[, 2])

  # Rename the second column of df1 to 'trait', then add the suffix Q_NORM
  colnames(PHENO)[2] <- trait
  colnames(PHENO)[2] <- paste(colnames(PHENO)[2], "Q_NORM", sep = "_")
  # base::save(PHENO, file = paste0(PHENO_OUTDIR, "/PHENOS_", trait, "_top", top, "_BRIT_NONBRIT_IMPUTED_aligned_QNORM.RData"))

  #      (2b) Adjust for age, sex and 20 PCs 
  # -------------------------------------------------------------------
  cat("PHENO step 2: Adj for age, sex, 20 PCs...\n")
  covar_traits<-c("eid","AGE","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
  base::load("/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/RARity_pipeline/COVARS_UKB_MERGED_PHENO_BRITISH_NONBRIT_CAUCASIAN_2023_05_22_FINAL_INCLUDE_PCs.RData")
  df<-merge(AgeSexPCA,PHENO[,1:2],by="eid") 
  lm_f <- function(x) {x<- residuals(lm(data=df,formula= x ~ AGE+as.factor(SEX)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20))}

  # Extract the trait name to use in the lm_f function
  trait_column <- names(df)[ncol(df)]
  resid <- lm_f(df[[trait_column]])
  resid_final <- data.frame(eid = df$eid)
  resid_final[[trait_column]] <- resid
  # base::save(resid_final, file = paste0(PHENO_OUTDIR, "/PHENOS_", trait, "_top", top, "_BRIT_NONBRIT_IMPUTED_aligned_QNORM_RESID_COVAR.RData"))
  rm(df, resid, AgeSexPCA, covar_traits)
  invisible(gc())

  #      (2c) Standardize mean and variance 
  # -------------------------------------------------------------------
  cat("PHENO step 3: Standardize mean and variance ...\n\n")
  resid_final <- as.matrix(resid_final)

  # Normalization function - (to have mean 0 and variance of 1):
  norm_trait_values <- normFunc(resid_final[, 2])
  second_col_name <- colnames(resid_final)[2] # Get the original second column name

  # Create a new dataframe with normalized trait values
  norm_df <- data.frame(eid = resid_final[, 1])
  norm_df[[second_col_name]] <- norm_trait_values
  cat("Final 'norm_df':", dim(norm_df)[1], "x", dim(norm_df)[2], "\n")
  print(head(norm_df, 2))
  cat("\n")

  # -------------------------------------------------------------------
  # Save a norm_df WITH "eid"  
  base::save(norm_df, file = norm_df_path)
  
  # Save a norm_df WITHOUT "eid" column
  norm_df <-as.matrix(norm_df[,-1]) # 173599 x 1 (remove eid column)
  base::save(norm_df, file=norm_df_path_no_eid) # 173599 x 1

  # clean up
  rm(GENE_DF)
  invisible(gc()) 
} else { 
    cat("Pheno files complete, skip to GENO processing...\n\n")
}

# End of all phenotype processing... skips to here
###################################################################################################
# STEP 3) GENOTYPE PROCESSING
#       (a)  NA-impute 
#       (b)  MINMAC: removes RV columns with a MAC < 2.
#       (c)  NORMALIZE to mean = 0 and sd = 1
#################################################################################
setwd(GENO_OUTDIR) # /3_NORM_MAC2_GENOMATRIX_0.1LD50kWIN_RRV_aligned

# Load aligned and mean-imputed phenotype matrix
base::load(norm_df_path)
setDT(norm_df)

ref_geno <- paste0(GENO_DIR, "/", anno, "_01range_CHR_all.clumped_00.raw.RData")

if (file.exists(ref_geno)) {
  cat("======== SCRIPT 3:", anno, trait, "all chromosomes in one block  ========\n")

  GENO_path <-paste0(GENO_DIR,"/",anno,"_01range_CHR_all.clumped_00.raw.RData")
  outfile <- sprintf(paste0(GENO_OUTDIR,"/NORM_MAC2_BRIT_NONBRIT_IMPUTE_CHR_all_00.RData"))

  if (!file.exists(GENO_path) || file.exists(outfile)) {
    next
  }
  start_time <- Sys.time()
  base::load(GENO_path)  # GENE_DF is loaded
  setDT(GENE_DF)
  colnames(GENE_DF)[1] <- "eid"

  # Merge GENO and norm_df by "eid" using data.table's efficient merge
  GENO_aligned <- merge(GENE_DF, norm_df[, .(eid)], by = "eid", allow.cartesian = TRUE)
  GENE_DF <- GENO_aligned[, names(GENE_DF), with = FALSE] # height GENO aligned:  173183 x 1356 
  if (names(GENE_DF)[ncol(GENE_DF)] == paste0(trait,"_Q_NORM")) {
        GENE_DF <- GENE_DF[, -ncol(GENE_DF), with = FALSE]
  }
  cat(paste("GENO: ", paste(dim(GENE_DF), collapse = " x "), "  PHENO: ", paste(dim(norm_df), collapse = " x "), "\n")) # 173183 x 2 

  # STEP 1: NA-impute
  GENE_DF <- as.data.frame(GENE_DF)
  for (i in 2:ncol(GENE_DF)) {# Apply NA2mean(x) to every column of geno_df starting from column 2
     GENE_DF[, i] <- NA2mean(GENE_DF[, i])
  }

  # STEP 2: Remove columns with MAC < 2
  GENE_DF <- as.matrix(GENE_DF[, colSums(GENE_DF, na.rm = TRUE) >= 2 & !is.nan(colSums(GENE_DF))])
  GENE_DF <- as.data.frame(GENE_DF) 
  cat("MAC=2 filtered:", dim(GENE_DF)[1], "x", dim(GENE_DF)[2], "\n") #  173183 x 1356 

  # STEP 3: Normalize
  for (i in 2:ncol(GENE_DF)) {
      GENE_DF[, i] <- normFunc(GENE_DF[, i])
  }

  # Save GENE_DF to .RData, including "eid' column (ONLY REMOVED IN SCRIPT 4)
  base::save(GENE_DF, file=outfile)
  cat("Duration:", round(as.numeric(Sys.time() - start_time, units = "mins"), 2), "mins\n\n")
  rm(GENE_DF)
  invisible(gc())

} else {

for (chr in 1:22) {
  for (blk in 0:300) {
    
    GENO_path <-paste0(GENO_DIR,"/",anno,"_01range_CHR_",chr,".clumped_", sprintf("%02d", blk), ".raw.RData")
    outfile <- sprintf(paste0(GENO_OUTDIR,"/NORM_MAC2_BRIT_NONBRIT_IMPUTE_CHR_",chr,"_",sprintf("%02d", blk),".RData"))

    if (!file.exists(GENO_path) || file.exists(outfile)) {
      next
    }

    if (file.exists(GENO_path)) {
      cat("======== SCRIPT 3:", anno, trait, "Chr", chr, "blk", sprintf("%02d", blk) ,"========\n")
      start_time <- Sys.time()
      base::load(GENO_path)  # GENE_DF is loaded
      setDT(GENE_DF)
      colnames(GENE_DF)[1] <- "eid"

      # Merge GENO and norm_df by "eid" using data.table's efficient merge
      GENO_aligned <- merge(GENE_DF, norm_df[, .(eid)], by = "eid", allow.cartesian = TRUE)
      GENE_DF <- GENO_aligned[, names(GENE_DF), with = FALSE] # height GENO aligned:  173183 x 1356 
      if (names(GENE_DF)[ncol(GENE_DF)] == paste0(trait,"_Q_NORM")) {
        GENE_DF <- GENE_DF[, -ncol(GENE_DF), with = FALSE]
      }
      cat(paste("GENO: ", paste(dim(GENE_DF), collapse = " x "), "  PHENO: ", paste(dim(norm_df), collapse = " x "), "\n")) # 173183 x 2 

      # STEP 1: NA-impute
      GENE_DF <- as.data.frame(GENE_DF)
      for (i in 2:ncol(GENE_DF)) {# Apply NA2mean(x) to every column of geno_df starting from column 2
        GENE_DF[, i] <- NA2mean(GENE_DF[, i])
      }

      # STEP 2: Remove columns with MAC < 2
      GENE_DF <- as.matrix(GENE_DF[, colSums(GENE_DF, na.rm = TRUE) >= 2 & !is.nan(colSums(GENE_DF))])
      GENE_DF <- as.data.frame(GENE_DF) 
      cat("MAC=2 filtered:", dim(GENE_DF)[1], "x", dim(GENE_DF)[2], "\n") #  173183 x 1356 

      # STEP 3: Normalize
      # cat("STEP 3: Normalized to mean = 0 and sd = 1...\n") # 173183 x 1356
      for (i in 2:ncol(GENE_DF)) {
        GENE_DF[, i] <- normFunc(GENE_DF[, i])
      }

      # Save GENE_DF to .RData, including "eid' column (ONLY REMOVED IN SCRIPT 4)
      base::save(GENE_DF, file=outfile)
      ncols <- ncol(GENE_DF)
      cat("Duration:", round(as.numeric(Sys.time() - start_time, units = "mins"), 2), "mins\n\n")
      rm(GENE_DF)
      invisible(gc())
    } # Done processing one blk 
  } # Done looping through all blks for ONE chromosome
} # Done looping through all chromosomes
}

cat("======== SCRIPT 3: Aligned geno and pheno matrices =========\n")
cat("Number of variants in the last GENE_DF: ", ncols, "\n")
end_time <- Sys.time()
time_taken <- as.numeric(end_time - total_start_time, units = "mins")
cat("Time:", round(time_taken, 2), "mins\n")
cat("Start:", format(total_start_time, "%B %d %H:%M:%S"), "   End:", format(end_time, "%B %d %H:%M:%S"), "\n\n")

