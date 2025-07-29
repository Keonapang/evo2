#!/usr/bin/Rscript
# SCRIPT 2
# OBJECTIVE: Loop through clumped genotype blocks, and extract only the RVs found in input variant list 
# Input: 
#     - PLINK_list consisting of clumped variant ids
#     - All clumped ~250 exome-wise (chr1-22) GENO blks from '/8_GENO_0.1LD50kWIN_RDATA'
# Output:
#     - Smaller/filtered set of ~5k genotype blks for RVs in each proportion in /2_GENO_0.1LD50kWIN_RDATA

# DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/MARS/refvar/RovHer_chr17_MARS_d3_np150_nv0.1_lowerAF0e+00_annono_embedrefvar_blk28"
# PLINK_list="${DIR}/4_CLUMP_RESULT/yhat_01range_CHR_17.clumped"
# prop_blks=(1,5,10)
# cores=10
# Rscript /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/h2/2_geno_blks.r $DIR $input_file $prop_blks $cores

############################################################################################################
args <- commandArgs(trailingOnly = TRUE)
 
DIR <- args[1] # DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/MARS/refvar/RovHer_chr17_MARS_d3_np150_nv0.1_lowerAF0e+00_annono_embedrefvar_blk28"
PLINK_list <- args[2] # input_file="${DIR}/4_CLUMP_RESULT/${anno_name}_01range_CHR_17.clumped"
prop_blks <- as.character(args[3]) # "1,5,10,15,20,25,30,35,40,50,60,70,80,90,100"
cores <- as.numeric(args[4]) # 20
# DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/MARS/refvar/RovHer_chr17_MARS_d3_np150_nv0.1_lowerAF0e+00_annono_embedrefvar_blk28"
# PLINK_list <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/MARS/refvar/RovHer_chr17_MARS_d3_np150_nv0.1_lowerAF0e+00_annono_embedrefvar_blk28/4_CLUMP_RESULT/yhat_01range_CHR_17.clumped"
# prop_blks <- c(1,5,10,15,20,25,30,35,40,50,60,70,80,90)
# cores <- 5

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(doParallel)
  library(parallel)
})
registerDoParallel(makeCluster(cores))

PLINK_list_base <- sub("\\.clumped$", "", PLINK_list)
  
############################# Directories  #######################

# Input dir (clumping dir)
INPUT_DIR <- paste0(DIR,"/8_GENO_0.1LD50kWIN_RDATA")
cat("\nInput dir:", INPUT_DIR, "\n\n")
if (!dir.exists(INPUT_DIR)) {
   stop("ERROR: No Input dir")
}

# Create root output dir 
root <- paste0(DIR, "/4_exome_h2"); if (!dir.exists(root)) { dir.create(root)}

######################################################################
# Iterate over each percentage
######################################################################

# Convert prop_blks to a numeric vector if it is a character
if (is.character(prop_blks)) {
  prop_blks <- as.numeric(unlist(strsplit(prop_blks, ",")))
}

for (top in prop_blks) {
  if (top == 100 || top == "100") {
    stop("Top 100% is not supported. Skip to SCRIPT 3.\n")
  }

  # Set output directory
  root2 <- paste0(root, "/top", top); if (!dir.exists(root2)) { dir.create(root2)}
  GENO_DIR <- paste0(root2, "/2_GENO_0.1LD50kWIN_RDATA")
  if (!dir.exists(GENO_DIR)) {
    dir.create(GENO_DIR, recursive = TRUE)
  } else {
    files_in_dir <- list.files(GENO_DIR, full.names = TRUE)
    if (length(files_in_dir) == 0) {
      cat("GENO_DIR exists but is empty.\n")
    } else {
      raw_data_files <- grep("\\.raw\\.Data$", files_in_dir, value = TRUE)
      if (length(raw_data_files) > 0) {
        stop("GENO_DIR contains files with the '.raw.Data' ending. Exiting the code.")
      }
    }
  }

  # Construct input PLINK_list for the current iteration
  PLINK_list <- paste0(PLINK_list_base, "_top", top, ".clumped")
  
  # Load input file
  var_list <- fread(PLINK_list) # input file
  num_top_variant <- nrow(var_list)

  # Convert var_list to a vector of SNP names
  PLINK_IDs <- var_list$PLINK_SNP_NAME
  required_RV_cols <- c("IID", PLINK_IDs) # Include the "IID" column

  # Initialize 
  temp_GENO <- NULL # empty dt that can be appended to
  col_count <- 0 # count for total number of columns across GENO blks
  new_blks <- 0

  # Get the chromosomes that this RV subset is in
  chromosomes <- sapply(strsplit(as.character(var_list[[1]]), ":"), "[", 1)
  unique_chr <- unique(chromosomes)
  unique_chr <- as.numeric(unique_chr)
  cat(paste0(num_top_variant, " RVS in top ", top, "% within chr: ", unique_chr, "\n"))

  ######################################################################
  #  LOOP: Use PLINK IDs to search each genotype matrix in diretory
  ######################################################################

  start_time <- Sys.time()

  for (chr in unique_chr) { # Searchiing .RData within ./RRV_PLINK_CLUMP/M_CAP_score/7_GENOMATRIX_0.1LD50kWIN_RRV_RDATA

    # List all geno matrices in INPUT_DIR:
    geno_files <- list.files(
      path = INPUT_DIR, 
      pattern = paste0("CHR_", chr, ".*\\.raw\\.RData$"), 
      full.names = TRUE
    )
    geno_files <- sort(geno_files) 

    if (length(geno_files) == 0) {
      stop("No genotype files found for chromosome", chr, "\n")
    } else {
      cat(length(geno_files), "clumped GENO files for chr", chr, "\n")
    }
    
    count <- 0 # Reset count for each chromosome
    
    for (geno_file in geno_files) {
      original_basename <- basename(geno_file)
      new_basename <- sub(
        "clumped_\\d{2}", 
        sprintf("clumped_%02d", as.numeric(count)), 
        original_basename
      )
      out_path <- file.path(GENO_DIR,new_basename)

      cat("\n[SCRIPT 2: Top", top, "%    ", basename(geno_file), "]\n")
      base::load(geno_file) # load object 'GENO'
      cat("Original:", dim(GENE_DF)[1], " x", dim(GENE_DF)[2], "\n")

      # re-write PLINK var column names
      colnames(GENE_DF) <- ifelse(
        colnames(GENE_DF) == "IID", 
        "IID", 
        gsub("_[A-Z]$", "", colnames(GENE_DF))
      )

      # ------------------------------------------------------
      # Search GENE_DF to keep only the variants from PLINK_LIST
      # ------------------------------------------------------
      pattern <- paste0("^", chr, ":")
      required_RV_cols_chr <- required_RV_cols[grepl(pattern, required_RV_cols)]
      required_RV_cols_chr <- c("IID", required_RV_cols_chr)  # Add "IID" to the vector
      common_columns <- intersect(names(GENE_DF), required_RV_cols_chr) # Subset the GENE_DF dataframe based on the relevant variants
      cat("Found", length(common_columns) - 1, "RVs\n")
      # print(head(GENE_DF[,1:5]))
      # cat("\n")

      if (length(common_columns) == 1) { # 1 = only IID column
        rm(GENE_DF)
        next
      }
      GENE_DF <- GENE_DF[, ..common_columns, with=FALSE]
      # Rename IID to eid
      if ("IID" %in% names(GENE_DF)) {
        GENE_DF <- rename(GENE_DF, eid = IID)
      }

      # Append filtered GENE_DF to existing temp_GENO
      if (is.null(temp_GENO)) {# If temp_GENO is empty, bind all columns from GENE_DF
          temp_GENO <- GENE_DF
      } else {                 # If not empty, bind all columns except for "eid"
          temp_GENO <- bind_cols(temp_GENO, GENE_DF[,-1])
          # temp_GENO <- bind_cols(temp_GENO, GENE_DF[ , !colnames(GENE_DF) %in% "eid"])
      }
      rm(GENE_DF); invisible(gc()) 

      # Current column count (Num of RVs) 
      col_count <- col_count + ncol(temp_GENO)
      cat("Variants found so far:", col_count, "\n\n")
      # print(head(temp_GENO[,1:5]))
      # cat("\n")

      # ----------------------------------------------------------------
      # If temp_GENO has >5k cols, save it, or else continue appending 
      # ----------------------------------------------------------------
      if (ncol(temp_GENO) >= 4999) {
        GENE_DF <- temp_GENO
        save(temp_GENO, file = out_path)
        cat("Blk has reached 5k size:", dim(temp_GENO)[1], "x", dim(temp_GENO)[2], "\n")
        cat("Saved:", new_basename, "\n\n")

      # Reinitialize
        temp_GENO <- NULL
        new_blks <- new_blks + 1
      } # ---- Saved blk if it has reached 5k size

    count <- count + 1 
    } # ------ Finished processing all geno files for a chromosome -----

    # Save any remaining data in temp_GENO
    if (!is.null(temp_GENO) && ncol(temp_GENO) > 0) {
      cat("Final temp_GENO:", dim(temp_GENO)[1], "x", dim(temp_GENO)[2], "\n")      
      GENE_DF <- temp_GENO
      print(head(GENE_DF[,1:5]))
      cat("\n")
      base::save(GENE_DF, file = out_path)
      cat("Saved :", out_path, "\n\n")

      # Reinitialize
      temp_GENO <- NULL
      GENE_DF <- NULL
      new_blks <- new_blks + 1
    }
  } 
  # ================== DONE: End of main loop across all Chr ==================
  num_participants <- nrow(temp_GENO) 
  cat("Total RVS found:", col_count, "\n")
}

end_time <- Sys.time()
time_taken <- as.numeric(end_time - start_time, units = "mins")
cat("Time:", round(time_taken, 2), "mins\n\n")

cat("Results:", GENO_DIR, "\n\n")
cat("\n============ SCRIPT 2:Top", top, "% RVs============","\n")
