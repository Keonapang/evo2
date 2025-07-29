# Date: July 27, 2025

# Rscript "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/h2/split_prop_blk.r" $input_file $anno_name $prop_blks $sort $score_file

# This script extracts the top percentages of rows from a clumped file based on a specified annotation column.
# It sorts the rows based on the annotation column in either ascending or descending order and then extracts the top specified percentages.

library(data.table)
suppressMessages(library(doParallel))
suppressMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
input_file <- as.character(args[1]) # 0.00003
anno_name <- as.character(args[2])
prop_blks <- as.character(args[3]) # "1,5,10,15,20,25,30,35,40,50,60,70,80,90,100"
sort <- as.character(args[4]) # "ascending" or "descending"
score_file <- as.character(args[5]) # Not used in the function, but included for compatibility

# check if input files exists
if (!file.exists(score_file)) {
  stop(paste("\n\nscore_file", score_file, "does not exist.\n\n"))
}
if (!file.exists(input_file)) {
  stop(paste("\n\ninput_file", input_file, "does not exist.\n\n"))
}

extract_top_percentages <- function(input_file, score_file, anno_name, prop_blks, sort) {

  df <- fread(input_file, header = TRUE, select = "SNP")
  score_df <- fread(score_file, header = TRUE)

  # rename SNP column from score_df to PLINK_SNP_NAME
  setnames(score_df, old = "SNP", new = "PLINK_SNP_NAME")
  # Check if the column exists
  if (!anno_name %in% colnames(score_df)) {
    stop(paste("Error: Column", anno_name, "not found in score_df."))
  }  
  # Perform a left join based on "PLINK_SNP_NAME"
  colnames(df) <- "PLINK_SNP_NAME"
  new_df <- merge(df, score_df[, .(PLINK_SNP_NAME, get(anno_name))], 
                  by = "PLINK_SNP_NAME", 
                  all.x = TRUE)
  setnames(new_df, old = "V2", new = anno_name)

  # Sort the data frame by the specified column and order
  if (sort == "ascending") {
    new_df <- new_df[order(new_df[[anno_name]]), ]
  } else if (sort == "descending") {
    new_df <- new_df[order(-new_df[[anno_name]]), ]
  } else {
    stop("Error: sort must be either 'ascending' or 'descending'.")
  }
  
  # Total number of rows
  total_rows <- nrow(new_df)

  # remove all columns except PLINK_SNP_NAME and the annotation column
  new_df <- new_df[, c("PLINK_SNP_NAME", anno_name), with = FALSE]
  print(head(new_df,2))
  cat("\n")
  # Convert prop_blks to a numeric vector if it is a character
  if (is.character(prop_blks)) {
    prop_blks <- as.numeric(unlist(strsplit(prop_blks, ",")))
  }
  
  # Process each proportion in prop_blks
  for (prop_blk in prop_blks) {
    # Calculate the number of rows to extract
    num_rows <- ceiling(total_rows * prop_blk / 100)
    
    # Extract the top rows
    top_rows <- new_df[1:num_rows, ]
    
    # Define the export file path
    export_file <- gsub("\\.clumped$", paste0("_top", prop_blk, ".clumped"), input_file)
    write.table(top_rows, file = export_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    cat(paste("Exported top", prop_blk, "%\n"))
  }
  cat(paste("\nLast exported file:", export_file, "\n\n"))
}
extract_top_percentages(input_file, score_file, anno_name, prop_blks, sort)
