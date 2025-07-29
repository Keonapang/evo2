#!/usr/bin/env Rscript
# Date: July 20, 2025
#####################################################################################################################################################
# Descrption: Training one MARS model with a single-trait's FDR using 5-fold CV, then saves all predictions/plots/scores into one folder.
    # - Allows adapting RovHer framework to include MORE predictor columns (i.e. evo2 embeddings etc..)
    # - enables pre-filtering step to remove poorly correlated predictors under a p-value threshold
# ####################################################################################################################################################
# trial_name="chr17" #  B C D E
# d=3
# np=1000 # 50, 20000, Reducing nprune reduces exhaustive search time; maximum number of permissible terms in the final pruned model.
# nv=0.1
# AF=0
# cores=18
# anno_cols="yes"
# embedding_col="refvar" # "delta" , "no" (no embeddings), "refvar"
# blks="28" # (not accounted for if embedding_col="no")
# reverse="no" # "yes or "no" (to include reverse complement embeddings)
# pval_thresh=0.1 # numeric(linreg pre-filtering to keep only sig. predictors to pval < pval_thresh)
# script_path="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/MARS"

# for t in $trial_name; do
# for blk in $blks; do
#   Rscript ${script_path}/train_cv_MARS.r $AF $d $np $nv $t $anno_cols $embedding_col $blk $cores $reverse $pval_thresh $VARWIN
# done
# done

#####################################################################################################################################################
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("earth", quietly = TRUE)) install.packages("earth")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(earth))
suppressMessages(library(doParallel))
suppressMessages(library(parallel))
suppressMessages(library(tibble))
suppressMessages(library(readr))
suppressMessages(library(openxlsx))

args <- commandArgs(trailingOnly = TRUE)
AF <- as.numeric(args[1]) # 0.00003
d <- as.numeric(args[2]) # 2 or 3
np <- as.numeric(args[3]) # 11, 13, 14
nv <- as.numeric(args[4]) # 0.1
trial_name <- args[5] # "A", "B", "C"
anno <- args[6] # "yes", "no"
embedding_col <- as.character(args[7]) # "delta" , "no" (no embeddings), "refvar"
blk <- as.character(args[8])
cores <- as.numeric(args[9]) # "yes", "no"
reverse <- as.character(args[10]) # "yes or "no" (to include reverse complement embeddings)
pval_thresh <- as.numeric(args[11])  # "yes or "no" (linreg pre-filtering to keep only sig. predictors to pval < pval_thresh)
VAR_WIN <- as.character(args[12]) # "1", "128" or "256"

# d <- 3
# AF <- 0 # Lower bound MAF threshold
# np <- 100
# nv <- 0.1
# trial_name <- "chr17" # "A", "B", "C", "D", "E"
# anno <- "yes" # "yes", "no"
# embedding_col <- "refvar"
# reverse <- "no" # "yes or "no" (to include reverse complement embeddings)
# blk <- "28" # "28", "29", "30", "31", "32", "33", "34", "35"
# VAR_WIN <- "1" # "1", "128" or "256"
# pval_thresh <- 0.01 # numeric(linreg pre-filtering to keep only sig. predictors to pval < pval_thresh)
# cores <- 30 # Number of cores to use for parallel processing

# date <- 20250723
# date <- format(Sys.Date(), "%Y%m%d")
# cat("date:", date, "\n")

registerDoParallel(makeCluster(cores))
upper_MAF <- 0.01
p <-2
t <- 0.00001
nk <- 400
if (np == "NULL") {
  np <- NULL
} else {
  np <- as.numeric(np)
}
upper_MAF <- format(as.numeric(upper_MAF), scientific = TRUE, digits = 4)
AF <- format(as.numeric(AF), scientific = TRUE, digits = 4)

# Output model name
if (embedding_col == "delta" || embedding_col == "refvar") {
    model_name <- paste0(trial_name, "_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_blk", blk, "_VARWIN", VAR_WIN) 
    if (reverse == "yes") {
        model_name <- paste0(model_name, "_rev")
    }
} else if (embedding_col == "no") {
    model_name <- paste0(trial_name, "_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_embed", embedding_col) 
} else {
    stop("Invalid 'embedding_col'. Must be either: delta, refvar, no")
}
cat("\nModel name: ", model_name, "\n\n")

# Input directory 
server <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk"
dir <- paste0(server,"/Training/height_FDR/GENEBASS_RV_train_RV_pred_Oct10") # old folder (no cross validation): allvar_predict_saved_models_NO_GENESETS 
DIR_SCORES <- paste0(server,"/evo2/June30")
DIR_EMBEDDING <- paste0(server,"/evo2/July20_embedding")

# Output directory
ROOTDIR <- paste0(server,"/evo2/MARS/RovHer") # <------ MODIFY!
OUTDIR1 <- paste0(ROOTDIR,"/", trial_name)
OUTDIR2 <- paste0(OUTDIR1,"/pval_", pval_thresh)
OUTDIR3 <- paste0(OUTDIR2, "/embedding_", embedding_col)
OUTDIR <- paste0(OUTDIR2,"/", model_name)
if (dir.exists(OUTDIR)) {
    if (length(list.files(OUTDIR)) > 0) {
        stop(paste0("Output dir ", OUTDIR, " already exists and is not empty.\n"), call. = FALSE)
    } else {
        cat(paste0("Output dir ", OUTDIR, " exists and is empty.\nProceeding...\n\n"))
    }
} else {
    cat(paste0("Output dir ", OUTDIR, " does not exist.\nProceeding...\n\n"))
}

if (!dir.exists(ROOTDIR)) {dir.create(ROOTDIR, recursive = TRUE)}
if (!dir.exists(OUTDIR1)) {dir.create(OUTDIR1, recursive = TRUE)}
if (!dir.exists(OUTDIR2)) {dir.create(OUTDIR2, recursive = TRUE)}
if (!dir.exists(OUTDIR3)) {dir.create(OUTDIR3, recursive = TRUE)}
if (!dir.exists(OUTDIR)) {dir.create(OUTDIR, recursive = TRUE)}

# Output files
pval_file <- paste0(OUTDIR2, "/sig_pval.txt") # only if pval_thresh is set
RData_file <- paste0(OUTDIR, "/", model_name, "_scores.RData")
scores_file <- paste0(OUTDIR, "/", model_name, "_scores.txt")
metrics_file <- paste0(OUTDIR, "/layer", blk, "_metrics.txt")

####################################################################################
# INPUT FILES
####################################################################################

# Script
pval_script <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Keona_scripts/linreg_pval_filter.r"

# Functional annotations
train_path <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Training/height_FDR/genebass_AF001_chrall_vars_gene_final.txt" # 4926639 x  80 
predict_path <- train_path

# Evo2 scores
file_score <- paste0(DIR_SCORES, "/evo2_chr17_win4096.txt")

# Evo2 embeddings
file_delta <- paste0(DIR_EMBEDDING, "/RovHer_", trial_name, "_blocks.", blk, ".mlp.l3_VARWIN",VAR_WIN, "_delta.csv")
file_var <- paste0(DIR_EMBEDDING, "/RovHer_", trial_name, "_blocks.", blk, ".mlp.l3_VARWIN",VAR_WIN, "_var.csv")
file_ref <- paste0(DIR_EMBEDDING, "/RovHer_", trial_name, "_blocks.", blk, ".mlp.l3_VARWIN",VAR_WIN, "_ref.csv")

# Reverse complement embeddings
file_delta_rev <- paste0(DIR_EMBEDDING, "/RovHer_", trial_name, "_blocks.", blk, ".mlp.l3_VARWIN",VAR_WIN, "_delta_rev.csv")
file_var_rev <- paste0(DIR_EMBEDDING, "/RovHer_", trial_name, "_blocks.", blk, ".mlp.l3_VARWIN",VAR_WIN, "_var_rev.csv")
file_ref_rev <- paste0(DIR_EMBEDDING, "/RovHer_", trial_name, "_blocks.", blk, ".mlp.l3_VARWIN",VAR_WIN, "_ref_rev.csv")

# all_data <- rbindlist(
#   lapply(1:6, function(i) {
#     file_embed <- paste0(DIR_EMBEDDING, "/RovHer_", trial_name, "_", i, "_blocks.", blk, ".mlp.l3_delta.csv")
#         fread(file_embed)}),use.names = TRUE,fill = TRUE)
# output_file <- paste0(DIR_EMBEDDING, "/RovHer_", trial_name, "_blocks.", blk, ".mlp.l3_delta.csv")
# fwrite(all_data, output_file)
# cat("All files have been combined and saved to:", output_file, "\n")

####################################################################################
# Load 75 functional annotations
####################################################################################
train_df <- as.data.frame(fread(train_path, header = TRUE))
cat("train_df: ", dim(train_df)[1], "RVs x", dim(train_df)[2], "\n\n")

# Evo2 scores
evo2_scores <- as.data.frame(fread(file_score, header = TRUE, 
                            select = c("PLINK_SNP_NAME", "evo2_delta_score"),
                            showProgress = FALSE))
####################################################################################
# Load embeddings
####################################################################################
if (embedding_col == "delta") {
  delta_embed <- fread(file_delta)
  delta_embed <- delta_embed %>% select(-layer, -input_file, -RovHer_score)
  cat("delta_embed:", dim(delta_embed)[1], "x", dim(delta_embed)[2], "\n\n")
} else if (embedding_col == "refvar") {
  var_embed <- fread(file_var)
  ref_embed <- fread(file_ref)
  var_embed <- var_embed %>% select(-layer, -input_file, -RovHer_score)
  ref_embed <- ref_embed %>% select(-layer, -input_file, -RovHer_score)
  cat("var_embed:",dim(var_embed)[1],"x",dim(var_embed)[2], "  ref_embed:",dim(ref_embed)[1],"x",dim(ref_embed)[2],"\n")
  if (nrow(var_embed) != nrow(ref_embed)) {
      stop("Number of rows in delta, var, and ref embeddings do not match.")
  }
}

# reverse complement embeddings
if (reverse == "yes") {
  cat("Adding reverse complement strand...\n")

  if (embedding_col == "delta") {
    delta_embed <- fread(file_delta_rev)
    delta_embed <- delta_embed %>% select(-layer, -input_file, -RovHer_score)

  } else if (embedding_col == "refvar") {
    var_embed <- fread(file_var_rev)
    ref_embed <- fread(file_ref_rev)
    var_embed <- var_embed %>% select(-layer, -input_file, -RovHer_score)
    ref_embed <- ref_embed %>% select(-layer, -input_file, -RovHer_score)

    if (nrow(var_embed) != nrow(ref_embed)) {
        stop("Number of rows in delta, var, and ref embeddings do not match.")
    }
  }
  colnames(csv1_df_rev) <- ifelse(grepl("^e[0-9]+$", colnames(csv1_df_rev)), paste0(colnames(csv1_df_rev), "_rev_delta"), colnames(csv1_df_rev))
  colnames(csv1ref_df_rev) <- ifelse(grepl("^e[0-9]+$", colnames(csv1ref_df_rev)), paste0(colnames(csv1ref_df_rev), "_rev_ref"), colnames(csv1ref_df_rev))
  colnames(csv1var_df_rev) <- ifelse(grepl("^e[0-9]+$", colnames(csv1var_df_rev)), paste0(colnames(csv1var_df_rev), "_rev_var"), colnames(csv1var_df_rev))

  # Remove columns "input_file", "layer", "PLINK_SNP_NAME", "RovHer_score" from csv1_df_rev
  columns_to_remove <- c("input_file", "layer", "PLINK_SNP_NAME", "RovHer_score")
  csv1_df_rev <- csv1_df_rev[, !..columns_to_remove, with = FALSE] 
  csv1ref_df_rev <- csv1ref_df_rev[, !..columns_to_remove, with = FALSE]
  csv1var_df_rev <- csv1var_df_rev[, !..columns_to_remove, with = FALSE]
  cat("csv1_df_rev:", dim(csv1_df_rev), "\n")

  # join csv1_df_rev to the end of csv1_df column-wise
  csv1_df <- cbind(csv1_df, csv1_df_rev)
  csv1ref_df <- cbind(csv1ref_df, csv1ref_df_rev)
  csv1var_df <- cbind(csv1var_df, csv1var_df_rev)
  cat("csv1_df + reverse complement:", dim(csv1_df), "\n")
}

# Keep only certain chromosomes RVs
train_df <- train_df %>% filter(grepl("^17:", PLINK_SNP_NAME)) 
cat("train_df Chr 17-only: ", dim(train_df)[1], "RVs x", dim(train_df)[2], "cols\n\n") # 'df' of size 3 million x 97 cols

####################################################################################
# If p-val threshold is set, do pre-filtering of predictor cols using linear regression
####################################################################################

if (pval_thresh != "NA" && pval_thresh > 0) {

  cat("Pre-filtering predictors to keep only p-values < ", pval_thresh, "\n")
  embedding_df <- train_df %>% select(PLINK_SNP_NAME, height_FDR, GENEBASS_AF)

  if (embedding_col == "delta") {
    embedding_df <- merge(embedding_df, delta_embed, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join
  } else if (embedding_col == "refvar") {
    colnames(var_embed) <-ifelse(grepl("^e[0-9]+$", colnames(var_embed)), paste0(colnames(var_embed), "_var"), colnames(var_embed))
    colnames(ref_embed) <-ifelse(grepl("^e[0-9]+$", colnames(ref_embed)), paste0(colnames(ref_embed), "_ref"),  colnames(ref_embed))
    embedding_df <- merge(embedding_df, var_embed, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join
    embedding_df <- merge(embedding_df, ref_embed, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join
  } else if (embedding_col == "no") {
    cat("No embeddings used.\n")
  } else {
    stop("Invalid 'embedding_col'. Must be either: delta, refvar, no")
  }

  if (!file.exists(pval_file)) {
    cat("Performing linear regression on each predictor col...\n")

    results <- data.frame(Predictor = character(), P_value = numeric(), stringsAsFactors = FALSE)
    exclude_cols <- c("PLINK_SNP_NAME", "height_FDR", "GENEBASS_AF")
    predictor_cols <- setdiff(colnames(embedding_df), exclude_cols)

    for (col in predictor_cols) {
      cat("Processing:", col, "\n")
      formula <- as.formula(paste("height_FDR ~", col))
      model <- lm(formula, data = embedding_df)
      p_value <- summary(model)$coefficients[2, 4]  # 2nd row corresponds to the predictor
      results <- rbind(results, data.frame(Predictor = col, P_value = p_value, stringsAsFactors = FALSE))
    }
    cat("Subset results to keep only significant p-value...\n")
    sig_results <- results[results$P_value < pval_thresh, ]
    cat("All predictors:", dim(results)[1],"     Sig. predictors:", dim(sig_results), "\n\n")
    fwrite(sig_results, file = pval_file, sep = "\t", quote = FALSE, row.names = FALSE)

    # Filter to keep only significant predictors
    embed_cols_to_keep <- sig_results$Predictor
    embedding_df <- embedding_df %>%
                    select(PLINK_SNP_NAME, height_FDR, GENEBASS_AF, all_of(embed_cols_to_keep))
    cat("Significant embedding predictors: ", dim(embedding_df)[1], "x", dim(embedding_df)[2], "\n\n")

  } else { # If significant predictors file exists already
      cat("Sig. predictors file already exists. Loading it in....\n")
      sig_results <- fread(pval_file, header = TRUE)
      embed_cols_to_keep <- sig_results$Predictor
      embedding_df <- embedding_df %>%
                    select(PLINK_SNP_NAME, height_FDR, GENEBASS_AF, all_of(embed_cols_to_keep))
    cat("Significant embedding: ", dim(embedding_df)[1], "x", dim(embedding_df)[2], "\n\n")
  }
} else { # if no p-value thresholding is needed, we just use all predictors

  if (embedding_col == "delta") {
      embedding_df <- merge(embedding_df, delta_embed, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join
  
  } else if (embedding_col == "refvar") {
      
      colnames(var_embed) <-ifelse(grepl("^e[0-9]+$", colnames(var_embed)), paste0(colnames(var_embed), "_var"), colnames(var_embed))
      colnames(ref_embed) <-ifelse(grepl("^e[0-9]+$", colnames(ref_embed)), paste0(colnames(ref_embed), "_ref"),  colnames(ref_embed))
      embedding_df <- merge(embedding_df, var_embed, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join
      embedding_df <- merge(embedding_df, ref_embed, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join
    
  } else if (embedding_col == "no") {
      cat("No embeddings used.\n")
  } else {
      stop("Invalid 'embedding_col'. Must be either: delta, refvar, no")
  }
}

####################################################################################
# If anno == "yes", merge embeddings with 75 functional annotations,
# then filter some columns and store 'extra_cols' for later use 
####################################################################################
if (anno == "yes") {
  if ("GENEBASS_AF" %in% colnames(embedding_df)) {
    embedding_df <- embedding_df %>% select(-height_FDR, -GENEBASS_AF)
  }
  train_df <- merge(train_df, embedding_df, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join
  train_df <- merge(train_df, evo2_scores, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join
  extra_cols <- train_df %>% select(PLINK_SNP_NAME, GENEBASS_AF, LOF_DISRUPTIVE, `Gene.refGene`)
  df <- dplyr::select(train_df, -PLINK_SNP_NAME, -`Gene.refGene`) 
} else {
  train_df <- embedding_df
  extra_cols <- train_df %>% select(PLINK_SNP_NAME)
  df <- dplyr::select(train_df, -PLINK_SNP_NAME) 
}
####################################################################################
# Filter training matrix
####################################################################################

# Keep RVs with MAF above the lower bound threshold 
count_rows <- df %>% filter(GENEBASS_AF < as.numeric(AF)) %>% nrow()
cat("Number of GENEBASS RVs with MAF <", AF, "(lower thresh):", count_rows, "\n")
df1 <- subset(df, GENEBASS_AF >= as.numeric(AF)) # Lower bound (4844229 - 2682377) = 2114607 x 76 
df1 <- subset(df1, GENEBASS_AF < upper_MAF) # Upper-bound 

cat("\ndf1:\n")
print(head(df1[, 1:4]),2)
cat("\n")

if (!"height_FDR" %in% colnames(df1)) {
    stop("Column 'height_FDR' not found in df1. Please check.")
}
X_df1 <- df1[, !colnames(df1) %in% "height_FDR"]
y_df1 = df1[["height_FDR"]]
cat(" X matrix after MAF >", AF, "removed:", dim(X_df1)[1], "x", dim(X_df1)[2], "\n")
print(head(X_df1[, 1:4]),2)
cat("\n y FDR:", length(y_df1), "x 1\n\n")  
print(head(y_df1),3)
cat("\n")

####################################################################################
# Load test matrix 
####################################################################################
if (!file.exists(RData_file)) {

# df <- as.data.frame(fread(predict_path, header = TRUE, showProgress = FALSE))
df_pred <- df
X_df_pred1 = df_pred[, -1] # 3mill x 92
y = df_pred[["height_FDR"]]

set.seed(42) # Set seed for reproducibility
k <- 5

# Assign random folds
X_df_pred1$fold <- sample(rep(1:k, length.out = nrow(X_df_pred1)))

# Initialize storage for predictions and models
all_predictions <- numeric(nrow(X_df_pred1))
cv_mars_list <- list() # To store models for each fold
start_time <- Sys.time()

# 5-fold cross-validation loop
cat("Training MARS for ", trial_name, " with anno=", anno, " (d:", d, " np:", np," p:", p, ")\n")
cat("pval_thresh :", pval_thresh, "\n")
for (i in 1:k) {
    cat("=============== [ CV Fold", i, "of", k, "] =================\n\n")

    # Split data into training and testing sets
    train_indices <- which(X_df_pred1$fold != i) # Rows not in the current fold
    test_indices <- which(X_df_pred1$fold == i) # Rows in the current fold

    df2_train <- X_df_pred1[train_indices, ]
    df2_test <- X_df_pred1[test_indices, ]

    # Extract X (features) and y (target) for training
    X_train <- df2_train[, -1] # Exclude the first column (assumed to be the dependent variable)
    y_train <- df2_train[, 1]  # The first column is the dependent variable

    # Get prediction data for the current fold
    pred_test_indices <- which(X_df_pred1$fold == i)
    X_pred <- X_df_pred1[pred_test_indices, ] # Features for prediction

    # Display fold information
    cat("Fold", i, "Training size:", length(train_indices), "rows\n")
    cat("Fold", i, "Testing size:", length(test_indices), "rows\n")
    cat("Fold", i, "Prediction size:", length(pred_test_indices), "rows\n\n")

    # Train the MARS model for the current fold
    cat("Training on X_train:", dim(X_train), "...\n")
    cv_mars <- earth(x = X_train, y = y_train, degree = d, penalty = p, nprune = np, thresh = t, newvar.penalty = nv, trace = 1)

    # Store the model
    cv_mars_list[[paste0("cv_mars_", i)]] <- cv_mars

    # Predict on unseen data (X_pred)
    cat("\nPredicting on X_pred:", dim(X_df_pred1), "...\n")
    predictions <- predict(cv_mars, newdata = X_pred)

    # Store predictions in the corresponding indices
    all_predictions[pred_test_indices] <- predictions
    cat("Stored", length(predictions), "predictions for fold", i, "\n\n")
}
# Convert the predictions from the test sets to a dataframe
yhat_dataframe <- data.frame(yhat = all_predictions) # 155069 x 1
yhat_dataframe <- as.data.frame(yhat_dataframe)
final_df <- cbind(extra_cols, yhat_dataframe) # 3128762 x 4
save(final_df, cv_mars_list, d, p, t, np, nv, train_path, predict_path, file = RData_file)
fwrite(final_df, file=scores_file, sep = "\t", quote = FALSE)

# Calculate duration
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")
cat("------- Duration: ", round(duration, 2), " mins -------\n")
} else {
    cat("\nRData file already exists. Loading it in...\n")
}

####################################################################################
# Plot variables
####################################################################################
font <- "Sans"
psize <- 15
red <- "#DE433C"
green <- "#76B900"
font <- "Open Sans"
royal_blue <- "#0B7C9A"
nice_red <- "#CC6a66"
yellow <- "#F2C46D"

plot_theme <- 
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
  axis.title = element_text(size = 14, face = "bold", family = font),  # Increase axis title size
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = font),  # Centered, bold title
  legend.position = "none"  # Remove the legend
)


####################################################################################
# PLOT: Variable Importance Plot
####################################################################################

load(RData_file)

if (exists("cv_mars_list")) {
    cat("\nSummarizing cross-validation models:\n")
    
    for (i in 1:length(cv_mars_list)) {
    
        cat("========= Trial", trial_name, " model ", i, "=========\n")
        model <- cv_mars_list[[i]] 
        variable_importance <- evimp(model)
        preds_list <- rownames(variable_importance)
        cat("Number of sig. predictors:", length(preds_list), "\n\n")
        gcv_list <- variable_importance[, "gcv"]

        plot_data <- data.frame(Predictor = names(gcv_list),GCV = as.numeric(gcv_list))

        #Sort the data frame by GCV in descending order
        plot_data <- plot_data[order(plot_data$GCV, decreasing = TRUE), ]
        
        # Remove backticks from predictor names
        plot_data$Predictor <- as.factor(gsub("`", "", as.character(plot_data$Predictor)))
        # Remove row if Predictor is "fold"
        plot_data <- plot_data[plot_data$Predictor != "fold", ]

        # Prepare predictor categories 
        anno_cols <- train_df
        gene_cols <- colnames(anno_cols)[(48 + 1):ncol(anno_cols)]
        variant_cols <- colnames(anno_cols)[1:48]

        plot_data$Category <- ifelse(
            grepl("^e", plot_data$Predictor),  # Predictors starting with "e_"
            "Embeddings (Evo2 7B)",
            ifelse(
                plot_data$Predictor %in% gene_cols,    # Predictors in gene_cols
                "Gene-level annotations",
                ifelse(
                    plot_data$Predictor %in% variant_cols,  # Predictors in variant_cols
                    "Variant-level annotations",
                    "Other"  # Fallback category if none of the above apply
                )
            )
        )

        # remove "_score" suffix from predictor names
        plot_data$Predictor <- as.factor(sub("_score$", "", as.character(plot_data$Predictor)))
        plot_data$Predictor <- gsub("PrimateAI_3D", "PrimateAI-3D", plot_data$Predictor)
        plot_data$Predictor <- gsub("BayesDel_addAF_rankscore", "BayesDel-addAF", plot_data$Predictor)
        plot_data$Predictor <- gsub("M-CAP_score", "M-CAP", plot_data$Predictor)
        plot_data$Predictor <- gsub("fathmm-XF_coding_score", "fathmm-XF coding", plot_data$Predictor)

        plot_data$Predictor <- factor(plot_data$Predictor, levels = plot_data$Predictor[order(plot_data$GCV, decreasing = TRUE)])
        p <- ggplot(plot_data, aes(x = GCV, y = Predictor, fill = Category)) +
        geom_bar(stat = "identity", color = "black", width = 0.7) +  # Horizontal bars
        labs(
            title = paste0("Variable Importance Plot for Fold ", i),
            x = "Generalized Cross-Validation (GCV)",
            y = "Predictors",
            fill = "Annotation type"  # Legend title
        ) +
        scale_fill_manual(
            values = c(
            "Embeddings (Evo2 7B)" = "#69b3a2", # green
            "Gene-level annotations" = "#F2C46D", # yellow
            "Variant-level annotations" = "#CC6a66", # red
            "Other" = "#A9A9A9"  # gray for other categories
            )
        ) +
        theme_minimal() +
        theme(
        axis.title.x = element_text(margin = margin(t = 15, b = 0.8), family = font),  # Padding for x-axis
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0.8), family = font),  # Padding for y-axis
        axis.text = element_text(size = 14, family = font),  # Increase axis text size
        axis.title = element_text(size = 14, face = "bold", family = font),  # Increase axis title size
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = font),
            panel.background = element_rect(fill = "white", color = NA),  # White background
            plot.background = element_rect(fill = "white", color = NA),   # White plot area
            legend.title = element_text(size = 14, family = font),
            legend.text = element_text(size = 12, family = font),
            plot.margin = margin(t = 15, r = 15, b = 20, l = 20),         # Adjust plot margins
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank()
        )
        plot_file <- paste0(OUTDIR, "/", model_name, "_gcv_",i,".png")
        
        if (anno == "no") {
          ggsave(plot_file, plot = p, width = 14, height = 18, dpi = 300)
        } else {
          ggsave(plot_file, plot = p, width = 14, height = 8, dpi = 300)
        }
    }
} else {
    variable_importance <- evimp(mars_model)
    preds_list <- rownames(variable_importance)
    cat("Number of sig. predictors:", length(preds_list), "\n")
    gcv_list <- variable_importance[, "gcv"]

    # Prepare the data
    plot_data <- data.frame(Predictor = names(gcv_list),GCV = as.numeric(gcv_list))
    plot_data <- plot_data[order(plot_data$GCV, decreasing = TRUE), ]

    # Remove backticks from predictor names
    plot_data$Predictor <- as.factor(gsub("`", "", as.character(plot_data$Predictor)))
    # Remove row if Predictor is "fold"
    plot_data <- plot_data[plot_data$Predictor != "fold", ]

    # Prepare predictor categories 
    anno_cols <- train_df
    gene_cols <- colnames(anno_cols)[(48 + 1):ncol(anno_cols)]
    variant_cols <- colnames(anno_cols)[1:48]

    plot_data$Category <- ifelse(
        grepl("^e", plot_data$Predictor),  # Predictors starting with "e" and ending with numbers
        "Embeddings (Evo2 7B)",
        ifelse(
            plot_data$Predictor %in% gene_cols,    # Predictors in gene_cols
            "Gene-level annotations",
            ifelse(
            plot_data$Predictor %in% variant_cols,  # Predictors in variant_cols
            "Variant-level annotations",
            "Other"  # Fallback category if none of the above apply
            )
        )
    )
      # remove "_score" suffix from predictor names
      plot_data$Predictor <- as.factor(sub("_score$", "", as.character(plot_data$Predictor)))
      plot_data$Predictor <- gsub("PrimateAI_3D", "PrimateAI-3D", plot_data$Predictor)
      plot_data$Predictor <- gsub("BayesDel_addAF_rankscore", "BayesDel-addAF", plot_data$Predictor)
      plot_data$Predictor <- gsub("M-CAP_score", "M-CAP", plot_data$Predictor)
      plot_data$Predictor <- gsub("fathmm-XF_coding_score", "fathmm-XF coding", plot_data$Predictor)

    # Ensure predictors are ordered by GCV (descending)
    plot_data$Predictor <- factor(plot_data$Predictor, levels = plot_data$Predictor[order(plot_data$GCV, decreasing = TRUE)])

    p <- ggplot(plot_data, aes(x = GCV, y = Predictor, fill = Category)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +  # Horizontal bars
    labs(
        title = "Variable Importance Plot",
        x = "Generalized Cross-Validation (GCV)",
        y = "Predictors",
        fill = "Annotation type"  # Legend title
    ) +
    scale_fill_manual(
            values = c(
            "Embeddings (Evo2 7B)" = "#69b3a2", # green
            "Gene-level annotations" = "#F2C46D", # yellow
            "Variant-level annotations" = "#CC6a66", # red
            "Other" = "#A9A9A9"  # gray for other categories
            )
    ) +
    theme_minimal() +
    theme(
    axis.title.x = element_text(margin = margin(t = 15, b = 0.8), family = font),  # Padding for x-axis
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0.8), family = font),  # Padding for y-axis
    axis.text = element_text(size = 14, family = font),  # Increase axis text size
    axis.title = element_text(size = 14, face = "bold", family = font),  # Increase axis title size
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = font),
        panel.background = element_rect(fill = "white", color = NA),  # White background
        plot.background = element_rect(fill = "white", color = NA),   # White plot area
        legend.title = element_text(size = 14, family = font),
        legend.text = element_text(size = 12, family = font),
        plot.margin = margin(t = 15, r = 15, b = 20, l = 20),         # Adjust plot margins
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
    )
    plot_file <- paste0(OUTDIR, "/", model_name, "_gcv.png")
        
    if (anno == "no") {
          ggsave(plot_file, plot = p, width = 14, height = 18, dpi = 300)
    } else {
          ggsave(plot_file, plot = p, width = 14, height = 8, dpi = 300)
    }
}

####################################################################################
# PLOT: Plot distribution of scores
####################################################################################

final_df <- fread(scores_file)
colnames(final_df)[ncol(final_df)] <- "y_df1"

p <- ggplot() +
  geom_density(
    data = final_df,
    aes(x = y_df1, y = after_stat(scaled), color = "#C3535E", fill = "#C3535E"), 
    alpha = 0.3, size = 1, adjust = 1.5
  ) +
  labs(x = "MARS scores", y = "Proportion", fill = "#C3535E") +
  
  # Adjust x-axis limits and breaks to fit y_df1's range
  scale_x_continuous(
    limits = c(-0.5, 1),  # Adjust based on the range of y_df1
    breaks = seq(-0.5, 1, by = 0.25), 
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 20, margin = margin(t = 10), family = font),
    axis.text.y = element_text(size = 20, margin = margin(l = 15, r = 3), family = font),
    axis.title.x = element_text(margin = margin(t = 20, b = 5), family = font),  # Further from axis
    axis.title.y = element_text(margin = margin(r = 20, l = 20), family = font),  # Further from axis
    axis.ticks.length = unit(0.5, "cm"),  # Longer ticks
    axis.ticks = element_line(size = 1.3),  # More visible ticks
    axis.line = element_line(size = 1.3),  # Increase thickness of axis lines
    plot.title = element_text(size = 20, hjust = 0, face = "bold", family = font),  # Center-align the plot title
    legend.spacing.y = unit(0.6, "cm"),  # Increase space between legend rows
    legend.spacing.x = unit(0.6, 'cm'), # Increase distance between the two legends
    legend.key.spacing.y = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'), # Increase width of the legend keys
    legend.key.height = unit(0.5, 'cm'), # Increase height of the legend keys
    legend.background = element_rect(fill = "white", colour = "white"),
    legend.position = c(0.05, 0.9),  # Top-left corner
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 20, margin = margin(b = 4, t = 4, l = 4, r = 4), family = font),  # Increase legend text size
    legend.title = element_blank()  # Remove legend title
  )
plot_file <- paste0(OUTDIR, "/", model_name, "_scores.png")
ggsave(plot_file, plot = p, width = 17, height = 6.5, dpi = 300)
cat("File:", plot_file, "\n")


# #################################################################################################
# # Plot ROC curve and calculate AUC
# #################################################################################################
# suppressMessages(library(pROC))
# final_df <- fread(scores_file)
# # rename last column to y_df1
# colnames(final_df)[ncol(final_df)] <- "y_df1"

# # Remove rows where clinvar_clnsig is 0.5
# final_df <- final_df[final_df$clinvar_clnsig != 0.5, ]

# # Step 1: Prepare ground truth as binary (P/LP = 1, B/LB = 0)
# final_df$binary_truth <- ifelse(final_df$clinvar_clnsig %in% c(1, 0.75), 1, 0)

# # Verify the results
# table(final_df$binary_truth)  # Check the counts of the binary classes
# head(final_df)  # Check the modified dataframe

# # Step 2: Calculate AUROC
# roc_obj <- roc(final_df$binary_truth, final_df$y_df1)
# auc_value <- auc(roc_obj)
# cat("AUC:", auc_value, "\n")

# # Step 3: Extract data for plotting the ROC curve
# roc_data <- data.frame(
#   TPR = rev(roc_obj$sensitivities), # True Positive Rate (Sensitivity)
#   FPR = rev(1 - roc_obj$specificities) # False Positive Rate (1 - Specificity)
# )

# auc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
#   geom_line(color = "black", size = 1.2) +  # ROC curve
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = red) +  # Random chance line
#   annotate(
#         "label",  # Use "label" instead of "text" to add a text box with a border
#         x = 0.75, y = 0.2, 
#         label = paste("AUC =", round(auc_value, 3)), 
#         size = 5, 
#         color = "black", 
#         fill = "white",  # Background color of the label
#         label.size = 0.2  # Border thickness (thin black border)
#     ) +
#     labs(
#     title = paste0("ROC Curve for ", test_gene, "  (MARS model: ", trial_name, ")"),
#     x = "False Positive Rate", #  (1 - Specificity)
#     y = "True Positive Rate" # (Sensitivity)
#   ) +
#   theme_minimal() +
#   plot_theme
# output_file <- paste0(OUTDIR,"/", model_name, "_auc.png")
# ggsave(output_file, auc_plot, width = 10, height = 6, dpi = 300)
# cat("Plot saved to:", output_file, "\n")


# ####################################################################################
# # PLOT: Distribution of Evo2 scores of two ClinVar classes (P and B)
# ####################################################################################
# final_df <- fread(scores_file)
# # rename last column to y_df1
# colnames(final_df)[ncol(final_df)] <- "y_df1"

# # Recode clinvar_clnsig into descriptive categories
# final_df$clinvar_clnsig <- factor(
#   final_df$clinvar_clnsig,
#   levels = c(0, 0.25, 0.5, 0.75, 1),  # Original numerical values
#   labels = c("Benign", "Likely Benign", "VUS", "Likely Pathogenic", "Pathogenic")  # New labels
# )

# table(final_df$clinvar_clnsig)
# head(final_df) #  116   5

# # Filter the dataframe to include only "Pathogenic" and "Benign" classes
# filtered_df <- final_df[final_df$clinvar_clnsig %in% c("Pathogenic", "Likely Pathogenic", "Likely Benign", "Benign"), ]

# red <- "#DE433C"
# green <- "#76B900"

# plot <- ggplot(filtered_df, aes(x = y_df1, y = clinvar_clnsig, color = clinvar_clnsig)) +
#   geom_jitter(height = 0.2, width = 0, size = 2, alpha = 0.8) +  # Jittered points
#   geom_boxplot(aes(group = clinvar_clnsig), width = 0.4, alpha = 0.2, outlier.shape = NA, color = "black") +  # Boxplot
#   scale_color_manual(
#     values = c(
#       "Pathogenic" = red,
#       "Likely Pathogenic" = red,
#       "Likely Benign" = green,
#       "Benign" = green
#     )
#   ) +
#   scale_x_continuous(breaks = seq(-1, 1.5, by = 0.25)) +  # Dynamically set x-axis ticks at every 0.25
#   labs(
#     title = paste0("Distribution of pathogenic and benign RVs in ", test_gene),
#     x = paste0("MARS scores (trial: ", trial_name, ")"), 
#     y = "ClinVar classification"
#   ) +
#   theme_minimal() +  # Minimal theme
#   plot_theme
# output_file <- paste0(OUTDIR,"/", model_name, "_dot.png")
# ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
# cat("Plot saved to:", output_file, "\n")

###############################################################################################
# write README file (model summary and variables chosen)
###############################################################################################

log_file_path <- paste0(OUTDIR, "/", model_name, "_log.txt")
if (file.exists(RData_file)) load(RData_file) else cat("ERROR: MODEL", RData_file, "NOT FOUND.\n")
cat("model_file: ", RData_file, "\n\n")
final_df <- fread(scores_file)
colnames(final_df)[ncol(final_df)] <- "y_df1"

# Open the file for writing
sink(log_file_path)

cat("Model:", model_name, "\n")
cat(paste0("Directory:", OUTDIR), "\n")
cat("Training:",train_path,"\n")
cat("Predict:",train_path, "\n\n")
cat(nrow(final_df), "scores (AF < 0.01) for genefbass variants\n")

if (exists("cv_mars_list")) {
    cat("\nSummarizing cross-validation models:\n")
    models <- c("cv_mars_1", "cv_mars_2", "cv_mars_3", "cv_mars_4", "cv_mars_5")
    for (i in 1:length(models)) {
    cat("=======================================", models[i], "==============================================\n")
    model <- cv_mars_list[[i]] # Access each model in the list
    print(summary(model))
    cat("\n")
    print(evimp(model))
    cat("\n")
    cat(format(model, style="C")) # new plot
    cat("\n")
    cat("Min score:", min(final_df$y_df1, na.rm = TRUE), "\n")
    cat("Max score:", max(final_df$y_df1, na.rm = TRUE), "\n\n")
    }
} else {
    cat("=======================================", "no CV", "==============================================\n")
    print(summary(mars_model))
    cat("\n")
    print(evimp(mars_model))
    cat("\n")
    cat(format(mars_model, style="C")) # new plot
    cat("\n")
    cat("Min score:", min(final_df$y_df1, na.rm = TRUE), "\n")
    cat("Max score:", max(final_df$y_df1, na.rm = TRUE), "\n\n")
}
sink() # Close the file
cat("\nLog:", log_file_path, "\n\n")


# ###############################################################################################
# # Write metrics file in excel to save AUROC
# ###############################################################################################

# initialize_metrics_file <- function(file_path) {
#   metrics_csv <- file_path
#     columns <- c(
#     "MODEL_SIZE", "input_file", "subset_method", "USE_RANDOM_SEED", 
#     "SEQ_LENGTH", "WINDOW_SIZE", "time_score_ref", "time_score_var", "AUROC"
#   )
#     if (!file.exists(metrics_csv)) {
#     df <- tibble::tibble(!!!setNames(vector("list", length(columns)), columns))
#     readr::write_csv(df, metrics_csv)
#   }
#   return(metrics_csv)
# }

# append_metrics_row <- function(metrics_file, row_data) {
#   df_csv <- readr::read_csv(metrics_file, show_col_types = FALSE)
#   new_row <- tibble::tibble(row_data)
#   updated_df <- dplyr::bind_rows(df_csv, new_row)
#   readr::write_csv(updated_df, metrics_file)
# }

# metrics_file <- initialize_metrics_file(metrics_file)
# new_row <- list(
#   TRIAL_NAME = trial_name,
#   LAYER = blk,
#   nv = nv,
#   AF = AF,
#   AUROC = auc_value,
#   model_name = model_name
# )
# append_metrics_row(metrics_file, new_row)
# cat(sprintf("Metrics appended to %s.\n\n", metrics_file))