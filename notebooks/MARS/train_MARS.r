# Train MARS using two genes (BRCA1 and LDLR) 
# Adding Evo2 7B model embeddings to train set + 75 functional annotations from RovHer study 
# June 26-30, 2025

# Trial name scenarios:
# A      Train on variants of both genes - 5 fold CV
# B      Train on variants of BRCA1, test on LDLR 
# C      Train on variants of LDLR, test on BRCA1 
# D      Train on BRCA1, test on held-out set of BRCA1  - 5 fold CV
# E      Train on LDLR, test on held-out set of LDLR  - 5 fold CV

# trial_name="B C" # 
# d=3
# np=2000 # Reducing nprune reduces exhaustive search time; maximum number of permissible terms in the final pruned model.
# nv=0.1
# AF=0.00003
# cores=4
# anno_cols="yes"
# embedding_cols="yes"
# # nk (maximum number of  forward pass terms)

# for t in $trial_name; do
#   Rscript /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/MARS/train_MARS.r $AF $d $np $nv $t $anno_cols $embedding_cols $cores
# done

if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("earth", quietly = TRUE)) install.packages("earth")
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
embedding_col <- as.character(args[7]) # | "delta" (delta embeddings), "no" (no embeddings), "refvar" (ref + var embeddings)
blk <- as.character(args[8])
cores <- as.numeric(args[9]) # "yes", "no"
date <- format(Sys.Date(), "%Y%m%d")
cat("date:", date, "\n")
# trial_name <- "A"
# anno <- "yes"
# d <- 3
# np <- 2000
# nv <- 0.1
# AF <- 0
# embedding_col <- "refvar"
# blk <- "28"
# cores <- 24
upper_MAF <- 0.01
p <-2
t <- 0.00001

registerDoParallel(makeCluster(cores))

# Output model name
upper_MAF <- format(as.numeric(upper_MAF), scientific = TRUE, digits = 4)
AF <- format(as.numeric(AF), scientific = TRUE, digits = 4)

if (embedding_col == "delta" || embedding_col == "refvar") {
    model_name <- paste0(trial_name, "_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_anno", anno, "_embed", embedding_col, "_blk", blk) 
} else if (embedding_col == "no") {
    model_name <- paste0(trial_name, "_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_anno", anno, "_embed", embedding_col) 
} else {
    stop("Invalid 'embedding_col'. Must be either: delta, refvar, no")
}
cat("Model name: ", model_name, "\n\n")

# Input directory 
server <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk"
DIR_EVO2 <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/July1_embedding"

# Output directory
ROOTDIR <- paste0("/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/MARS/",date, "_", embedding_col) # <------ MODIFY!
OUTDIR1 <- paste0(ROOTDIR,"/", trial_name)
OUTDIR <- paste0(OUTDIR1,"/", model_name)

if (!dir.exists(ROOTDIR)) {dir.create(ROOTDIR, recursive = TRUE)}
if (!dir.exists(OUTDIR1)) {dir.create(OUTDIR1, recursive = TRUE)}
if (!dir.exists(OUTDIR)) {dir.create(OUTDIR, recursive = TRUE)}

# Output files
RData_file <- paste0(OUTDIR, "/", model_name, "_scores.RData")
scores_file <- paste0(OUTDIR, "/", model_name, "_scores.txt")
metrics_file <- paste0(OUTDIR, "/layer", blk, "_metrics.txt")

####################################################################################
# INPUT FILES
####################################################################################
# 1. Variant and gene level anno (old training set)
train_path <- paste0(server, "/Training/height_FDR/genebass_AF001_chrall_vars_gene_final.txt")  # 4926639 x  80 

# Evo2 scores
DIR_EVO <- paste0("/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June23")
file_evo2a <- paste0(DIR_EVO, "/RovHer_BRCA1_win8192_seq812_all.xlsx")
file_evo2b <- paste0(DIR_EVO, "/RovHer_LDLR_win4096_seq470_all.xlsx")

# 2. csv file with embeddings
csv1 <- paste0(DIR_EVO2,"/RovHer_BRCA1_blocks.", blk,".mlp.l3_delta.csv")
csv2 <- paste0(DIR_EVO2,"/RovHer_LDLR_blocks.", blk,".mlp.l3_delta.csv")
csv1_df <- fread(csv1, header = TRUE)
csv2_df <- fread(csv2, header = TRUE)
colnames(csv1_df) <- ifelse(grepl("^e", colnames(csv1_df)), 
                            paste0(colnames(csv1_df), "_delta"), 
                            colnames(csv1_df))

colnames(csv2_df) <- ifelse(grepl("^e", colnames(csv2_df)), 
                            paste0(colnames(csv2_df), "_delta"), 
                            colnames(csv2_df))

csv1ref <- paste0(DIR_EVO2,"/RovHer_BRCA1_blocks.",blk,".mlp.l3_ref.csv")
csv2ref <- paste0(DIR_EVO2,"/RovHer_LDLR_blocks.",blk,".mlp.l3_ref.csv")
csv1ref_df <- fread(csv1ref, header = TRUE)
csv2ref_df <- fread(csv2ref, header = TRUE)
colnames(csv1ref_df) <-ifelse(grepl("^e", colnames(csv1ref_df)), 
                            paste0(colnames(csv1ref_df), "_ref"), 
                            colnames(csv1ref_df))

colnames(csv2ref_df) <-ifelse(grepl("^e", colnames(csv2ref_df)), 
                            paste0(colnames(csv2ref_df), "_ref"), 
                            colnames(csv2ref_df))

csv1var <- paste0(DIR_EVO2,"/RovHer_BRCA1_blocks.", blk,".mlp.l3_var.csv")
csv2var <- paste0(DIR_EVO2,"/RovHer_LDLR_blocks.", blk,".mlp.l3_var.csv")
csv1var_df <- fread(csv1var, header = TRUE)
csv2var_df <- fread(csv2var, header = TRUE)
colnames(csv1var_df) <-ifelse(grepl("^e", colnames(csv1var_df)), 
                            paste0(colnames(csv1var_df), "_var"), 
                            colnames(csv1var_df))
colnames(csv2var_df) <-ifelse(grepl("^e", colnames(csv2var_df)), 
                            paste0(colnames(csv2var_df), "_var"), 
                            colnames(csv2var_df))


# if column "layer" exists, remove it
if ("layer" %in% colnames(csv1_df)) {
  csv1_df <- csv1_df[, !("layer"), with = FALSE]
}
if ("layer" %in% colnames(csv2_df)) {
  csv2_df <- csv2_df[, !("layer"), with = FALSE]
}
if ("region" %in% colnames(csv1_df)) {
  setnames(csv1_df, old = "region", new = "input_file")
}
if ("region" %in% colnames(csv2_df)) {
  setnames(csv2_df, old = "region", new = "input_file")
}
if (!all(colnames(csv1_df) == colnames(csv2_df))) {
    stop("CSV files have different columns. Please check the input files.")
}

# 3. predict path
predict_path <- train_path

# Dependent variable - ClinGen gold standard annotations
ACMG_file1 <- paste0(server,"/RARity_monogenic_benchmark/BRCAexchange/BRCA1_clinvar_cleaned.txt") 
ACMG_file2 <- paste0(server, "/RARity_monogenic_benchmark/LOVD_LDLR/LDLR_clinvar_curated.txt") # all British heart foundation-classified variants published on LOVD

####################################################################################
# LOAD FILES
####################################################################################

# Load embedding columns
if (embedding_col == "delta" || embedding_col == "no") {
  if (trial_name == "A" || trial_name == "B" || trial_name == "C") {
      train_gene <- "BRCA1 and LDLR"
      test_gene <- "BRCA1 and LDLR"
      embed_cols <- rbind(csv1_df, csv2_df)
      # embed_cols <- rbind(fread(csv1, header = TRUE), fread(csv2, header = TRUE))
  } else if (trial_name == "D") {
      cat("Trial D: Training on held-out BRCA1 variants, testing on held-out BRCA1 variants\n")
      embed_cols <- csv1_df
      train_gene <- "BRCA1"
      test_gene <- "BRCA1"
  } else if (trial_name == "E") {
      cat("Trial E: Training on held-out LDLR variants, testing on held-out LDLR variants\n")
      embed_cols <- csv2_df
      train_gene <- "LDLR"
      test_gene <- "LDLR"
  } else {
      stop("Invalid trial_name. Must be one of: A, B, C, D, E.")
  }
}

if (embedding_col == "refvar") {
  if (trial_name == "A" || trial_name == "B" || trial_name == "C") {
      train_gene <- "BRCA1 and LDLR"
      test_gene <- "BRCA1 and LDLR"
      ref_embed <- rbind(csv1ref_df, csv2ref_df) # reference embeddings from both genes 
      var_embed <- rbind(csv1var_df, csv2var_df)  # reference embeddings from both genes 
      embed_cols <- cbind(ref_embed, var_embed) # combine reference and variant embeddings
  } else if (trial_name == "D") {
      cat("Trial D: Training on held-out BRCA1 variants, testing on held-out BRCA1 variants\n")
      embed_cols <- cbind(csv1ref_df, csv1var_df)
      train_gene <- "BRCA1"
      test_gene <- "BRCA1"
  } else if (trial_name == "E") {
      cat("Trial E: Training on held-out LDLR variants, testing on held-out LDLR variants\n")
      embed_cols <- cbind(csv2ref_df, csv2var_df)
      train_gene <- "LDLR"
      test_gene <- "LDLR"
  } else {
      stop("Invalid trial_name. Must be one of: A, B, C, D, E.")
  }
}

# check if there are duplicate columns for  [PLINK_SNP_NAME, input_file, RovHer_score, layer].
if (any(duplicated(colnames(embed_cols)))) {
    embed_cols <- embed_cols[, !duplicated(colnames(embed_cols)), with = FALSE]
    cat("Removed duplicate columns from embed_cols. New dimensions:", dim(embed_cols), "\n\n")
} else {
    cat("No duplicate columns found in embed_cols.\n\n")
}

cat("embed_cols for", embedding_col, ":",dim(embed_cols), "\n\n")
head(embed_cols[,1:5],2)
cat("\n")

#  dependent variable - scores
ACMG_col1 <- fread(ACMG_file1, header = TRUE, select = c("PLINK_SNP_NAME", "ACMG_final")) 
ACMG_col1 <- rename(ACMG_col1, clinvar_clnsig = ACMG_final)
ACMG_col2 <- fread(ACMG_file2, header = TRUE, select = c("PLINK_SNP_NAME", "clinvar_clnsig"))
ACMG_cols <- rbind(ACMG_col1, ACMG_col2) 
cat("ACMG_cols merged:",dim(ACMG_cols), "\n")

# Load Evo2 delta likelihood scores
data_evo2 <- rbind(read.xlsx(file_evo2a), read.xlsx(file_evo2b))
data_evo2 <- data_evo2 %>% select(PLINK_SNP_NAME,evo2_delta_score)

# 1. Build Training set, then add y-dependent variable (ClinVar classifications for each RV)
anno_cols <- fread(train_path, header = TRUE)
anno_cols <- merge(anno_cols, data_evo2, by = "PLINK_SNP_NAME", all.x = TRUE)
cat("anno cols:\n\n")
colnames(anno_cols)
cat("\n")

if (anno == "yes" && embedding_col == "delta" || anno == "yes" && embedding_col == "refvar") {
    ("Loaded and merged annotation and embedding columns....\n")
    df <- merge(embed_cols, anno_cols, by = "PLINK_SNP_NAME", all.x = TRUE)
} else if (anno == "no" && embedding_col == "delta" || anno == "no" && embedding_col == "refvar") {
    ("Loaded embedding columns only....\n")
    anno_cols <- anno_cols[, .(PLINK_SNP_NAME, GENEBASS_AF, `Gene.refGene`)] # subset version
    df <- merge(embed_cols, anno_cols, by = "PLINK_SNP_NAME", all.x = TRUE)
} else if (anno == "yes" && embedding_col == "no") {
    ("Loaded annotation columns only....\n")
    df <- anno_cols
} else {
    stop("Invalid combination of anno and embedding_col")
}

df <- merge(df, ACMG_cols, by = "PLINK_SNP_NAME", all.x = TRUE)
cat("df:\n\n")
head(df[,1:5],2)
cat("\n")

# 2. Recode the "clinvar_clnsig" column
df <- df %>% select(PLINK_SNP_NAME, clinvar_clnsig, everything()) # 1282 4175
df$clinvar_clnsig[is.na(df$clinvar_clnsig)] <- "NA" # 1286 x 4171
# head(df[,1:5])
# tail(df[,1:5])

# 3. Filter out rows where clinvar_clnsig is "", "NA", or "CCP"

df <- df[!df$clinvar_clnsig %in% c("", "NA", "CCP"), ]
# print(table(df$clinvar_clnsig, useNA = "ifany"))

df$clinvar_clnsig <- recode(df$clinvar_clnsig,
  "P" = 1,
  "B" = 0,
  "LB" = 0.25,
  "LP" = 0.75,
  "LP,P" = 0.75,  # LP and P are considered pathogenic
  "B/LB" = 0.25,   # LP and B are considered likely
  .default = 0.5  # Assign 0.5 to all other values
)
cat("Recoded categories:\n")
print(table(df$clinvar_clnsig, useNA = "ifany"))
cat("\n")
# 4. Get prediction set
df_pred <- df #  1282 4174
if ("layer" %in% colnames(df_pred)) {
  df_pred <- df_pred[, !("layer"), with = FALSE]
}
####################################################################################
# FILTER TRAINING SET
####################################################################################

# Apply MAF threshold (0.000003 < MAF < 0.01)
# count_rows <- df %>% filter(GENEBASS_AF < as.numeric(AF)) %>% nrow()
# df1 <- subset(df, GENEBASS_AF >= as.numeric(AF)) # Lower bound MAF  
# cat(count_rows, "RVs with MAF <", AF, "removed: ", dim(df1), "\n")  # (4844229 - 2682377) = 2114607 x 76 
# df1 <- subset(df1, GENEBASS_AF < upper_MAF) # Upper-bound MAF threshold
# dim(df1) 

# Train on one gene, test on another gene
if  (trial_name == "B") {
    df <- df[df$Gene.refGene == "BRCA1", ]
    train_gene <- "BRCA1"
    test_gene <- "LDLR"
} else if (trial_name == "C") {
    df <- df[df$Gene.refGene == "LDLR", ]
    train_gene <- "LDLR"
    test_gene <- "BRCA1"
}

# Extract and store a column of just PLINK_SNP_NAME
extra_cols <- data.frame(PLINK_SNP_NAME = df$PLINK_SNP_NAME,
                            Gene.refGene = df$`Gene.refGene`,
                            GENEBASS_AF = df$GENEBASS_AF,
                            clinvar_clnsig = df$`clinvar_clnsig`) # 1286 x 1

# Remove columns
# columns_to_remove <- c("PLINK_SNP_NAME", "Gene.refGene", "input_file", "height_FDR", "RovHer_score")
columns_to_remove <- c("PLINK_SNP_NAME", "Gene.refGene", "input_file", "height_FDR", "RovHer_score",
                       "VEST4_score","ClinPred_score","gMVP_score","fathmm-XF_coding_score",
                       "REVEL_score","M-CAP_score","MetaLR_score","MetaRNN_score","BayesDel_noAF_score","BayesDel_addAF_rankscore") 
                # fathmm-XF_coding_score" trained on HGMD and 1000G

df2 <- df %>% select(-any_of(columns_to_remove))
cat("df2:",dim(df2)) # [1] 1286 4171
if ("layer" %in% colnames(df2)) {
  df2 <- df2[, !("layer"), with = FALSE]
}
# Get x and y matrices
X_df1 <- dplyr::select(df2, -clinvar_clnsig) # or "evo2_delta_score",'height_FDR'
y_df1 = df2[["clinvar_clnsig"]] # or "evo2_delta_score"

cat("\n\nTraining set:\n\n")
head(X_df1[,1:4],2)
cat("\nX_df1:", dim(X_df1), "\n\n") # 1286 4170


####################################################################################
# FILTER PREDICTION SET
####################################################################################

# Train on one gene, test on another gene
if (trial_name == "B") {
    df_pred <- df_pred[df_pred$Gene.refGene == "LDLR", ]
} else if (trial_name == "C") {
    df_pred <- df_pred[df_pred$Gene.refGene == "BRCA1", ]
}

if (trial_name == "B" || trial_name == "C") {
    extra_cols <- data.frame(PLINK_SNP_NAME = df_pred$PLINK_SNP_NAME,
                            Gene.refGene = df_pred$`Gene.refGene`,
                            GENEBASS_AF = df_pred$GENEBASS_AF,
                            clinvar_clnsig = df_pred$`clinvar_clnsig`) # 1286 x 1
}

# Filter out columns
X_df_pred1 <- df_pred %>% select(-any_of(columns_to_remove))
X_df_pred1 <- dplyr::select(X_df_pred1, -clinvar_clnsig) # or "evo2_delta_score",'height_FDR'
cat("Test set:\n\n")
head(X_df_pred1[,1:4],2)
cat("X_df_pred1:", dim(X_df_pred1), "\n\n") # 1286 4170

##################################################################################### 
# No cross validation
####################################################################################
start_time <- Sys.time()

if (trial_name == "B" || trial_name == "C") {

cat("Training trial", trial_name, "on ", train_gene, " (n =", nrow(X_df1), " variants)...\n\n")
mars_model <- earth(x = X_df1, y = y_df1, degree = d, penalty = p, nprune = np, thresh = t,  newvar.penalty = nv, trace = 1)

# Predict on LDLR data
cat("\nPredicting on ", test_gene, "...\n\n")
predictions <- predict(mars_model, newdata = X_df_pred1)

# OUTPUT RESULTS
cat("Saving predictions...\n\n")
yhat_dataframe <- data.frame(yhat = as.numeric(predictions))
final_df <- cbind(extra_cols, yhat_dataframe) # 3128762 x 4
head(final_df,2)
cat("\n")
fwrite(final_df, file = scores_file, sep = "\t", quote = FALSE)
save(mars_model, final_df, d, p, t, np, nv, train_path, predict_path, file = RData_file)
}

####################################################################################
# 5-fold cross validation
####################################################################################

if (trial_name == "A" || trial_name == "D" || trial_name == "E") {

set.seed(42) # Set seed for reproducibility
k <- 5

# Assign random folds
df2$fold <- sample(rep(1:k, length.out = nrow(df2)))
X_df_pred1$fold <- sample(rep(1:k, length.out = nrow(X_df_pred1)))

# Initialize storage for predictions and models
all_predictions <- numeric(nrow(X_df_pred1))
cv_mars_list <- list() # To store models for each fold

# 5-fold cross-validation loop
cat("Training MARS for trial", trial_name, " with annocols=", anno, " (deg:", d, " np:", np," penalty:", p, ")\n")
for (i in 1:k) {
    cat("=============== [ CV Fold", i, "of", k, "] =================\n\n")

    # Split data into training and testing sets
    train_indices <- which(df2$fold != i) # Rows not in the current fold
    test_indices <- which(df2$fold == i) # Rows in the current fold

    df2_train <- df2[train_indices, ]
    df2_test <- df2[test_indices, ]

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
head(final_df,2)

save(final_df, cv_mars_list, d, p, t, np, nv, train_path, predict_path,
    file = RData_file)
fwrite(final_df, file=scores_file, sep = "\t", quote = FALSE)
}

# Calculate duration
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")
cat("------- Total Duration: ", round(duration, 2), " mins -------\n")

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
    
        cat("====================== Trial", trial_name, " model ", i, "=======================\n")
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
        anno_cols <- fread(train_path, header = TRUE)
        gene_cols <- colnames(anno_cols)[(48 + 1):ncol(anno_cols)]
        variant_cols <- colnames(anno_cols)[1:48]

        plot_data$Category <- ifelse(
            grepl("^e_", plot_data$Predictor),  # Predictors starting with "e_"
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
    anno_cols <- fread(train_path, header = TRUE)
    gene_cols <- colnames(anno_cols)[(48 + 1):ncol(anno_cols)]
    variant_cols <- colnames(anno_cols)[1:48]

    plot_data$Category <- ifelse(
        grepl("^e[0-9]+$", plot_data$Predictor),  # Predictors starting with "e" and ending with numbers
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


#################################################################################################
# Plot ROC curve and calculate AUC
#################################################################################################
suppressMessages(library(pROC))
final_df <- fread(scores_file)
# rename last column to y_df1
colnames(final_df)[ncol(final_df)] <- "y_df1"

# Remove rows where clinvar_clnsig is 0.5
final_df <- final_df[final_df$clinvar_clnsig != 0.5, ]

# Step 1: Prepare ground truth as binary (P/LP = 1, B/LB = 0)
final_df$binary_truth <- ifelse(final_df$clinvar_clnsig %in% c(1, 0.75), 1, 0)

# Verify the results
table(final_df$binary_truth)  # Check the counts of the binary classes
head(final_df)  # Check the modified dataframe

# Step 2: Calculate AUROC
roc_obj <- roc(final_df$binary_truth, final_df$y_df1)
auc_value <- auc(roc_obj)
cat("AUC:", auc_value, "\n")

# Step 3: Extract data for plotting the ROC curve
roc_data <- data.frame(
  TPR = rev(roc_obj$sensitivities), # True Positive Rate (Sensitivity)
  FPR = rev(1 - roc_obj$specificities) # False Positive Rate (1 - Specificity)
)

auc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "black", size = 1.2) +  # ROC curve
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = red) +  # Random chance line
  annotate(
        "label",  # Use "label" instead of "text" to add a text box with a border
        x = 0.75, y = 0.2, 
        label = paste("AUC =", round(auc_value, 3)), 
        size = 5, 
        color = "black", 
        fill = "white",  # Background color of the label
        label.size = 0.2  # Border thickness (thin black border)
    ) +
    labs(
    title = paste0("ROC Curve for ", test_gene, "  (MARS model: ", trial_name, ")"),
    x = "False Positive Rate", #  (1 - Specificity)
    y = "True Positive Rate" # (Sensitivity)
  ) +
  theme_minimal() +
  plot_theme
output_file <- paste0(OUTDIR,"/", model_name, "_auc.png")
ggsave(output_file, auc_plot, width = 10, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")


####################################################################################
# PLOT: Distribution of Evo2 scores of two ClinVar classes (P and B)
####################################################################################
final_df <- fread(scores_file)
# rename last column to y_df1
colnames(final_df)[ncol(final_df)] <- "y_df1"

# Recode clinvar_clnsig into descriptive categories
final_df$clinvar_clnsig <- factor(
  final_df$clinvar_clnsig,
  levels = c(0, 0.25, 0.5, 0.75, 1),  # Original numerical values
  labels = c("Benign", "Likely Benign", "VUS", "Likely Pathogenic", "Pathogenic")  # New labels
)

table(final_df$clinvar_clnsig)
head(final_df) #  116   5

# Filter the dataframe to include only "Pathogenic" and "Benign" classes
filtered_df <- final_df[final_df$clinvar_clnsig %in% c("Pathogenic", "Likely Pathogenic", "Likely Benign", "Benign"), ]

red <- "#DE433C"
green <- "#76B900"

plot <- ggplot(filtered_df, aes(x = y_df1, y = clinvar_clnsig, color = clinvar_clnsig)) +
  geom_jitter(height = 0.2, width = 0, size = 2, alpha = 0.8) +  # Jittered points
  geom_boxplot(aes(group = clinvar_clnsig), width = 0.4, alpha = 0.2, outlier.shape = NA, color = "black") +  # Boxplot
  scale_color_manual(
    values = c(
      "Pathogenic" = red,
      "Likely Pathogenic" = red,
      "Likely Benign" = green,
      "Benign" = green
    )
  ) +
  scale_x_continuous(breaks = seq(-1, 1.5, by = 0.25)) +  # Dynamically set x-axis ticks at every 0.25
  labs(
    title = paste0("Distribution of pathogenic and benign RVs in ", test_gene),
    x = paste0("MARS scores (trial: ", trial_name, ")"), 
    y = "ClinVar classification"
  ) +
  theme_minimal() +  # Minimal theme
  plot_theme


output_file <- paste0(OUTDIR,"/", model_name, "_dot.png")
ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")


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


###############################################################################################
# Write metrics file in excel to save AUROC
###############################################################################################

initialize_metrics_file <- function(file_path) {
  metrics_csv <- file_path
    columns <- c(
    "MODEL_SIZE", "input_file", "subset_method", "USE_RANDOM_SEED", 
    "SEQ_LENGTH", "WINDOW_SIZE", "time_score_ref", "time_score_var", "AUROC"
  )
    if (!file.exists(metrics_csv)) {
    df <- tibble::tibble(!!!setNames(vector("list", length(columns)), columns))
    readr::write_csv(df, metrics_csv)
  }
  return(metrics_csv)
}

append_metrics_row <- function(metrics_file, row_data) {
  df_csv <- readr::read_csv(metrics_file, show_col_types = FALSE)
  new_row <- tibble::tibble(row_data)
  updated_df <- dplyr::bind_rows(df_csv, new_row)
  readr::write_csv(updated_df, metrics_file)
}

metrics_file <- initialize_metrics_file(metrics_file)
new_row <- list(
  TRIAL_NAME = trial_name,
  LAYER = blk,
  nv = nv,
  AF = AF,
  AUROC = auc_value,
  model_name = model_name
)
append_metrics_row(metrics_file, new_row)
cat(sprintf("Metrics appended to %s.\n\n", metrics_file))