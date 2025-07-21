#!/usr/bin/env Rscript
# NEW: saves all 5 cv_mars models!!!
# Date: Oct 8, 2024
#####################################################################################################################################################
# Descrption: Training one 'universal' MARS model with a single-trait's FDR
#              5-fold CV (per chr region). Then saving FDR predictions on all ~ 3 million RVs from the UKB 200K WES. 
#             This subset of RV will be the same across all cont' traits. 
# Model usage: Can also be used for multi-trait MARS model (for height_FDR, LDL_FDR, BMI_FDR etc...)
# Note: these scores were not performing well on the Mendelian analysis, these scores work well for the h2 curve only as they enrich for functional RV but not necessarily the ultra-rare highly penetrance Mendelian causal variants
# ####################################################################################################################################################
# trait="height" # LDL_direct, BMI, height, Cystatin_C, Alkaline_phosphatase, alanine_aminotransferase, Creatinine, Aspartate_aminotransferase
# train_path="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Training/height_FDR/genebass_AF001_chrall_vars_gene_final.txt"  # 4926639 x  80 
# d=3
# np=2000            # np=14, np=16, np=18, np=NULL
# nv=0.1
# AF=0.00003         # lower MAF (AF=0.000003,0.000004, 0.000006, 0.00001)
# upper_MAF=0.01     # upper MAF
# model="her"        # "her"" or "men""
# /usr/bin/time -v Rscript "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/20231101_annotation/SCRIPT_12_GENEBASS_RV_train_RV_pred_Oct10.r" $trait $AF $train_path $d $np $nv $model $upper_MAF > "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Training/height_FDR/log/cv_MARS_d${d}_np${np}_nv${nv}_lowerAF${AF}_upperAF${upper_MAF}_Feb15_${model}.log" 2>&1 &

#####################################################################################################################################################
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
registerDoParallel(makeCluster(cores))

# AF <- 0 # Lower bound MAF threshold
# d <- 3
# np <- 2000
# nv <- 0
# trial_name <- "Chr17" # "A", "B", "C", "D", "E"
# anno <- "yes" # "yes", "no"
# embedding_col <- "no"
# blk <- "28" # "28", "29", "30", "31", "32", "33", "34", "35"
# cores <- 30 # Number of cores to use for parallel processing
# date <- 20250717

upper_MAF <- 0.01
p <-2
t <- 0.00001
nk <- 600
if (np == "NULL") {
  np <- NULL
} else {
  np <- as.numeric(np)
}

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
cat("\nModel name: ", model_name, "\n\n")

# Input directory 
server <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk"
dir <- paste0(server,"/Training/height_FDR/GENEBASS_RV_train_RV_pred_Oct10") # old folder (no cross validation): allvar_predict_saved_models_NO_GENESETS 
DIR_SCORES <- paste0(server,"/evo2/June30")

# Output directory
ROOTDIR <- paste0(server,"/evo2/MARS/",date, "_", embedding_col) # <------ MODIFY!
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
# Functional annotations
train_path <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Training/height_FDR/genebass_AF001_chrall_vars_gene_final.txt"  # 4926639 x  80 
predict_path <- train_path

# Evo2 scores
file_score <- paste0(DIR_SCORES, "/evo2_chr17_win4096.txt")

####################################################################################
# Load training matrix
####################################################################################
train_df <- as.data.frame(fread(train_path, header = TRUE))
cat("train_df: ", dim(train_df)[1], "RVs x", dim(df)[2], "cols\n\n") # 'df' of size 3 million x 97 cols

# Load Evo2 scores
evo2_scores <- as.data.frame(fread(file_score, header = TRUE, 
                            select = c("PLINK_SNP_NAME", "evo2_delta_score"),
                            showProgress = FALSE))
train_df <- merge(train_df, evo2_scores, by = "PLINK_SNP_NAME", all.x = TRUE) #  left join

####################################################################################
# Filter training matrix
####################################################################################

# Keep only certain chromosomes RVs
train_df <- train_df %>% filter(grepl("^17:", PLINK_SNP_NAME)) 
cat("train_df Chr 17-only: ", dim(train_df)[1], "RVs x", dim(df)[2], "cols\n\n") # 'df' of size 3 million x 97 cols
colnames(train_df)
cat("\n")

# Subset to get columns for final prediction file 
extra_cols <- train_df %>% select(PLINK_SNP_NAME, GENEBASS_AF, LOF_DISRUPTIVE, `Gene.refGene`)

# Subset training data 
df <- dplyr::select(train_df, -PLINK_SNP_NAME, -`Gene.refGene`) 

# Keep RVs with MAF above the lower bound threshold 
count_rows <- df %>% filter(GENEBASS_AF < as.numeric(AF)) %>% nrow()
cat("Number of GENEBASS RVs with MAF <", AF, "(lower thresh):", count_rows, "\n")
df1 <- subset(df, GENEBASS_AF >= as.numeric(AF)) # Lower bound (4844229 - 2682377) = 2114607 x 76 
df1 <- subset(df1, GENEBASS_AF < upper_MAF) # Upper-bound 

print(head(df1[, 1:4]),3)
cat("\n")

X_df1 = df1[, -1] 
y_df1 = df1[["height_FDR"]]
cat("Training preds:\n\n")
colnames(X_df1)
cat("\n")
cat("'X' matrix MAF >", AF, "removed:", dim(X_df1)[1], "x", dim(X_df1)[2], "\n")
cat("'y FDR:", length(y_df1), "x 1\n")  
cat("\n")
print(head(X_df1,2))
cat("\n")
cat("----------------------------------------------------------\n\n")

####################################################################################
# Load test matrix 
####################################################################################
# df <- as.data.frame(fread(predict_path, header = TRUE, showProgress = FALSE))
df_pred <- df

X_df_pred1 = df_pred[, -1] # 3mill x 92
cat("\n")
print(head(X_df_pred1,2))
cat("\n")
y = df_pred[["height_FDR"]]
cat("----------------------------------------------------------\n\n")

set.seed(42) # Set seed for reproducibility
k <- 5

# Assign random folds
X_df_pred1$fold <- sample(rep(1:k, length.out = nrow(X_df_pred1)))

# Initialize storage for predictions and models
all_predictions <- numeric(nrow(X_df_pred1))
cv_mars_list <- list() # To store models for each fold
start_time <- Sys.time()

# 5-fold cross-validation loop
cat("Training MARS for trial", trial_name, " with annocols=", anno, " (deg:", d, " np:", np," penalty:", p, ")\n")
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
        anno_cols <- train_df
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
    anno_cols <- train_df
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