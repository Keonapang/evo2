if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("earth", quietly = TRUE)) install.packages("earth")

library(data.table)
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(earth))
suppressMessages(library(doParallel))
suppressMessages(library(parallel))
library(openxlsx)
library(ggpubr)


args <- commandArgs(trailingOnly = TRUE)
AF <- as.numeric(args[1]) # 0.00003
d <- as.numeric(args[2]) # 2 or 3
np <- as.numeric(args[3]) # 11, 13, 14
nv <- as.numeric(args[4]) # 0.1
trial_name <- args[5] # "A", "B", "C"
anno <- args[6] # "yes", "no"
cores <- args[7] # "yes", "no"

upper_MAF <- 0.01
p <-2
t <- 0.00001
# Directory 
server <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk"
DIR_EVO2 <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June25_embedding"
ROOTDIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/MARS"
OUTDIR <- paste0(ROOTDIR,"/", trial_name)

if (!dir.exists(ROOTDIR)) {dir.create(ROOTDIR, recursive = TRUE)}
if (!dir.exists(OUTDIR)) {dir.create(OUTDIR, recursive = TRUE)}

# Model name
upper_MAF <- format(as.numeric(upper_MAF), scientific = TRUE, digits = 4)
AF <- format(as.numeric(AF), scientific = TRUE, digits = 4)
model_name <- paste0(trial_name, "_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_anno", anno) 
cat("Model name: ", model_name, "\n\n")


RData_file <- paste0(OUTDIR, "/", model_name, "_scores.RData")
scores_file <- paste0(OUTDIR, "/", model_name, "_scores.txt")


# Train on one gene, test on another gene
if (trial_name == "A") {
    train_gene <- "BRCA1 and LDLR"
} else if (trial_name == "B") {
    train_gene <- "BRCA1"
} else if (trial_name == "C") {
    train_gene <- "LDLR"
} else if (trial_name == "D") {
    train_gene <- "BRCA1"
    test_gene <- "BRCA1"
} else if (trial_name == "E") {
    train_gene <- "LDLR"
    test_gene <- "LDLR"
} else {
    stop("Invalid trial_name. Must be one of: A, B, C, D, E.")
}


# Train on one gene, test on another gene
if (trial_name == "B") {
    test_gene <- "LDLR"
} else if (trial_name == "C") {
    test_gene <- "BRCA1"
}

# variables

font <- "Sans"
psize <- 15
teal <- "#0B7C9A"
red <- "#DE433C"
green <- "#76B900"
font <- "Open Sans"
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
    
        cat("======================Trial", trial_name, " model ", i, "=======================\n")
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

cat("Saved:", plot_file, "\n\n")

####################################################################################
# PLOT: Plot distribution of scores
####################################################################################

final_df <- fread(scores_file)
colnames(final_df)[ncol(final_df)] <- "y_df1"

p <- ggplot() +
  geom_density(
    data = final_df,
    aes(x = y_df1, y = after_stat(scaled)), 
    alpha = 0.3, size = 1, adjust = 1.5, 
    color = "#C3535E", fill = "#C3535E", 
    show.legend = FALSE  # Remove legend
  ) +

  labs(x = "MARS scores", y = "Proportion") +  # Removed the fill label
  
  # Adjust x-axis limits and breaks to fit y_df1's range
  scale_x_continuous(
    limits = c(-0.5, 1),  # Adjust based on the range of y_df1
    breaks = seq(-0.5, 1, by = 0.25), 
    expand = c(0, 0)
  ) +
  theme_minimal() +
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

    legend.position = "none"  # Remove the legend completely
  )
plot_file <- paste0(OUTDIR, "/", model_name, "_scores.png")
ggsave(plot_file, plot = p, width = 17, height = 6.5, dpi = 300)
cat("File:", plot_file, "\n")

p <- ggplot(final_df, aes(x = y_df1)) +
  geom_histogram(
    binwidth = 0.05,  # Adjust bin width as needed
    color = "#C3535E",  # Stick color
    fill = "#C3535E",   # Stick fill
    alpha = 0.8         # Transparency
  ) +
  labs(x = "MARS scores", y = "Count") +  # Update y-axis label
  
  # Adjust x-axis limits and breaks to fit y_df1's range
  scale_x_continuous(
    limits = c(-0.5, 1.25),  # Adjust based on the range of y_df1
    breaks = seq(-0.5, 1.25, by = 0.25), 
    expand = c(0, 0)
  ) +
  theme_minimal() +
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
    plot.title = element_text(size = 20, hjust = 0, face = "bold", family = font)  # Center-align the plot title
  )

# Save the plot
plot_file <- paste0(OUTDIR, "/", model_name, "_scores2.png")
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
    title = paste0("ROC Curve for ", test_gene, "  (model: ", trial_name, ")"),
    x = "False Positive Rate", #  (1 - Specificity)
    y = "True Positive Rate" # (Sensitivity)
  ) +
  theme_minimal() +
  plot_theme
output_file <- paste0(OUTDIR,"/", model_name, "_auc.png")
ggsave(output_file, auc_plot, width = 6, height = 6, dpi = 350)
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
cat("\n")

# Filter the dataframe to include only "Pathogenic" and "Benign" classes
# filtered_df <- final_df[final_df$clinvar_clnsig %in% c("Pathogenic", "Benign"), ]
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
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +  # Dynamically set ticks, approx. 8 intervals
  labs(
    title = paste0("Distribution of pathogenic and benign RVs in ", test_gene),
    x = paste0("MARS scores (trial: ", trial_name, ")"), 
    y = "ClinVar classification"
  ) +
  theme_minimal() +  # Minimal theme
  plot_theme

output_file <- paste0(OUTDIR,"/", model_name, "_dot.png")
ggsave(output_file, plot, width = 10, height = 4.5, dpi = 300)
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

font <- "Open Sans"
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





###################################################################################################
# VALIDATION: how does the embedding-trained MARS model scores perform compared to 
#             reference evo2 scores?
# Plot scatter plot of these two scores on x and y-axis!
###################################################################################################
d=3
np=2000
nv=0.1 # 0.1 or 0
trial_name="B" # "A", "B", "C", "D", "E"
anno="yes" # "yes" or "no"

model_name <- paste0(trial_name, "_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF0e+00_anno", anno) 
cat("Model name: ", model_name, "\n\n")

# variables
if (trial_name == "A") {
    test_gene <- "BRCA1 and LDLR"
} else if (trial_name == "B") {
    test_gene <- "LDLR"
} else if (trial_name == "C") {
    test_gene <- "BRCA1"
} else if (trial_name == "D") {
    test_gene <- "BRCA1"
} else if (trial_name == "E") {
    test_gene <- "LDLR"
} else {
    stop("Invalid trial_name. Must be one of: A, B, C, D, E.")
}

# Directories
DIR_MARS <- paste0("/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/MARS/", trial_name)
DIR_EVO <- paste0("/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June23")

# output files
outfile <- paste0(DIR_MARS, "/", model_name, "_embed_vs_evo2.png") # clinvar_clnsig	, yhat

# input files
file_mars <- paste0(DIR_MARS, "/", model_name, "_scores.txt") # clinvar_clnsig	, yhat
file_evo2a <- paste0(DIR_EVO, "/RovHer_BRCA1_win8192_seq812_all.xlsx")
file_evo2b <- paste0(DIR_EVO, "/RovHer_LDLR_win4096_seq470_all.xlsx")

# load
data_mars <- fread(file_mars) # embedding MARS="yhat"
colnames(data_mars)[ncol(data_mars)] <- "embedding_yhat"
data_mars <- data_mars %>% select(PLINK_SNP_NAME, Gene.refGene, clinvar_clnsig, embedding_yhat)

data_evo2 <- rbind(read.xlsx(file_evo2a), read.xlsx(file_evo2b))
data_evo2 <- data_evo2 %>% select(PLINK_SNP_NAME, Gene.refGene, GENEBASS_AF, class, evo2_delta_score)

dim(data_mars) # 681   5 #     PLINK_SNP_NAME clinvar_clnsig      yhat
dim(data_evo2) # 1282  14; "evo2_delta_score"

# merge data_mars with data_evo2
data <- merge(data_mars, data_evo2, by = "PLINK_SNP_NAME", all.x = TRUE)
dim(data)
head(data,2)
colnames(data)

# Calculate Spearman's rank correlation
spearman_corr <- cor(data$embedding_yhat, data$evo2_delta_score, method = "spearman")

x_label <- max(data$embedding_yhat, na.rm = TRUE)
y_label <- max(data$evo2_delta_score, na.rm = TRUE)

font <- "Sans"
p <- ggplot(data, aes(x = embedding_yhat, y = evo2_delta_score)) +
  geom_point(color = teal, size = 2, alpha = 0.65) + 
  labs(
    title = paste0(test_gene, " rare variants"),
    x = "Embedding-trained MARS scores",
    y = "Evo2 7B scores"
  ) +
annotate(
    "label", 
    x = x_label, y = y_label, 
    label = paste0("Spearman's ρ = ", round(spearman_corr, 4)),
    hjust = 1.1, vjust = 1.1,  # Adjust alignment
    size = 5, 
    family = font, 
    label.size = 0.5,  # Border thickness
    fill = "white",    # Background color
    color = "black"    # Text color
  ) +
  theme_minimal() +
  theme(
  panel.background = element_rect(fill = "white", color = NA),  # White background
  plot.background = element_rect(fill = "white", color = NA),   # White plot area
  plot.margin = margin(t = 15, r = 15, b = 20, l = 20),         # Adjust plot margins
  axis.title.x = element_text(margin = margin(t = 15, b = 0.8), family = font),  # Padding for x-axis
  axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0.8), family = font),  # Padding for y-axis
  axis.text = element_text(size = 14, family = font),  # Increase axis text size
  axis.title = element_text(size = 14, face = "bold", family = font),  # Increase axis title size
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = font),  # Centered, bold title
)
ggsave(outfile, p, width = 18, height = 6, dpi = 300)
cat("Plot:", outfile, "\n")
