
library(data.table)
library(ggplot2)
library(showtext)
library(sysfonts)
library(readxl)
library(openxlsx)
library(dplyr)
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# Intialize the plot variables
font <- "Open Sans"
plot_theme <- 
  theme_minimal(base_family = font) +
  theme(
  panel.background = element_rect(fill = "white", color = NA),  # White background
  plot.background = element_rect(fill = "white", color = NA),   # White plot area
  plot.margin = margin(t = 15, r = 15, b = 20, l = 20),         # Adjust plot margins
  axis.title.x = element_text(margin = margin(t = 15, b = 0.8)),  # Padding for x-axis
  axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0.8)),  # Padding for y-axis
  axis.text = element_text(size = 14),  # Increase axis text size
  axis.title = element_text(size = 16, face = "bold"),  # Increase axis title size
  plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Centered, bold title
  legend.position = "none"  # Remove the legend
)

font <- "Open Sans"
red <- "#DE433C"
nice_red <- "#CC6a66"
nvidia_green <- "#76B900"
teal <- "#0B7C9A"
lightblue <- "#9ec1da"
grey_blue <- "#668298"
yellow <- "#F2C46D"
yellow <- "#F6D874"
pale_orange <- "#F9B987"
color3 <-"#C89BBD"
color4 <- "#EA7F85"
color6 <-"#C7D59E"

###########################################################################
# Scatter plot: Input sequence length vs time for different window sizes
# gauge runtime for different parameters
###########################################################################

# Input directory
dir <- "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/results"
# Input file
file <- "metrics.xlsx"

# Output plot files
outfile1 <- "plot1.png"
outfile2 <- "plot2.png"
outfile3 <- "plot3.png"
outfile4 <- "plot_3D.png"

# Set dir
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE)}
setwd(dir)
dt <- read_excel(file)
# remove rows where the SEQ_LENGTH is 3890
dt <- dt[dt$SEQ_LENGTH != 3890, ]

# Calculate the time (sum of score_ref and score_var)
dt$Time <- dt$score_ref + dt$score_var

# Plot 1: 
ggplot(data = dt, aes(x = SEQ_LENGTH, y = Time, color = as.factor(WINDOW_SIZE), shape = as.factor(WINDOW_SIZE))) +
  geom_point(size = 3) +  # Add scatter plot points
  geom_line(aes(group = WINDOW_SIZE), size = 1) +  # Add lines grouped by WINDOW_SIZE
  scale_x_continuous(name = "Sequence Length", limits = c(10, 2000)) +  # Set x-axis limits
  scale_y_continuous(name = "Time (s)", limits = c(70, 1480)) +  # Set y-axis limits
  labs(color = "Window Size", shape = "Window Size") +  # Add legend titles
  plot_theme
ggsave(outfile1, width = 10, height = 6, dpi = 300)

# Plot 2
ggplot(data = dt, aes(x = SEQ_LENGTH, y = score_ref, color = as.factor(WINDOW_SIZE), shape = as.factor(WINDOW_SIZE))) +
  geom_point(size = 3) +  # Add scatter plot points
  geom_line(aes(group = WINDOW_SIZE), size = 1) +  # Add lines grouped by WINDOW_SIZE
  scale_x_continuous(name = "Sequence Length", limits = c(10, 2000)) +  # Set x-axis limits
  scale_y_continuous(name = "Time taken to score reference sequence (s)", limits = c(35, 1085)) +  # Set y-axis limits
  labs(color = "Window Size", shape = "Window Size") +  # Add legend titles
  plot_theme
ggsave(outfile2, width = 10, height = 6, dpi = 300)

# Plot 3
ggplot(data = dt, aes(x = SEQ_LENGTH, y = score_var, color = as.factor(WINDOW_SIZE), shape = as.factor(WINDOW_SIZE))) +
  geom_point(size = 3) +  # Add scatter plot points
  geom_line(aes(group = WINDOW_SIZE), size = 1) +  # Add lines grouped by WINDOW_SIZE
  scale_x_continuous(name = "Sequence Length", limits = c(10, 2000)) +  # Set x-axis limits
  scale_y_continuous(name = "Time taken to score variant sequence (s)", limits = c(35, 1085)) +  # Set y-axis limits
  labs(color = "Window Size", shape = "Window Size") +  # Add legend titles
  plot_theme
ggsave(outfile3, width = 10, height = 6, dpi = 300)


# June 23
###################################################################################################
# 1. Distribution of Evo2 scores of high/low RovHer scores for Chr 17, chr19 or a gene (e.g. BRCA1)
# OUTPUT: RovHer_chr19.txt_win4096_seq1000_dot.png and RovHer_chr17.txt_win4096_seq1000_dot.png
###################################################################################################
# Variables
seq = 100 # 50,84,100,250,500,1000,2000,5000,8000,10000
WINDOW_SIZE = 2048 # 1024,2048,4096,8192
chr="BRCA1"

# Input dir
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June23"
# Output Dir
OUTDIR <- paste0(DIR,"/plot")
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR, recursive = TRUE)
}

# Input file - a "bottom" and "top" pair for
bottom <- read_excel(paste0("RovHer_",chr,"_win", WINDOW_SIZE, "_seq",seq, "_bottom.xlsx"))
top <- read_excel(paste0("RovHer_",chr,"_win", WINDOW_SIZE, "_seq", seq, "_top.xlsx"))

# Combine the datasets into one and add a new column for category
top$Category <- paste0("Top ", seq)
bottom$Category <- paste0("Bottom ", seq)
combined_data <- bind_rows(top, bottom)

dot <- ggplot(combined_data, aes(x = evo2_delta_score, y = Category, color = Category)) +
  geom_jitter(height = 0.2, width = 0, size = 2, alpha = 0.8) +  # Jittered points
  geom_boxplot(aes(group = Category), width = 0.4, alpha = 0.2, outlier.shape = NA, color = "black") +  # Boxplot
  scale_color_manual(values = c("Top 100" = "#76B900", "Bottom 100" = "red")) +  # Custom colors
  labs(
    title = paste0("Comparing Evo 2 scores across high and low\nRovHer scores for rare variants in ", chr),
    x = "Delta Likelihood Score, Evo 2 7B Model",
    y = "RovHer scores"
  ) +
  theme_minimal() +  # Minimal theme
  plot_theme
output_file <- paste0(OUTDIR,"/", input_file, "_win", WINDOW_SIZE, "_seq", seq, "_dot.png")
ggsave(output_file, dot, width = 10, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")

#June23
####################################################################################
# 2. Distribution of Evo2 scores of two genes (LDLR and NPR2)
####################################################################################
# Variables
seq = 432
WINDOW_SIZE = 8192 # 1024,2048,4096,8192
gene = "NPR2" # LDLR NPR2

data <- read_excel(paste0("RovHer_", gene, ".txt_win", WINDOW_SIZE, "_seq",seq, "_random.xlsx"))
data$class <- factor(data$class, levels = c("FUNC/INT", "LOF"))# Ensure the 'class' column is treated as a factor

# Create the plot
plot <- ggplot(data, aes(x = evo2_delta_score, y = class, color = class)) +
  geom_jitter(height = 0.2, width = 0, size = 2, alpha = 0.8) +  # Jittered points
  geom_boxplot(aes(group = class), width = 0.4, alpha = 0.2, outlier.shape = NA, color = "black") +  # Boxplot
  scale_color_manual(values = c("FUNC/INT" = "#76B900", "LOF" = "red")) +  # Custom colors
  labs(
    title = paste0("Comparing Evo 2 scores for different RovHer-scored missense rare variants in ", gene),
    x = "Delta Likelihood Score, Evo 2 7B Model",
    y = "Variant Class"
  ) +
  theme_minimal() +  # Minimal theme
  plot_theme
output_file <- paste0(OUTDIR,"/plot/", gene, "_win", WINDOW_SIZE,"_dot.png")
ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")

#June23
####################################################################################
# 3. linear regression of RovHer vs Evo2 scores for rare variants
####################################################################################

lm_model <- lm(evo2_delta_score ~ yhat, data = combined_data)  # Regression between yhat and evo2_delta_score
r_squared <- summary(lm_model)$r.squared  # Extract R-squared
p_value <- summary(lm_model)$coefficients[2, 4]  # Extract p-value for the slope term
cat("R-squared:", r_squared, "  P-value:", p_value, "\n")

plot <- ggplot(combined_data, aes(x = yhat, y = evo2_delta_score)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Linear regression line
  labs(
    title = paste0("RovHer vs Evo2 scores for missense rare variants on Chr", chr),
    x = "RovHer score",
    y = "Evo2 delta likelihood score"
  ) +  # Axis labels
  theme_minimal() +  # Minimal theme
  xlim(0.5, 1) +  # X-axis range
  annotate(
    "text",
    x = 0.55, y = max(combined_data$evo2_delta_score, na.rm = TRUE),  # Position the annotation
    label = paste0("R² = ", round(r_squared, 3), "\nP = ", signif(p_value, 3)),
    size = 5,
    hjust = 0,
    color = "black"
  ) +
  plot_theme
output_file <- paste0(OUTDIR, "/", input_file, "_win", WINDOW_SIZE, "_seq", seq, "_linreg.png")
ggsave(output_file, plot, width = 8, height = 6, dpi = 300)



# June 24
####################################################################################
# 4. How different winindow sizes affect AUROC on Chr 17 RVs (balanced dataset)
####################################################################################
# Input dir
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June24"
OUTDIR <- paste0(DIR, "/plot")
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR, recursive = TRUE)
}

data <- read_excel(paste0("metrics.xlsx"))

x_ticks <- c(128,256,512,1024, 2048, 4096, 8192)
plot <- ggplot(data, aes(x = WINDOW_SIZE, y = AUROC)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Scatter plot points
  geom_hline(yintercept = mean(data$AUROC, na.rm = TRUE), 
             color = "red", linetype = "dotted", size = 1) +  # Mean line
  labs(
    title = "Influence of window sizes on AUROC\nEvaluated 1040 Chr17 rare variants (balanced number of LoF and missense)",
    x = "Window Size (number of nucleotides)",
    y = "AUROC (Evo2 7B model)"
  ) +
  scale_x_continuous(breaks = x_ticks) +  # Customize x-axis tick labels
  theme_minimal(base_family = font) +
  plot_theme
output_file <- paste0(OUTDIR,"/chr17_WIN_vs_AUROC.png")
ggsave(output_file, plot, width = 15, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")


####################################################################################
# Test ONE variant, e.g row=10, and compare the RovHer score for different window sizes
# row_num = 10, 21, 50,168340,64892
####################################################################################
chr <- "19"
row_num <- 64892 

# Load excel files for each WINDOW_SIZE run
final_data <- data.frame()
x_ticks <- c(128,256,512,1024, 2048, 4096, 8192) # window size

for (WINDOW_SIZE in x_ticks) {
  excel_file <- paste0("RovHer_chr",chr, "_win", WINDOW_SIZE, "_seq", row_num, "_row",row_num, ".xlsx")
  temp_data <- read_excel(excel_file)
  temp_data$WINDOW_SIZE <- WINDOW_SIZE
  final_data <- bind_rows(final_data, temp_data)
}
cat("Final Dimensions:", dim(final_data), "\n")

# Extract the first value from the "PLINK_SNP_NAME" column and store it in a new variable
plink_snp_name <- final_data$PLINK_SNP_NAME[1]
cat("plink_snp_name:", plink_snp_name, "\n")

# Calculate the mean of evo2_delta_score
mean_evo2_delta <- mean(final_data$evo2_delta_score, na.rm = TRUE)
cat("mean_evo2_delta:", mean_evo2_delta, "\n")

plot <- ggplot(final_data, aes(x = WINDOW_SIZE, y = evo2_delta_score)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Scatter plot points
  geom_hline(yintercept = mean(final_data$evo2_delta_score, na.rm = TRUE), 
             color = "red", linetype = "dotted", size = 1) +  # Mean line
  labs(
    title = paste0("Varying window sizes for a single variant (",plink_snp_name, ") on Chr",chr,"\n\nHow does it affect Evo2 scores?"),
    x = "Window Size (number of nucleotides)",
    y = "Evo2 Delta Likelihood Score"
  ) +
  scale_x_continuous(breaks = x_ticks) +  # Customize x-axis tick labels
  theme_minimal(base_family = font) +
  plot_theme

output_file <- paste0(OUTDIR,"/chr",chr,"_diff_win_size_row",row_num,".png")
ggsave(output_file, plot, width = 16, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")

# June 24 embeddings
####################################################################################
#  Plot: Absolute delta embedding scores across embedding positions
# embeddng dimensions: 1 x 4096
####################################################################################

# Input dir
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June24_embedding"
OUTDIR <- paste0(DIR, "/plot")
if (!dir.exists(OUTDIR)) {dir.create(OUTDIR, recursive = TRUE)}

# Variables
region <- "BRCA1"
var_pos_to_plot <- 10
red <- "#DE433C"
blue <- "#3C5FDE"
green <- "#30E35A"
yellow <- "#E3BD30"
c <- yellow

# Outfile 
OUTFILE <- paste0(OUTDIR, "/",region,"_delta_embedding_varpos_",var_pos_to_plot,".png")

# Load excel files for each WINDOW_SIZE run
final_data <- data.frame()
pos_list <- c(1,10,100,200, 500, 1000, 2000, 3000, 4000, 4097, 5000, 6000, 7000, 8191)

for (pos in pos_list) {
  input_name <- paste0(region,"_pos",pos,".xlsx")
  
  if (!file.exists(input_name)) {
    cat("File does not exist:", input_name, "\n")
    next
  }
  temp_data <- read_excel(input_name)
  temp_data$seq_position <- seq_len(nrow(temp_data))
  temp_data$abs_delta_embedding <- abs(temp_data$delta_embedding)
  counter <- counter + 1
  final_data <- bind_rows(final_data, temp_data)

}
c <- blue
pos <- "8191"
var_pos_to_plot <- pos
input_name <- paste0(region,"_pos",pos,".xlsx")
final_data <- read_excel(input_name)
final_data$seq_position <- seq_len(nrow(final_data))
final_data$abs_delta_embedding <- abs(final_data$delta_embedding)
cat("Final Dimensions:", dim(final_data), "\n")
OUTFILE <- paste0(OUTDIR, "/",region,"_delta_embedding_varpos_",var_pos_to_plot,".png")

plot <- ggplot(final_data, aes(x = seq_position, y = abs_delta_embedding)) +
  geom_line(color = c, size = 1) +  # Use lines instead of points
  scale_x_continuous(
    limits = c(1, 4096),  # Set x-axis range
    breaks = seq(0, 4096, by = 256),  # Set x-axis ticks every 256
    labels = seq(0, 4096, by = 256)  # Ensure tick labels match
  ) +
  scale_y_continuous(
    limits = c(0, NA),  # Set the bottom y-axis limit to 0 and adjust the upper limit automatically
    expand = expansion(mult = c(0, 0.05)),  # Add a small margin above the highest values
    labels = scales::scientific_format(digits = 2)  # Format y-axis labels in scientific notation with 2 decimals
  ) +
  labs(
    title = paste0("Input sequence length = 8192\nVariant position: ", var_pos_to_plot),
    x = "Embedding sequence position",
    y = "Delta embedding score (abs)"
  ) +
  theme_minimal(base_size = 14) +  # Minimal theme with readable font size
  theme(
    axis.title.x = element_text(hjust = 0.5),  # Center x-axis label
    axis.title.y = element_text(hjust = 1),  # Right-align y-axis label
    plot.background = element_rect(fill = "white", color = NA)  # White background
) + 
plot_theme
ggsave(OUTFILE, plot, width = 16, height = 6, dpi = 300)
cat("Plot saved to:", OUTFILE, "\n")


# June 25 embedding
####################################################################################
# Scatter plot: Does sequence length impact Evo2 delta likelihood scores?
####################################################################################
# - Compare concordance of Evo2 scores generated from a single batch run versus a multi-batch run
# 1. we selected 1000 random variants
# 2. single batch run is where SEQ_LENGTH=1000 and WIN=4096
# 3. multi batch run is where we iterate the script 20 times, each time involving a random selection of 50 variants per iteration done without replacement, so that all rows are selected exactly once

# Input dir
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June25_embedding"
setwd(DIR)

# Output dir
OUTDIR <- paste0(DIR, "/plot")
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR, recursive = TRUE)}

# Input files (single batch evo2 scores vs multi-batch scores)
WIN_SIZE <- "8192" # "4096" or "8192" <--- MODIFY!!
file1 <- paste0("RovHer_chr17_win", WIN_SIZE, "seq1000random.csv") # fixed seq length

# Read the data 
data <- read.csv(file1) # replace

p <- ggplot(data, aes(x = evo2_score_singlebatch, y = evo2_delta_score)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "solid", size = 0.5) +  # 45-degree line  geom_point(color = green, size = 2, alpha = 0.8) + 
  geom_point(color = nvidia_green, size = 2, alpha = 0.7) + 
  labs(
    title = paste0("Evo2 Single batch vs multi-batch scores"),
    subtitle = paste0("Window size: ", WIN_SIZE, "   Sequence length (nt): 1000"),
    x = "Single batch score\n(1 iteration x 1000nt seq)",
    y = "Multi-batch score\n(20 iterations x 50nt seq w/o replacement)"
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
  subtitle = element_text(size = 12, hjust = 0.5, family = font)  # Centered subtitle
)
output_file <- paste0(OUTDIR,"/single_vs_multibatch_chr17_win",WIN_SIZE,"_seq1000_random.png")
ggsave(output_file, p, width = 18, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")


####################################################################################
# Are Evo2 scores monotonous?
# Objective: Investigate the impact of window size (2048 vs 8192) 
# on the Evo2 scores for a fixed input sequence of 1000 nucleotides
####################################################################################

# Input files for diff window sizes of fixed seq length of 1000nt
file1 <- paste0("RovHer_chr17_win8192_seq1000_random.xlsx") 
file2 <- paste0("RovHer_chr17_win2048_seq1000_random.xlsx")

# Read the data 
data1 <- read_excel(file1)
data2 <- read_excel(file2)

# rename data1 evo2_delta_score column to score_win8192
data1 <- data1 %>% rename(evo2_delta_score_win8192 = evo2_delta_score)
data2 <- data2 %>% rename(evo2_delta_score_win2048 = evo2_delta_score)

# remove columns that aren't needed
data1 <- data1 %>% select(PLINK_SNP_NAME, evo2_delta_score_win8192)
data2 <- data2 %>% select(PLINK_SNP_NAME, evo2_delta_score_win2048)

# merge data1 and data2 on PLINK_SNP_NAME
data <- merge(data1, data2, by = "PLINK_SNP_NAME")

# Calculate Spearman's rank correlation
library(ggpubr)
spearman_corr <- cor(data$evo2_delta_score_win8192, data$evo2_delta_score_win2048, method = "spearman")

p <- ggplot(data, aes(x = evo2_delta_score_win8192, y = evo2_delta_score_win2048)) +
  geom_point(color = royal_blue, size = 2, alpha = 0.7) + 
  labs(
    title = paste0("Does window size (win) affects Evo2 scores for a fixed input sequence of 1000nt"),
    x = "Evo2 scores (win=8192)",
    y = "Evo2 scores (win=2048)"
  ) +
  geom_label(
    aes(x = Inf, y = Inf, label = paste0("Spearman's ρ = ", round(spearman_corr, 4))),
    hjust = 1.1, vjust = 1.1,  # Adjust alignment for top-left
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
  subtitle = element_text(size = 12, hjust = 0.5, family = font)  # Centered subtitle
)
output_file <- paste0(OUTDIR,"/win2048_vs_8192_chr17_seq1000_random.png")
ggsave(output_file, p, width = 18, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")
