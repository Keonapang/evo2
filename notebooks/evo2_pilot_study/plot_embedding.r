
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
  plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
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


# July 2 embeddings
####################################################################################
#  Plot: Absolute delta embedding scores across embedding positions
# embeddng dimensions: 1 x 4096
####################################################################################

gene <- "BRCA1" # gene of interest

# Input dir
DIR_EVO2 <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/July1_embedding"
# output dir
OUTDIR <- paste0(DIR_EVO2, "/plot"); if (!dir.exists(OUTDIR)) {dir.create(OUTDIR, recursive = TRUE)}

# Load excel files for each WINDOW_SIZE run
final_data <- data.frame()
blks <- c(10,15,20,25,28,31)
for (blk in blks) {
  input_name <- paste0(DIR_EVO2,"/RovHer_",gene,"_blocks.", blk,".mlp.l3.csv")

  if (!file.exists(input_name)) {
    cat("File does not exist:", input_name, "\n")
    next
  }
  temp_data <- fread(input_name)
  temp_data$seq_position <- seq_len(nrow(temp_data)) #  812 4101

  # get the absolute sum of each of the 4096 embedding columns 
  embedding_cols <- grep("^e\\d+$", colnames(temp_data), value = TRUE)
  abs_sum_df <- data.table(
    embedding_col = embedding_cols, # Column names
    blk = blk, # Block number
    seq_position = seq_len(length(embedding_cols)), # Sequence position
    abs_sum = sapply(embedding_cols, function(col) sum(abs(temp_data[[col]]))) # Absolute sums
  ) # 4096 rows
  final_data <- rbind(final_data, abs_sum_df)
}

# Plot the data
plot <- ggplot(final_data, aes(x = seq_position, y = abs_sum, color = as.factor(blk), group = blk)) +
  geom_line(size = 1) +  # Plot lines for each block
  scale_x_continuous(
    limits = c(1, 4096),  # Set x-axis range
    breaks = seq(0, 4096, by = 256),  # Set x-axis ticks every 256
    labels = seq(0, 4096, by = 256)  # Ensure tick labels match
  ) +
  scale_y_continuous(
    limits = c(0, NA),  # Automatically set the upper y-axis limit
    expand = expansion(mult = c(0, 0.05)),  # Add a small margin above the highest values
    labels = scales::scientific_format(digits = 2)  # Format y-axis labels in scientific notation
  ) +
  labs(
    title = paste0("Sum of embedding column values across ", x, " RVs in ", gene),
    subtitle = paste0("Each RV has an embedding vector of 1 x 4096"),
    x = "Embedding sequence position",
    y = "Sum of absolute delta embeddings",  # Updated y-axis title
    color = "Evo2 Embedding Block"  # Legend title
  ) +
  theme_minimal(base_size = 14) +  # Minimal theme with readable font size
  theme(
    axis.title.x = element_text(hjust = 0.5),  # Center x-axis label
    axis.title.y = element_text(hjust = 0.5),  # Center y-axis label
    plot.background = element_rect(fill = "white", color = NA),  # White background
    legend.position = "right"  # Position the legend on the right
  ) + 
  plot_theme
OUTFILE <- paste0(OUTDIR, "/",gene,"_sum_abs_delta_embedding.png")
ggsave(OUTFILE, plot, width = 16, height = 6, dpi = 300)
cat("Plot saved to:", OUTFILE, "\n")


data <- fread(input_name)
x <- nrow(data)

# filter final_data to only keep rows where the value in "blk" is equal to 10
for (b in c(10,15,20,25,28,31)) {

final_data2 <- final_data[blk == b]
dim(final_data2)

# print mean of abs_sum_df column from final_data
mean_abs_sum <- mean(final_data2$abs_sum, na.rm = TRUE)
cat("Mean of absolute sums across all blocks:", mean_abs_sum, "\n")


input_name <- paste0(DIR_EVO2,"/RovHer_",gene,"_blocks.", b,".mlp.l3.csv")

plot <- ggplot(final_data2, aes(x = seq_position, y = abs_sum, color = as.factor(blk), group = blk)) +
  geom_line(size = 1) +  # Plot lines for each block
  scale_x_continuous(
    limits = c(1, 4096),  # Set x-axis range
    breaks = seq(0, 4096, by = 256),  # Set x-axis ticks every 256
    labels = seq(0, 4096, by = 256)  # Ensure tick labels match
  ) +
  scale_y_continuous(
    limits = c(0, NA),  # Automatically set the upper y-axis limit
    expand = expansion(mult = c(0, 0.05)),  # Add a small margin above the highest values
    labels = scales::scientific_format(digits = 2)  # Format y-axis labels in scientific notation
  ) +
  labs(
    title = paste0("Sum of embedding column values across ", x, " RVs in ", gene, "\n\n7B model embedding layer: ", blks),
    subtitle = paste0("Each RV has an embedding vector of 1 x 4096"),
    x = "Embedding sequence position",
    y = "Sum of absolute delta embeddings",  # Updated y-axis title
    color = "Evo2 Embedding Block"  # Legend title
  ) +
  theme_minimal(base_size = 14) +  # Minimal theme with readable font size
  theme(
    axis.title.x = element_text(hjust = 0.5),  # Center x-axis label
    axis.title.y = element_text(hjust = 0.5),  # Center y-axis label
    plot.background = element_rect(fill = "white", color = NA),  # White background
    legend.position = "right"  # Position the legend on the right
  ) + 
  plot_theme
OUTFILE <- paste0(OUTDIR, "/",gene,"_sum_abs_delta_embedding_blk",b,".png")
ggsave(OUTFILE, plot, width = 16, height = 6, dpi = 300)
cat("\nPlot saved to:", OUTFILE, "\n")
}


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



