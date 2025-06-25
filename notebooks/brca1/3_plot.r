
library(data.table)
library(ggplot2)
library(showtext)
library(sysfonts)
library(readxl)
library(openxlsx)
library(dplyr)

# Read the Excel file
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# Input directory
dir <- "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/results"

# Input file
file <- "metrics.xlsx"

# Output plot files
outfile1 <- "plot1.png"
outfile2 <- "plot2.png"
outfile3 <- "plot3.png"
outfile4 <- "plot_3D.png"

###########################################################################
# Set dir
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE)}
setwd(dir)

dt <- read_excel(file)
colnames(dt)
dim(dt)

# remove rows where the SEQ_LENGTH is 3890
dt <- dt[dt$SEQ_LENGTH != 3890, ]

###########################################################################
# PLOT 1: input sequence length vs time for different window sizes
###########################################################################

# Calculate the time (sum of score_ref and score_var)
dt$Time <- dt$score_ref + dt$score_var
font <- "sans"

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
  plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Centered, bold title
)


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

###########################################################################
# PLOT AUROC 
###########################################################################

if (!requireNamespace("raster", quietly = TRUE)) {
  install.packages("raster")
}
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}
if (!requireNamespace("plot3D", quietly = TRUE)) {
  install.packages("plot3D")
}
if (!requireNamespace("rgl", quietly = TRUE)) {
  install.packages("rgl")
}
library(plotly)
library(tidyr)
library(plot3D)
library(raster)
library(rgl)
install.packages("reticulate")
library(reticulate)
virtualenv_create("r-reticulate")  # Creates a virtual environment

# Prepare the dataset for a surface plot
# Pivot the data to create a grid of z-values (AUROC) across x (SEQ_LENGTH) and y (WINDOW_SIZE)
dt_surface <- dt %>%
  pivot_wider(names_from = SEQ_LENGTH, values_from = AUROC)

# Extract the grid data
z <- as.matrix(dt_surface[, -1])  # AUROC values (z-axis)
x <- unique(dt$SEQ_LENGTH)       # SEQ_LENGTH (x-axis)
y <- unique(dt$WINDOW_SIZE)      # WINDOW_SIZE (y-axis)

# Create the 3D surface plot
plot_ly(
  x = x, y = y, z = z, 
  type = "surface",
  colorscale = "Viridis"  # Gradient color scheme
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Sequence Length"),
      yaxis = list(title = "Window Size"),
      zaxis = list(title = "AUROC")
    ),
    title = "3D Surface Plot of SEQ_LENGTH, WINDOW_SIZE, and AUROC"
  )
library(htmlwidgets)
saveWidget(plot, "3d_plot.html")
save_image(plot, "3d_plot.png")


###########################################################################
# PLOT 3: 
###########################################################################

# 3D AUC plot based on the excel file data
#
import os
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

# Set the directory
dir = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/results"

# Create the directory if it doesn't exist
if not os.path.exists(dir):
    os.makedirs(dir)

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d

# Creating the values for x and y axis
x = [1, 2, 3, 4, 5, 6, 7]
y = [1, 2, 3, 4]
x, y = np.meshgrid(x, y)

# AUC values
t1024 = np.array([0.8425, 0.886948529, 0.879032258, 0.7328, 0.744653518, 0.743377282, 0.739842661])
t2048 = np.array([0.945, 0.953125, 0.913151365, 0.8235, 0.824025794, 0.806951729, 0.804288383])
t4096 = np.array([0.9525, 0.953125, 0.910049628, 0.8346, 0.836570191, 0.811037053, 0.806454966])
t8192 = np.array([0.9525, 0.953125, 0.925248139, 0.8739, 0.873321746, 0.870751719, 0.866670692])

# t1024 =np.array([0.555555556,0.8425,0.886948529,0.879032258,0.7328,0.744653518,0.743377282,0.739842661])
# t2048 =np.array([0.888888889,0.945,0.953125,0.913151365,0.8235,0.824025794,0.806951729,0.804288383])
# t4096 =np.array([0.888888889,0.9525,0.953125,0.910049628,0.8346,0.836570191,0.811037053,0.806454966])
# t8192 =np.array([0.888888889,0.9525,0.953125,0.925248139,0.8739,0.873321746,0.870751719,0.866670692])

z = np.array([t1024, t2048, t4096, t8192])

# Create the 3D plot
fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection='3d')

# Use a surface plot for better visualization
surf = ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none')

# Customize ticks and labels
new_tick_locations = [1, 2, 3, 4, 5, 6, 7]
new_tick_labels = ['50', '84', '150', '250', '500', '1000', '2000']
ax.set_xticks(new_tick_locations)
ax.set_xticklabels(new_tick_labels)

new_tick_locations2 = [1, 2, 3, 4]
new_tick_labels2 = ['1024', '2048', '4096', '8192']
ax.set_yticks(new_tick_locations2)
ax.set_yticklabels(new_tick_labels2)

ax.set_xlabel('Input seq length')
ax.set_ylabel('Window Size')
ax.set_zlabel('AUC')
ax.set_title('AUC Performance of Evo2 7B Model')

# Adjust viewing angle
ax.view_init(50, 70)

# Add a color bar to indicate AUC values
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)

# Save the plot
output_file = "auc_3d_plot.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')  # Save with high resolution

# Display the plot
plt.show()

print(f"The plot has been saved to {output_file}")

# Save the plot before showing it
output_file = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/results/auc_3d_plot.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')  # Save with high resolution


####################################################################################
# Plotting variables
####################################################################################

OUTDIR <- paste0(DIR, "/plot")
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR, recursive = TRUE)
}
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June24"
setwd(DIR)
red <- "#ad5c5f"
grey_blue <- "#668298"
yellow <- "#F6D874"
pale_orange <- "#F9B987"
color3 <-"#C89BBD"
color4 <- "#EA7F85"
lightblue <- "#9ec1da"
color6 <-"#C7D59E"

font <- "sans"
plot_theme <- 
  theme_minimal(base_family = font) +
  theme(
  panel.background = element_rect(fill = "white", color = NA),  # White background
  plot.background = element_rect(fill = "white", color = NA),   # White plot area
  plot.margin = margin(t = 15, r = 15, b = 20, l = 20),         # Adjust plot margins
  axis.title.x = element_text(margin = margin(t = 15, b = 0.8), family = font),  # Padding for x-axis
  axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0.8), family = font),  # Padding for y-axis
  axis.text = element_text(size = 12, family = font),  # Increase axis text size
  axis.title = element_text(size = 14, face = "bold", family = font),  # Increase axis title size
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5, family = font),  # Centered, bold title
  legend.position = "none"  # Remove the legend
)

###################################################################################################
# 1. Distribution of Evo2 scores of high and low RovHer scores for Chr 17 or for a gene (e.g. BRCA1)
###################################################################################################
# Variables
seq = 100 # 50,84,100,250,500,1000,2000,5000,8000,10000
WINDOW_SIZE = 2048 # 1024,2048,4096,8192
chr="BRCA1"
input_file = paste0("RovHer_",chr)

# Input dir
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June24"
setwd(DIR)
bottom <- read_excel(paste0(input_file,"_win", WINDOW_SIZE, "_seq",seq, "_bottom.xlsx"))
top <- read_excel(paste0(input_file,"_win", WINDOW_SIZE, "_seq", seq, "_top.xlsx"))

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
output_file <- paste0(DIR,"/plot/", input_file, "_win", WINDOW_SIZE, "_seq", seq, "_dot.png")
ggsave(output_file, dot, width = 10, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")

####################################################################################
# 2. linear regression of RovHer vs Evo2 scores for rare variants
####################################################################################

# Perform linear regression
lm_model <- lm(evo2_delta_score ~ yhat, data = combined_data)  # Regression between yhat and evo2_delta_score
r_squared <- summary(lm_model)$r.squared  # Extract R-squared
p_value <- summary(lm_model)$coefficients[2, 4]  # Extract p-value for the slope term

cat("Linear regression results:\n")
cat("R-squared:", r_squared, "\n")
cat("P-value:", p_value, "\n")

# Create the plot
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
output_file <- paste0(DIR, "/plot/", input_file, "_win", WINDOW_SIZE, "_seq", seq, "_linreg.png")
ggsave(output_file, plot, width = 8, height = 6, dpi = 300)


####################################################################################
# 3. Distribution of Evo2 scores of two genes (LDLR and NPR2)
####################################################################################
# Variables
seq = 432
WINDOW_SIZE = 8192 # 1024,2048,4096,8192
gene = "NPR2" # LDLR NPR2

# Input dir
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2"
setwd(DIR)

data <- read_excel(paste0("RovHer_", gene, ".txt_win", WINDOW_SIZE, "_seq",seq, "_random.xlsx"))
# Ensure the 'class' column is treated as a factor
data$class <- factor(data$class, levels = c("FUNC/INT", "LOF"))

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
output_file <- paste0(DIR,"/plot/", gene, "_win", WINDOW_SIZE,"_dot.png")
ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")



# June 24
####################################################################################
# 4. Compare how differing win sizes affect performance on Chr 17 rare variants 
# balanced dataset
####################################################################################
# Input dir
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June24"
setwd(DIR)
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

output_file <- paste0(DIR,"/plot/chr17_WIN_vs_AUROC.png")
ggsave(output_file, plot, width = 15, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")


####################################################################################
# Test ONE variant, e.g row=10, and compare the RovHer score for different window sizes
# row_num = 10, 21, 50,168340,64892
####################################################################################
chr <- "19"
input_name <- paste0("RovHer_chr",chr)
row_num <- 64892 

# Load excel files for each WINDOW_SIZE run
final_data <- data.frame()
x_ticks <- c(128,256,512,1024, 2048, 4096, 8192)
for (WINDOW_SIZE in x_ticks) {
  excel_file <- paste0(input_name, "_win", WINDOW_SIZE, "_seq", row_num, "_row",row_num, ".xlsx")
  temp_data <- read_excel(excel_file)

  # add new column for WINDOW_SIZE
  temp_data$WINDOW_SIZE <- WINDOW_SIZE
  final_data <- bind_rows(final_data, temp_data)
}
cat("Final Dimensions:", dim(final_data), "\n")
head(final_data)

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

output_file <- paste0(DIR,"/plot/chr",chr,"_diff_win_size_row",row_num,".png")
ggsave(output_file, plot, width = 16, height = 6, dpi = 300)
cat("Plot saved to:", output_file, "\n")

