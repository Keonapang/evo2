
library(data.table)
library(ggplot2)
library(showtext)
library(sysfonts)


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

# LOAD exel file
library(readxl)
# Read the Excel file
library(openxlsx)
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# load excel spreadsheet "metrics.xlsx"
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



import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Creating the values for x and y axis
x = [1, 2, 3, 4, 5, 6, 7]
y = [1, 2, 3, 4]
x, y = np.meshgrid(x, y)

# AUC values
t1024 = np.array([0.8425, 0.886948529, 0.879032258, 0.7328, 0.744653518, 0.743377282, 0.739842661])
t2048 = np.array([0.945, 0.953125, 0.913151365, 0.8235, 0.824025794, 0.806951729, 0.804288383])
t4096 = np.array([0.9525, 0.953125, 0.910049628, 0.8346, 0.836570191, 0.811037053, 0.806454966])
t8192 = np.array([0.9525, 0.953125, 0.925248139, 0.8739, 0.873321746, 0.870751719, 0.866670692])
z = np.array([t1024, t2048, t4096, t8192])

# Create the 3D plot
fig = plt.figure(figsize=(25, 30))  # Increase figure size for more white space
ax = plt.axes(projection='3d')

# Use a surface plot for better visualization
surf = ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none')

# Customize ticks and labels
new_tick_locations = [1, 2, 3, 4, 5, 6, 7]
new_tick_labels = ['50', '84', '150', '250', '500', '1000', '2000']
ax.set_xticks(new_tick_locations)
ax.set_xticklabels(new_tick_labels, fontsize=15)  # Increase font size of tick labels

new_tick_locations2 = [1, 2, 3, 4]
new_tick_labels2 = ['1024', '2048', '4096', '8192']
ax.set_yticks(new_tick_locations2)
ax.set_yticklabels(new_tick_labels2, fontsize=15)

# Set labels and title
ax.set_xlabel('Input seq length', labelpad=25, fontsize=20)  # Add padding to move the label farther
ax.set_ylabel('Window Size', labelpad=25, fontsize=20)      # Add padding to move the label farther
ax.set_zlabel('AUC', labelpad=25, fontsize=20)              # Add padding to move the label farther
ax.set_title('AUC Performance of Evo2 7B Model', fontweight='bold', pad=30, fontsize=20)  # Bold title and padding

ax.tick_params(axis='z', labelsize=14)  # Make z-axis tick labels larger
ax.tick_params(axis='x', pad=11)        # Move x-axis tick labels farther from the plot
ax.tick_params(axis='y', pad=8)        # Move x-axis tick labels farther from the plot
ax.tick_params(axis='z', pad=11)        # Move x-axis tick labels farther from the plot

# Adjust viewing angle
ax.view_init(50, 70)
# Add a color bar to indicate AUC values with reduced size
cbar = fig.colorbar(surf, ax=ax, shrink=0.4, aspect=10)  # Shrink and reduce aspect ratio
cbar.ax.tick_params(labelsize=15)  # Smaller font size for color bar ticks

# Adjust margins to add more white space
plt.subplots_adjust(left=0.25, right=0.8, top=0.85, bottom=0.2)

# Save the plot
output_file = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/results/auc_3d_plot_adjusted.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')  # Save with high resolution
plt.show()
print(f"The plot has been saved to {output_file}")

