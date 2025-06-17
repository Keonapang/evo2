
library(data.table)
library(ggplot2)

# Input directory
dir <- "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data"

# Input file
file <- "metrics.xlsx"

# Output plot files
outfile1 <- "plot1.png"
outfile2 <- "plot2.png"
outfile3 <- "plot3.png"

###########################################################################
# Set dir
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE)}
setwd(dir)

# LOAD files
dt <- fread(paste0(file))
colnames(dt)
dim(dt)

###########################################################################
# PLOT 1: input sequence length vs time for different window sizes
###########################################################################

# Calculate the time (sum of score_ref and score_var)
dt$Time <- dt$score_ref + dt$score_var

# Define the custom x-axis breaks
x_breaks <- c(50, 84, 150, 250, 500, 1000, 2000)

# Create the scatter plot
ggplot(dt, aes(x = factor(SEQ_LENGTH, levels = x_breaks), y = Time, 
                    color = as.factor(WINDOW_SIZE), shape = as.factor(WINDOW_SIZE))) +
  geom_point(size = 3) +  # Scatter plot points
  geom_line(aes(group = WINDOW_SIZE), size = 1) +  # Lines connecting the points
  scale_x_discrete(name = "Input sequence length (nucleotides)") +  # Custom x-axis title and breaks
  scale_y_continuous(name = "Time (s)") +  # y-axis title
  labs(color = "Window Size", shape = "Window Size") +  # Legend titles
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis text for better readability

##save plot
ggsave(outfile1, width = 10, height = 6, dpi = 300)

###########################################################################
# PLOT 2: 
###########################################################################



###########################################################################
# PLOT 3: 
###########################################################################