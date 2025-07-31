# Plot AUC (test BRCA1 set) per block for Evo 2 7B
# AUC was estimated using the NN versus MARS model
# July 27, 2025

# model="NN" # "NN" or "MARS"
# REGION="BRCA1_DATA" # "RovHer_BRCA1"  or "BRCA1_DATA"
# y_layer="class" # "class" or "clinvar"
# EMBED_COL="refvar" # "delta", "refvar", or "no"
# folder_name="July25_embedding" # MODIFY
# anno="no"
# VAR_WINS=("16", "32","64","128","256") # "128" or "256"
# reverse="yes"

# VAR_WIN="256"
# for REGION in "BRCA1_DATA" "RovHer_BRCA1"; do
# for y_layer in "clinvar" "class"; do
# Rscript /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/NN/barplot_AUC_layers.r \
# $model $REGION $y_layer $EMBED_COL $folder_name $anno $VAR_WIN $reverse
# done
# done

###############################################################################

if (!requireNamespace("extrafont", quietly = TRUE)) install.packages("extrafont")
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
suppressMessages(library(extrafont))

args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
REGION <- args[2] # "BRCA1_DATA" or "RovHer_BRCA1"
y_layer <- args[3] # "class", "clinvar"
EMBED_COL <- args[4] # "delta", "refvar", "no"
folder_name <- args[5] # "July25_embedding" or "Aug01_embedding"
anno <- args[6] # "yes" or "no"
VAR_WINS <- args[7] # "128" or "256"
reverse <- args[8] # "yes" or "no"

# model <- "NN"
# REGION <- "BRCA1_DATA"
# y_layer <- "class"
# EMBED_COL <- "refvar"
# anno <- "no"
# folder_name <- "July25_embedding" # MODIFY
# VAR_WINS <- c(16,32,64,128,256)
# reverse <- "yes"

##################################################################
# variables
##################################################################

psize <- 15
teal <- "#0B7C9A"
red <- "#DE433C"
green <- "#76B900"
font <- "sans"
plot_theme <- 
  theme_minimal() +
  theme(
  panel.background = element_rect(fill = "white", color = NA),  # White background
  plot.background = element_rect(fill = "white", color = NA),   # White plot area
  plot.margin = margin(t = 15, r = 15, b = 20, l = 20),         # Adjust plot margins
  panel.grid.minor = element_blank(),   # Removes minor gridlines
  axis.title.x = element_text(margin = margin(t = 15, b = 0.8), family = font),  # Padding for x-axis
  axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0.8), family = font),  # Padding for y-axis
  axis.text = element_text(size = 12, family = font),  # Increase axis text size
  axis.title = element_text(size = 12, face = "bold", family = font),  # Increase axis title size
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5, family = font),  # Centered, bold title
  plot.subtitle = element_text(size = 10, hjust = 0.5, family = font),  # Centered subtitle
  legend.position = "none"  # Remove the legend
)

# Convert prop_blks to a numeric vector if it is a character
if (is.character(VAR_WINS)) {
  VAR_WINS <- as.numeric(unlist(strsplit(VAR_WINS, ",")))
}

##################################################################
# Input directory 
##################################################################
server <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk"
ROOTDIR <- paste0(server,"/evo2/",model,"/",folder_name) # <------ MODIFY!

##################################################################
# Plot AUC per layer for a single VAR_WIN
##################################################################

if (length(VAR_WINS) == 1) {
  VAR_WIN <- VAR_WINS[1]

# Output model name
if (model == "MARS") {
  if (EMBED_COL == "delta" || EMBED_COL == "refvar") {
      MARS_name <- paste0(trial_name,"_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_VARWIN", VAR_WIN, "_blk", blk) 
      if (reverse == "yes") {
          MARS_name <- paste0(MARS_name,"_rev")
      }
  } else if (EMBED_COL == "no") {
      MARS_name <- paste0(trial_name, "_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_anno", anno) 
  } else {
      stop("Invalid 'EMBED_COL'. Must be either: delta, refvar, no")
  }
  cat("MARS: ", MARS_name, "\n\n")
}

# Input file (csv)
if (model == "MARS") {
    input_file <- paste0(ROOTDIR,"/D/", EMBED_COL, "/", MARS_name,"/", MARS_name)
} else {
  input_file <- paste0(ROOTDIR,"/",REGION,"_VARWIN",VAR_WIN,"_",y_layer,"_ANNO",anno,"_AUC.csv")
}
# Load input
data <- fread(input_file, header = TRUE, sep = ",", data.table = FALSE)

# Step 1: Compute the mean AUC_test for each block
mean_data <- data %>%
  group_by(LAYER) %>%
  summarize(mean_AUC_test = mean(AUC_test), .groups = "drop")

# Step 2: Extract block indices and rename LAYER
mean_data <- mean_data %>%
  mutate(Block = gsub("blocks\\.(\\d+)\\.mlp\\.l3", "\\1", LAYER)) %>% # Extract the block number
  select(Block, mean_AUC_test) # Keep only Block and mean AUC columns

if (REGION == "BRCA1_DATA") {
  title <- paste0("Findlay et al (2018) BRCA1 SNVs using ", model)
} else if (REGION == "RovHer_BRCA1") {
  title <- paste0("RovHer study BRCA1 RVs using ", model)
} else {
  stop("Invalid REGION. Must be either: BRCA1_DATA or RovHer_BRCA1")
}

if (y_layer == "clinvar") {
  label <- "ClinVar"
} else if (y_layer == "class") {
  label <- "LOF vs FUNC/INT"
} else {
  stop("Invalid y_label. Must be either: height_FDR, clinvar, class")
}

# Get the lowest AUC_test value for the y-axis lower limit
lowest_AUC <- min(mean_data$mean_AUC_test)

# Step 3: Compute the lower_cutoff
lower_cutoff <- floor(lowest_AUC / 0.05) * 0.05  # Round down to the nearest multiple of 0.05
if (lower_cutoff == lowest_AUC) {
  lower_cutoff <- lower_cutoff - 0.05  # Adjust if it's exactly on a multiple of 0.05
}

# Set the upper cutoff
upper_cutoff <- 1

# Step 3: Plot the bar chart
ggplot(mean_data, aes(x = as.numeric(Block), y = mean_AUC_test)) +
  geom_bar(stat = "identity", fill = teal, color = "black", width = 0.6) +
  labs(
    title = title,
    subtitle = paste0("Av. embed window: ", VAR_WIN, "nt , Functional annos: ", anno, ", Reverse complement: ", reverse),
    x = "Evo 2 7B Block",
    y = paste0("Mean AUC (", label, "), BRCA1 Test Set")
  ) +
  scale_x_continuous(
    breaks = as.numeric(mean_data$Block)   # Ensure all block indices are shown
  ) +
  scale_y_continuous(
    breaks = seq(lower_cutoff, upper_cutoff, by = 0.05),   # Set y-axis tick intervals explicitly
  ) +
  coord_cartesian(ylim = c(lower_cutoff, upper_cutoff)) +
  plot_theme 
output_file <- paste0(ROOTDIR, "/", REGION, "_VARWIN", VAR_WIN, "_", y_layer, "_ANNO", anno, "_AUC.png")
ggsave(output_file, width = 8, height = 6, dpi = 300)
cat("\nPlot:", output_file, "\n")
}



##################################################################
# Plot AUC for multiple VAR_WINS
##################################################################

if (length(VAR_WINS) > 1) {

for (VAR_WIN in VAR_WINS) { # VAR_WIN <- "32"

# Output model name
if (model == "MARS") {
  if (EMBED_COL == "delta" || EMBED_COL == "refvar") {
      MARS_name <- paste0(trial_name,"_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_VARWIN", VAR_WIN, "_blk", blk) 
      if (reverse == "yes") {
          MARS_name <- paste0(MARS_name,"_rev")
      }
  } else if (EMBED_COL == "no") {
      MARS_name <- paste0(trial_name, "_MARS_d",d,"_np",np,"_nv", nv,"_lowerAF",AF,"_anno", anno) 
  } else {
      stop("Invalid 'EMBED_COL'. Must be either: delta, refvar, no")
  }
  cat("MARS: ", MARS_name, "\n\n")
}

# Input file (csv)
if (model == "MARS") {
    MARSDIR <- paste0(ROOTDIR,"/D/", EMBED_COL, "/", MARS_name) # <------ MODIFY!
    input_file <- paste0(MARSDIR,"/", MARS_name)
} else {
  input_file <- paste0(ROOTDIR,"/",REGION,"_VARWIN",VAR_WIN,"_",y_layer,"_ANNO",anno,"_AUC.csv")
}
# Load input
data <- fread(input_file, header = TRUE, sep = ",", data.table = FALSE)
dim(data)

# Subset data to keep rows where the EMBED_COLS column matches the specified EMBED_COLS
data <- data[data$EMBED_COLS == EMBED_COL, ]
data <- data[data$ANNO_COLS == anno, ]

# Step 1: Compute the mean AUC_test for each block
mean_data <- data %>%
  group_by(LAYER) %>%
  summarize(mean_AUC_test = mean(AUC_test), .groups = "drop")

# Step 2: Extract block indices and rename LAYER
mean_data <- mean_data %>%
  mutate(Block = gsub("blocks\\.(\\d+)\\.mlp\\.l3", "\\1", LAYER)) %>% # Extract the block number
  select(Block, mean_AUC_test) # Keep only Block and mean AUC columns

if (REGION == "BRCA1_DATA") {
  title <- paste0("Findlay et al (2018) BRCA1 SNVs using ", model)
} else if (REGION == "RovHer_BRCA1") {
  title <- paste0("RovHer study BRCA1 RVs using ", model)
} else {
  stop("Invalid REGION. Must be either: BRCA1_DATA or RovHer_BRCA1")
}

if (y_layer == "clinvar") {
  label <- "ClinVar"
} else if (y_layer == "class") {
  label <- "LOF vs FUNC/INT"
} else {
  stop("Invalid y_label. Must be either: height_FDR, clinvar, class")
}

# Get the lowest AUC_test value for the y-axis lower limit
lowest_AUC <- min(mean_data$mean_AUC_test)

# Step 3: Compute the lower_cutoff
lower_cutoff <- floor(lowest_AUC / 0.05) * 0.05  # Round down to the nearest multiple of 0.05
if (lower_cutoff == lowest_AUC) {
  lower_cutoff <- lower_cutoff - 0.05  # Adjust if it's exactly on a multiple of 0.05
}

# Set the upper cutoff
upper_cutoff <- 1

# Step 3: Plot the bar chart
ggplot(mean_data, aes(x = as.numeric(Block), y = mean_AUC_test)) +
  geom_bar(stat = "identity", fill = teal, color = "black", width = 0.6) +
  labs(
    title = title,
    subtitle = paste0("Av. embed window: ", VAR_WIN, "nt , Functional annos: ", anno, ", Reverse complement: ", reverse),
    x = "Evo 2 7B Block",
    y = paste0("Mean AUC (", label, "), BRCA1 Test Set")
  ) +
  scale_x_continuous(
    breaks = as.numeric(mean_data$Block)   # Ensure all block indices are shown
  ) +
  scale_y_continuous(
    breaks = seq(lower_cutoff, upper_cutoff, by = 0.05),   # Set y-axis tick intervals explicitly
  ) +
  coord_cartesian(ylim = c(lower_cutoff, upper_cutoff)) +
  plot_theme 
output_file <- paste0(ROOTDIR, "/", REGION, "_VARWIN", VAR_WIN, "_", y_layer, "_ANNO", anno, "_AUC.png")
ggsave(output_file, width = 8, height = 6, dpi = 300)
cat("\nPlot:", output_file, "\n")

  }
}
