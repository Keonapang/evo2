

library(data.table)
suppressMessages(library(doParallel))
suppressMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
AF <- as.numeric(args[1]) # 0.00003
d <- as.numeric(args[2]) # 2 or 3
np <- as.numeric(args[3]) # 11, 13, 14
nv <- as.numeric(args[4]) # 0.1
trial_name <- args[5] # "A", "B", "C"
anno <- args[6] # "yes", "no"
cores <- args[7] # "yes", "no"


# Directory 
DIR_EVO2 <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June25_embedding"
ROOTDIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/MARS"
OUTDIR <- paste0(ROOTDIR,"/", trial_name)

# Input dir
DIR <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June23"
# Output Dir
OUTDIR <- paste0(DIR,"/plot")
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR, recursive = TRUE)
}