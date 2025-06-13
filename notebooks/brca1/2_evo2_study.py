# Run Evo2 for pilot study
# -------- Variables ----------
# gene
# chromosome
# input varaint file
# output directory 
# how many variants
# window size
!pip install matplotlib pandas seaborn scikit-learn openpyxl

# Required imports
from Bio import SeqIO
import gzip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns
from sklearn.metrics import roc_auc_score

# Set root path
os.chdir('/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2')
# sort data table by the yhat column in ascending order
PLINK_AF_LOF_yhat_sorted <- PLINK_AF_LOF_yhat[order(PLINK_AF_LOF_yhat$yhat), ]
top <- PLINK_AF_LOF_yhat_sorted[1:100, ]
bottom <- PLINK_AF_LOF_yhat_sorted[(nrow(PLINK_AF_LOF_yhat_sorted)-99):nrow(PLINK_AF_LOF_yhat_sorted), ]

brca1_df = pd.read_excel(
    os.path.join('41586_2018_461_MOESM3_ESM.xlsx'),
    header=2, # Skip the first 2 rows (0-based index) and use the 3rd row as the header
)
brca1_df = brca1_df[[
    'chromosome', 'position (hg19)', 'reference', 'alt', 'function.score.mean', 'func.class',
]]
brca1_df.head(3)

# Extract rows where Gene.refGene is exactly "NPR2"
gene_name <- "BRCA1"
 paste0(dir,"/", model_name, "_", gene_name,".txt")
# write the top and bottom 1000 rows to a file 
dir <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Softwares/evo2/data"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
model_name <- "RovHer"
write.table(top, paste0(dir,"/", model_name, "_top_100.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(bottom, paste0(dir,"/", model_name, "_bottom_100.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(gene, paste0(dir,"/", model_name, "_", gene_name,".txt"), sep="\t", row.names=FALSE, quote=FALSE)

# Extract rows where PLINK_SNP_NAME starts with "17:"
PLINK_AF_LOF_yhat_chr17 <- PLINK_AF_LOF_yhat[grep("^17:", PLINK_AF_LOF_yhat$PLINK_SNP_NAME), ]
fwrite(PLINK_AF_LOF_yhat_chr17, 
       paste0(model_dir, "/", model_name, "_chr17.txt"), 
       sep="\t", row.names=FALSE, quote=FALSE)
    

