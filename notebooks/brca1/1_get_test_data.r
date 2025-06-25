
library(data.table)

# Output dir 
dir <- "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data"

train_path <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Training/height_FDR/genebass_AF001_chrall_vars_gene_final.txt"  # 4926639 x  80 
dt <- fread(paste0(train_path))
colnames(dt)
dim(dt)

d <-3
np <-15 # NULL
nv <-0.1
AF <- 0.00003
upper_MAF <- "0.01"
model_dir <- "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Training/height_FDR/GENEBASS_RV_train_RV_pred_Oct10"
model_name <- paste0("cv_MARS_d3_np15_nv0.1_lowerAF3e-05_upperAF1e-02") 
PLINK_AF_LOF_yhat <- fread(paste0(model_dir,"/", model_name, "_PLINK_AF_LOF_yhat.txt"))
dim(PLINK_AF_LOF_yhat)
colnames()
head(PLINK_AF_LOF_yhat)
tail(PLINK_AF_LOF_yhat)

# Extract rows where PLINK_SNP_NAME starts with "17:"
PLINK_AF_LOF_yhat_chr <- PLINK_AF_LOF_yhat[grep("^19:", PLINK_AF_LOF_yhat$PLINK_SNP_NAME), ]
fwrite(PLINK_AF_LOF_yhat_chr, # 301074 x 5
       paste0(dir, "/RovHer_chr19.txt"), 
       sep="\t", row.names=FALSE, quote=FALSE)
    
# sort data table by the yhat column in ascending order
PLINK_AF_LOF_yhat_sorted <- PLINK_AF_LOF_yhat[order(PLINK_AF_LOF_yhat$yhat), ]


# Extract rows where Gene.refGene is exactly "NPR2"
gene_name <- "BRCA1"
gene <- PLINK_AF_LOF_yhat[PLINK_AF_LOF_yhat$Gene.refGene == gene_name, ]
write.table(gene, paste0(dir,"/", model_name, "_", gene_name,".txt"), sep="\t", row.names=FALSE, quote=FALSE)


# obtain ref genome file from NCBI
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/