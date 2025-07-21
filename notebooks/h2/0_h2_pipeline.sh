#!/bin/bash
# OBJECTIVE: Calculate h2 explained by RVs from each genotype block (can be per-chromosome blocks or ONE block containing all RVs)
#            This script takes in LD clumped RVs, extracted as genotype blocks.
# NOTE: LD clumping MUSt be repeated if a different set of RVs were used as input.
# July 1 2025
# Steps:
#       Script 3. Create genotype matrix and align with phenotype matrix in /3_NORM_MAC2_GENO_0.1LD50kWIN_aligned_trait 
#       Script 4. Calculate adj_R2 (h2) in /4_H2_RESULTS
#       Script 5. Get h2 summary /mnt/nfs/rigenenfs/Softwares/evo2/notebooks/h2/5_summarize_h2.r

# UPDATE LMUTILS PACKAGE:
# devtools::install_github("GMELab/lmutils.r")
# file.remove("/home/pangk/R/x86_64-redhat-linux-gnu-library/3.6/00LOCK-lmutils.r-master")
# unlink("/home/pangk/R/x86_64-redhat-linux-gnu-library/3.6/00LOCK-lmutils.r-master", recursive = TRUE)
# Convert .RData to .mat.gz
# mat1$save("matrix1.mat.gz")
# mat2$save("matrix2.mat.gz")
# results <- lmutils::calculate_r2(
#                     c(matrix1.mat.gz, matrix2.mat.gz),     
#                     "pheno_outcomes.RData")

###########################################################################################################
# SCRIPT 3 & 4 - RV exome blocks
###########################################################################################################

# Variables to define
WIN="4096" # 4096
chr="17"
subset="top" # top or bottom?
variant_subset="all" # subset rows from the top; "none" means no subsetting needed

score_col="yhat" # "yhat" or "evo2_delta_score"
anno="rovher"

score_col="evo2_delta_score" # "evo2_delta_score" or "yhat"
anno="evo2"

# -------------------------------------------------------------
traits=("height" "LDL_direct" -0"Alkaline_phosphatase" "BMI" "ApoB" "ApoA" "Glycated_haemoglobin" "Triglycerides" "IGF_1" "alanine_aminotransferase" "C_reactive_protein" "Phosphate")
traits=("height" "LDL_direct" "Alkaline_phosphatase" "BMI")
traits=("ApoB" "ApoA" "Glycated_haemoglobin" "Triglycerides")
traits=("IGF_1" "alanine_aminotransferase" "C_reactive_protein" "Phosphate")

DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/${anno}_chr${chr}_win${WIN}_${variant_subset}RV_${subset}"
[ "$variant_subset" == "all" ] && DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/chr${chr}_win${WIN}_${variant_subset}"
[ -d "${DIR_WORK}" ] && echo "${DIR_WORK} exists. Please rewrite."
if ["$variant_type" == "missense" ]; then
    DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/missense/${anno}_chr${chr}_win${WIN}_${variant_subset}RV_${subset}"
    [ "$variant_subset" == "all" ] && DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/missense/chr${chr}_win${WIN}_${variant_subset}"
    [ -d "${DIR_WORK}" ] && echo "${DIR_WORK} exists. Please rewrite."
fi
DIR_CLUMP="${DIR_WORK}"
script="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Keona_scripts/RARity/exome_blk"
echo $DIR_WORK

# -------------------------------------------------------------
# Scripts 3 and 4 (RARity pipeline)
# -------------------------------------------------------------

for trait in "${traits[@]}"; do
  # Script 3
  Rscript "${script}/3_align_geno_pheno.r" $score_col $trait $DIR_WORK $DIR_CLUMP
  jobs=2
  export TMPDIR=/tmp 
  # Conversion from RData to rkyv
  dir=${DIR_WORK}/4_exome_h2/3_NORM_MAC2_GENO_0.1LD50kWIN_aligned_${trait}
  if [ $(ls -1 $dir | grep ".rkyv.gz" | wc -l) -ne $(ls -1 $dir | grep ".RData" | wc -l) ]; then
    bash "${script}/rdata_dir_to_rkyv.sh" $dir $jobs
  fi
  # Script 4: calc. h2 
  threads=1
  if [ $(ls -1 $dir | grep ".rkyv.gz" | wc -l) -eq $(ls -1 $dir | grep ".RData" | wc -l) ]; then
    Rscript -e 'Sys.setenv(TMPDIR = "/tmp"); source(paste0("'$script'", "/4_exome_blk_h2.r"))' $score_col $trait $threads $DIR_WORK
  fi
done

# -------------------------------------------------------------
# Scrip to get h2 summary
# -------------------------------------------------------------

WIN="4096" # 4096
chr="17"
subsets=("top" "bottom") # top or bottom?
variant_subsets=("1000") # top and bottom XX subset of variants 
DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/missense"

for variant_subset in "${variant_subsets[@]}"; do
for subset in "${subsets[@]}"; do
  Rscript "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/h2/5_summarize_h2.r" $DIR_WORK $WIN $chr $subset $variant_subset
done
done