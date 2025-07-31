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

MODEL="NN" # MARS or NN or RovHer
EMBED_COLS="delta" # delta, refvar embedno
LAYER="28"
ANNO_COLS="yes"
y_label="height_FDR"

if [ "$MODEL" == "MARS" ]; then
    np=200
    anno_name="yhat"
    REGION="RovHer_chr17"
    name="${REGION}_${MODEL}_d3_np${np}_nv0.1_lowerAF0e+00_anno${ANNO_COLS}_embed${EMBED_COLS}_blk${LAYER}"
    DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/${MODEL}/${EMBED_COLS}/${name}"
    DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/${MODEL}/${EMBED_COLS}/${name}"
elif [ "$MODEL" == "NN" ]; then
    EMBED_COLS="refvar" # delta, refvar embedno
    anno_name="rovher_NN_score"
    REGION="chr17"
    name="RovHer_${REGION}_blocks.${LAYER}.mlp.l3_${y_label}_anno${ANNO_COLS}"
    DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/${MODEL}/${EMBED_COLS}/${name}"
    DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/${MODEL}/${EMBED_COLS}/RovHer_${REGION}_blocks.${LAYER}.mlp.l3_height_FDR_anno${ANNO_COLS}"
else
    EMBED_COLS="embedno" # delta, refvar embedno
    anno_name="yhat"
    REGION="chr17"
    DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/${MODEL}/${EMBED_COLS}/RovHer"
    DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2/${MODEL}/${EMBED_COLS}/RovHer"
fi

########################################################################
# Step 1: Split clumped RVs into lists of top 1%, 5%, 10%, ..., 100% 
#   - input: /4_CLUMP_RESULT/.clumped RV list
#   - output: /4_CLUMP_RESULT/_top1.clumped, _top5.clumped..., _top100.clumped
#         - list of clumped RVs split by proportion of blocks
########################################################################

prop_blks=(1,5,10,15,20,25,30,35,40,50,60,70,80,90)
input_file="${DIR}/4_CLUMP_RESULT/${anno_name}_01range_CHR_17.clumped"
score_file="${DIR}/1_convert_to_pval/${anno_name}_01range.txt"
sort="ascending"
ls ${input_file}
Rscript "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/h2/split_prop_blk.r" $input_file $anno_name $prop_blks $sort $score_file

########################################################################
# 2: Get smaller num of geno blocks
#   - input: 
#         - PLINK_list consisting of clumped variant ids from /4_CLUMP_RESULT/_top1.clumped, _top5.clumped..., _top100.clumped
#         - clumped genotype matrices (all clumped RVs) from /8_GENO_0.1LD50kWIN_RDATA
#   - output: smaller set of genotype blocks for RVs in each proportion in /2_GENO_0.1LD50kWIN_RDATA
########################################################################
prop_blks=(15,10,5,1)
cores=15
Rscript /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/h2/2_geno_blks.r $DIR $input_file $prop_blks $cores


traits=("height" "LDL_direct" "Alkaline_phosphatase" "BMI" "ApoB" "ApoA" "Glycated_haemoglobin" "Triglycerides" "IGF_1" "alanine_aminotransferase" "C_reactive_protein" "Phosphate")

traits=("height" "LDL_direct" "Alkaline_phosphatase" "BMI")
traits=("ApoB" "ApoA" "Glycated_haemoglobin" "Triglycerides")
traits=("alanine_aminotransferase" "Alkaline_phosphatase")
traits=("Urate" "Cystatin_C")
# -------------------------------------------------------------
# Scripts 3 and 4 (RARity)
# -------------------------------------------------------------

root="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2"
DIR_CLUMP="${DIR_WORK}"
echo ${DIR_WORK}
script="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Keona_scripts/RARity/exome_blk"

for top in 1 5 10; do
  for trait in "${traits[@]}"; do
    cores=8
    Rscript "${script}/3_align_geno_pheno.r" $anno_name $trait $DIR_WORK $DIR_CLUMP $top $cores

    threads=2
    export TMPDIR=/tmp 
    Rscript -e 'Sys.setenv(TMPDIR = "/tmp"); source(paste0("'$script'", "/4_exome_blk_h2.r"))' $anno_name $trait $threads $DIR_WORK $top
  done
done

  # export TMPDIR=/tmp 
  # # Conversion from RData to rkyv
  # dir=${DIR_WORK}/4_exome_h2/top${top}/3_NORM_MAC2_GENO_0.1LD50kWIN_aligned_${trait}
  # if [ $(ls -1 $dir | grep ".rkyv.gz" | wc -l) -ne $(ls -1 $dir | grep ".RData" | wc -l) ]; then
  #   bash "${script}/rdata_dir_to_rkyv.sh" $dir 2 10
  # fi

# -------------------------------------------------------------
# Plot: h2 summary for top vs bottom 1000 RVs
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


# -------------------------------------------------------------
# Plot: h2 summary for top vs bottom 1000 RVs
# -------------------------------------------------------------
top_list=(5)
root="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/h2"
NN_DIR="${root}/NN/refvar/RovHer_chr17_blocks.28.mlp.l3_height_FDR_annoyes"
MARS_DIR="${root}/MARS/delta/RovHer_chr17_MARS_d3_np200_nv0.1_lowerAF0e+00_annoyes_embeddelta_blk28"
rovher_DIR="${root}/RovHer/embedno/RovHer"

Rscript "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/h2/get_h2_curve_excel_summary.r" $top_list $NN_DIR $MARS_DIR $rovher_DIR
