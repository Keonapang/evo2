#!/usr/bin/bash
# LD Clumping script

#  Purpose: LD clumping based on Evo2 scores and Rovher scores
#  Start with a score file of 1000 RVs at the top and bottom scores

# This pipeline handles the following situations:
#       - method="one_block" - RVs from all chromosomes are grouped into a SINGLE RARITY BLOCK
#       - method="per_chr" - RVs from each chromosome are grouped into separate RARITY BLOCKS
#####################################################################################################

# Variables to define
WIN="4096" # 4096
chr="17"
subset="top" # extract top or bottom XX subset of RVs
var_type="missense" # missense or "all" (all variant classes)
num_var_subset="all" # subset rows from the top; "all" means no subsetting needed - you want to evaluate the TOTAL heritability of all RVs

anno_name="evo2_delta_score" # "evo2_delta_score" or "yhat"
anno="evo2"
input_file="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30/evo2_chr${chr}_win${WIN}.txt" # all RVs in Chromosome 17

# or 
anno_name="yhat" # "evo2_delta_score" or "yhat"
anno="rovher"
input_file="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data/RovHer_chr${chr}.txt" # all RVs in Chromosome 17

sort=$([ "$subset" == "top" ] && echo "descending" || echo "ascending") # "descending"(highest to lowest), (lowest to highest;low=pathogenic)

# check $num_var_subset, and $var_type
DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/${anno}_chr${chr}_win${WIN}_${num_var_subset}RV_${subset}"
[ "$num_var_subset" == "all" ] && DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/chr${chr}_win${WIN}_${num_var_subset}"
[ -d "${DIR_WORK}" ] && echo "${DIR_WORK} exists. Please rewrite."

if ["$var_type" == "missense" ]; then
    DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/missense/${anno}_chr${chr}_win${WIN}_${num_var_subset}RV_${subset}"
    [ "$num_var_subset" == "all" ] && DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30_h2/missense/chr${chr}_win${WIN}_${num_var_subset}"
    [ -d "${DIR_WORK}" ] && echo "${DIR_WORK} exists. Please rewrite."
fi

echo ${DIR_WORK}

# Initialize
LOG_DIR="${DIR_WORK}/log"
mkdir -p "${LOG_DIR}"
root="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Keona_scripts/LD_CLUMP"

# -------------------------------------------------------------------------------------------------
# STEP 1) Convert variant scores to p-values (0-1); takes in any format input file 
own="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/h2"
Rscript ${own}/1_convert_score_to_pval.r $input_file $anno_name $sort $DIR_WORK $var_type $num_var_subset

# STEP 2) Using PLINK ids from step 1), extract BAM,BIM,FAM files with MAC=2 filter
bash ${root}/2_extract_bfiles.sh $anno_name $DIR_WORK 2 # >> "${LOG_DIR}/2_extract_bfiles_${anno_name}.txt" 2>&1 &

# STEP 3) Ensure that SNPs in 1 and 2 are in the same order (reordering 1 according to 2). Identical order= faster clumping
Rscript ${root}/3_common_pids.r $anno_name $DIR_WORK

# STEP 4) LD clump genetic variants (per chromosome!). Executes chromosome-wise clumping in parallel.
# ps aux | grep plink
chr_start=1
chr_end=22
r2=0.1 # clump parameter
nohup bash "${root}/4_clump_per_chr.sh" $anno_name $DIR_WORK $chr_start $chr_end $r2 >> ${LOG_DIR}/4_${anno_name}_chr${chr_start}_${chr_end}.txt # 2>&1 &
wait

# STEP 5) Clean LD clumped RV list; and merging them into one file 
bash ${root}/5_clean_and_merge.sh $anno_name $DIR_WORK 

# STEP 6) Split clumped variants into separate Chr and blks (.clumped_00, .clumped_01, etc.)
bash ${root}/6_split_by_blks.sh $anno_name $DIR_WORK 

# STEP 7) Genotype extraction: get individual-level (participant) data for the variants 
# OUTPUT: ${DIR_OUT}/${anno_name}_01range_CHR_${ch}.clumped_00.raw
chr_start=17
chr_end=17
method="per_chr" # "per_chr" or "one_blk" (all Chr1-22 Rvs within the same block)
bash "${root}/7_extract_geno.sh" $anno_name $DIR_WORK $chr_start $chr_end $method >> ${LOG_DIR}/7_${anno_name}_chr${chr_start}_${chr_end}.txt # 2>&1
wait

# STEP 8) Genotype conversion: Convert .raw to .RData
chr_start=1
chr_end=22
method="per_chr" # "per_chr" or "one_blk" (all Chr1-22 Rvs within the same block)
bash "${root}/8_convert_to_RData.sh" $anno_name $DIR_WORK $chr_start $chr_end $method >> ${LOG_DIR}/8_${anno_name}_chr${chr_start}_${chr_end}.txt # 2>&1
wait


# -----------------------------------------------------------
# remove directory
rm -r ${DIR_WORK}/2_BFILES
rm -r ${DIR_WORK}/3_COMMON_PIDS
rm -r ${DIR_WORK}/7_GENO_0.1LD50kWIN
