#!/bin/bash
# To run Evo2 pilot study with different parameters
# June - July 2025
# SCRIPT TAKES IN INPUT VARIANT FILE FROM THE SAME CHROMOSOME ONLY!

###########################################################
# Instructions:
###########################################################

# 1. Initiate NVIDIA VM by selecting launchable with H200 x 1 GPU
# 2. Install packages:
pip install matplotlib pandas seaborn scikit-learn openpyxl biopython

# 3. Set parameters and run scripts: 
#     split_varlist.py --> if you want to split the input file into smaller chunks
#     run_evo2_rovher.py --> score input seq with Evo2 scores 
#     run_evo2_multibatch_2.py


###########################################################
# Example of how to copy files in/out of docker container
# to verb workspace
###########################################################
# ---- copy INPUT files from verb workspace to local environment ----
input_files=(
    "GRCh38_chr17.fasta"
    "BRCA1_DATA.xlsx"
    "RovHer_chr17.txt"
    "run_evo2_rovher.py"
    "RovHer_chr17.txt"
    "GRCh37_chr17.fasta"
)
verb_workspace="."
docker="b5b5a0146c31" # Replace with Docker container ID
for input_file in "${input_files[@]}"; do
    echo "Copying ${input_file} to Docker container..."
    docker cp "${verb_workspace}/${input_file}" "${docker}:/workspace/bionemo2/${input_file}"
    if [ $? -eq 0 ]; then
        echo "Successfully copied ${input_file}."
    else
        echo "Failed to copy ${input_file}. Check if the file exists and try again."
    fi
done

# ---- copy results file from local environment to verb workspace ----
file="RovHer_LDLR_embedding.csv"
verb_workspace="."
docker="4dc075169e6a"  
docker cp ${docker}:/workspace/${file} ${verb_workspace}/${file}
ls ${verb_workspace}/${file}

file="metrics.xlsx"
verb_workspace="."
docker="4dc075169e6a"  
docker cp ${docker}:/workspace/${file} ${verb_workspace}/${file}
ls ${verb_workspace}/${file}

###########################################################
# June 22-24th Exploration of Evo2 with RovHer dataset
# Scoring rare variants from the UK Biobank
###########################################################
MODEL_SIZE="7b" # 1.5 minutes to download

# Run through various combinations of parameters (WINDOW_SIZE, SEQ_LENGTH)
WINDOW_SIZE="2048"
SEQ_LENGTH="10 50 100 150 250 500 2000"
input_file="RovHer_chr17.txt"
REF_CHR="17"
SUBSET_METHOD="balanced" # random, top, bottom, balanced

# Compare how differing win sizes affect performance on Chr 17 rare variants 
WINDOW_SIZE="128 256 512 1024 2048 4096 8192"
SEQ_LENGTH="1000"
input_file="RovHer_chr17.txt"
REF_CHR="17"
SUBSET_METHOD="balanced" # random, top, bottom, balanced

# Compare how differing win sizes affect performance on Chr 17 rare variants 
WINDOW_SIZE="128 256 512 1024 2048 4096 8192"
SEQ_LENGTH="500"
input_file="RovHer_BRCA1.txt"
REF_CHR="17"
SUBSET_METHOD="balanced" # random, top, bottom, balanced

# Compare distribution of RovHer score for top and bottom 100 of the SNVs in BRCA1
WINDOW_SIZE="2048"
SEQ_LENGTH="100"
input_file="RovHer_BRCA1.txt"
REF_CHR="17"
SUBSET_METHOD="top" # random, top, bottom, balanced
SUBSET_METHOD="bottom" # random, top, bottom, balanced

# Test ONE variant, e.g row=10, and compare the RovHer score for different window sizes
WINDOW_SIZE="128 256 512 1024 2048 4096 8192"
SEQ_LENGTH="64892" # row_num = SEQ_LENGTH | row_num = 10, 21, 50,168340,64892
input_file="RovHer_chr19.txt"
REF_CHR="19"
SUBSET_METHOD="row" # # row_num = 10, 21, 50,168340,64892

# --------------- June 25 -----------------

# 1. Run BALANCED dataset, varying SEQ-LENGTH (not useful)
WINDOW_SIZE="2048"
SEQ_LENGTH="10 50 100 150 250 500 1000 2000 3000"
input_file="RovHer_chr17.txt"
REF_CHR="17"
SUBSET_METHOD="balanced" # done

WINDOW_SIZE="2048"
SEQ_LENGTH="10 50 100 150 250 500 1000 2000 3000"
input_file="BRCA1_DATA.xlsx" # done

# 2. Generate Evo2 scores for all RVs in the gene 
# OUTPUT: RovHer_BRCA1_win8192_seq812_all.xlsx
WINDOW_SIZE="8192"
SEQ_LENGTH="812"
REF_CHR="17"
input_file="RovHer_BRCA1.txt"
SUBSET_METHOD="all" # done

# 3. Discussion: Evo2 single batch vs multiple batches
#   OUTPUT: RovHer_chr17_win4096seq1000random.csv and RovHer_chr17_win8192seq1000random.csv 

WINDOW_SIZE="4096 8192"
SEQ_LENGTH="1000"
input_file="RovHer_chr17.txt"
REF_CHR="17"
SUBSET_METHOD="random" 
USE_RANDOM_SEED="True" # Set to True for reproducibility, False for randomness
for ws in $WINDOW_SIZE; do
    python run_evo2_rovher.py --SEQ_LENGTH $SEQ_LENGTH --WINDOW_SIZE $ws --MODEL_SIZE "7b" --INPUT_FILE $input_file \
                              --REF_CHR $REF_CHR --SUBSET_METHOD $SUBSET_METHOD --USE_RANDOM_SEED $USE_RANDOM_SEED
done

# Pre-sets variants for this parameter: REF_CHR=17, WINDOW_SIZE=4096, SEQ_LENGTH=1000, subset_method=random
for ws in $WINDOW_SIZE; do
  python run_evo2_multibatch.py --SEQ_LENGTH $SEQ_LENGTH --WINDOW_SIZE $ws --MODEL_SIZE "7b" --INPUT_FILE $input_file \
                              --REF_CHR $REF_CHR --SUBSET_METHOD $SUBSET_METHOD --USE_RANDOM_SEED $USE_RANDOM_SEED
done

# --------------- June 25 -----------------

# 4. NO RANDOM SEED - RANDOM RESULTS! need to run, but not yet done
WINDOW_SIZE="8192" # 8192 and 2048
SEQ_LENGTH="10 50 100 150 250 500 1000 2000 3000"
input_file="RovHer_chr17.txt"
REF_CHR="17"
SUBSET_METHOD="random" 
USE_RANDOM_SEED="False"

WINDOW_SIZE="4096"
SEQ_LENGTH="10 50 100 150 250 500 1000 2000 3000"
input_file="RovHer_chr17.txt"
USE_RANDOM_SEED="False"


for ws in $WINDOW_SIZE; do
  for sl in $SEQ_LENGTH; do
    python run_evo2_rovher.py --SEQ_LENGTH $sl --WINDOW_SIZE $ws --MODEL_SIZE "7b" --INPUT_FILE $input_file \
                              --REF_CHR $REF_CHR --SUBSET_METHOD $SUBSET_METHOD --USE_RANDOM_SEED $USE_RANDOM_SEED
  done
done

# ----------------------- June 30 -------------------------

# 1. Score all variants in ONE chromosome (~30k variants)
#     - first running split_varlist.py to split the input file into smaller chunks
#     - run run_evo2_rovher.py for each chunk (RovHer_chr17_1.txt, RovHer_chr17_2.txt, ..., RovHer_chr17_101.txt)

WINDOW_SIZE="4096"
SEQ_LENGTH="1000"
INPUT_FILE="RovHer_chr17.txt"
REF_CHR="17"
SUBSET_METHOD="all" 
USE_RANDOM_SEED="True" # Set to True for reproducibility, False for randomness

chunk_size="3000" # Number of variants per chunk
python3 split_varlist.py --INPUT_FILE $INPUT_FILE --chunk_size $chunk_size

for i in {1..101}; do
    INPUT_FILE="RovHer_chr17_${i}.txt"
    SEQ_LENGTH=${chunk_size}
    python run_evo2_rovher.py --SEQ_LENGTH $SEQ_LENGTH --WINDOW_SIZE $WINDOW_SIZE --MODEL_SIZE "7b" --INPUT_FILE $INPUT_FILE \
                              --REF_CHR $REF_CHR --SUBSET_METHOD $SUBSET_METHOD --USE_RANDOM_SEED $USE_RANDOM_SEED
done
# When running on BREV cloud, you can remove certain files [CAUTION!!]
for i in {5..8}; do
    rm -v RovHer_chr17_${i}.txt
done

# This script reads in multiple Excel files (SEQ_LENGTH number of chunks) and combines them into a
# single text file by rbind() 
DIR_WORK="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/June30"
WIN="4096" # 4096
CHR="17"
SEQ_LENGTH="3000"
Rscript "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/evo2_pilot_study/rbind_evo2_scores.r" $DIR_WORK $WIN $CHR $SEQ_LENGTH

# ----------------------- July1-------------------------
# For Guilherme's analysis 
# path: /mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/evo2_AF_gene
# SCRIPT TAKES IN INPUT VARIANT FILE FROM THE SAME CHROMOSOME ONLY!
docker pull nvcr.io/nvidia/clara/bionemo-framework:2.6.2
docker run --rm -it --gpus all \
  nvcr.io/nvidia/clara/bionemo-framework:2.6.2 \
  /bin/bash
pip install matplotlib pandas seaborn scikit-learn openpyxl biopython


# TRANSFER FILES
DOCKERS="d74be1f63a8b 53edea3966c2"
verb_workspace="."
input_files=(
    "run_evo2_rovher.py"
    "ACTC1.txt"
    "ACTN2.txt"
    "ALPK3.txt"
    "BAG3.txt"
    "CSRP3.txt"
    "DES.txt"
    "DSC2.txt"
    "DSG2.txt"
    "DSP.txt"
    "FHOD3.txt"
    "FLNC.txt"
    "JPH2.txt"
    "LMNA.txt"
    "MYBPC3.txt"
    "MYH7.txt"
    "MYL2.txt"
    "MYL3.txt"
    "NEXN.txt"
    "PKP2.txt"
    "PLN.txt"
    "RBM20.txt"
    "SCN5A.txt"
    "TMEM43.txt"
    "TNNC1.txt"
    "TNNT2.txt"
    "TPM1.txt"
    "TTN.txt"
    "VCL.txt"
    "GRCh38_chr1.fasta"
    "GRCh38_chr2.fasta"
    "GRCh38_chr3.fasta"
    "GRCh38_chr6.fasta"
    "GRCh38_chr7.fasta"
    "GRCh38_chr10.fasta"
    "GRCh38_chr11.fasta"
    "GRCh38_chr12.fasta"
    "GRCh38_chr14.fasta"
    "GRCh38_chr15.fasta"
    

)
for docker in $DOCKERS; do
for file in "${input_files[@]}"; do
    docker cp "${verb_workspace}/${file}" "${docker}:/workspace/${file}" && echo "Moved ${file}." || echo "Failed to copy ${file}"
done
done

# Modify arrays for gene set 
INPUT_FILES=("ALPK3.txt" "BAG3.txt" "DES.txt" "DSC2.txt" "DSG2.txt" "DSP.txt" "FHOD3.txt" "FLNC.txt" "JPH2.txt" "LMNA.txt" "MYBPC3.txt" "MYH7.txt" "MYL2.txt" "MYL3.txt" "NEXN.txt" "PKP2.txt"
             "RBM20.txt" "SCN5A.txt" "TMEM43.txt" "TNNC1.txt" "TNNT2.txt" "TTN.txt" "VCL.txt")
REF_CHRS=("15" "10" "2" "18" "18" "6" "18" "7" "20" "1" "11" "14" "12" "3" "1" "12"
        "10" "3" "3" "3" "1" "2" "10")
# chromosomes 1,2,3,6,7,10,11,12,14,15,18,20

WINDOW_SIZE="8192"
USE_RANDOM_SEED="True"
SUBSET_METHOD="all"

for i in "${!INPUT_FILES[@]}"; do
    INPUT_FILE=${INPUT_FILES[$i]}
    python run_evo2_rovher.py --INPUT_FILE $INPUT_FILE --WINDOW_SIZE $WINDOW_SIZE --MODEL_SIZE "7b" --SUBSET_METHOD $SUBSET_METHOD --USE_RANDOM_SEED $USE_RANDOM_SEED
done
                           # --REF_CHR $REF_CHR --SUBSET_METHOD $SUBSET_METHOD --SEQ_LENGTH $SEQ_LENGTH


# ---- Move to ALL .csv from docker container to workspace ----
verb_workspace="."
for docker in $DOCKERS; do
  csv_files=$(docker exec "${docker}" bash -c "ls /workspace/*.csv 2>/dev/null")
  if [ -n "$csv_files" ]; then
      for file in $csv_files; do
          filename=$(basename "$file")
          docker cp "${docker}:${file}" "${verb_workspace}/"
          echo "Copied: $filename"
      done
  else
      echo "No .csv FOUND."
  fi
done
verb_workspace="."
for docker in $DOCKERS; do
  csv_files=$(docker exec "${docker}" bash -c "ls /workspace/*.xlsx 2>/dev/null")
  if [ -n "$csv_files" ]; then
      for file in $csv_files; do
          filename=$(basename "$file")
          docker cp "${docker}:${file}" "${verb_workspace}/"
          echo "Copied: $filename"
      done
  else
      echo "No .xlsx FOUND."
  fi
done
################################ useful code! ################################
# Objective: Check which files from $INPUT_FILES exist in directory $path 
missing_files=()
for file in "${INPUT_FILES[@]}"; do
    if [[ ! -f "$path/$file" ]]; then
        missing_files+=("$file")
    fi
done
if [[ ${#missing_files[@]} -eq 0 ]]; then
    echo "All files are present in $path."
else
    echo "The following files are missing from $path:"
    for missing_file in "${missing_files[@]}"; do
        echo "$missing_file"
    done
fi
################################################################################