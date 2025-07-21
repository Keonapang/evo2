#!/bin/bash
# To run Evo2 pilot study to extract embeddings from 7B model 
# get_embeddings(dataset, embedding_layer=None) processes the input data and extracts embeddings
# which embedding layer is best?

# June - July 2025
# Brev user ID: user-2v8YcB5gWRLymdF7cLs8f7tMJVG
# Install CLI:  brew install brevdev/homebrew-brev/brev

# nvidia@inst-2zI5Uf67UkfRCzqgVZvFf80WzOV-gpu-01:~/helical/helical/models/evo_2$ 


###########################################################
# Login to brev
###########################################################

# 1. Go to browser NVIDIA Profile --> "Login via CLI (paste this in your terminal)"":
brev login --token ***** # <---------------- CHANGE THIS
brev ls
ID="jqdwu8yiy" # <---------------- CHANGE THIS

# 2. Refresh
brev refresh
curl ifconfig.me
IP="216.185.73.133" # <---------------- CHANGE THIS

# 3. Generating public/private rsa key pair (Private Key: ~/.ssh/brev_key, Public: ~/.ssh/brev_key.pub).
# Run in local terminal:
ls -l ~/.ssh/brev_key # does it exist?
ssh-keygen -t rsa -b 4096 -f ~/.ssh/brev_key -C "pangk3@mcmaster.ca"
chmod 700 ~/.ssh
chmod 600 ~/.ssh/brev_key
cat ~/.ssh/brev_key.pub # <---------- COPY THE ENTIRE LINE!!!!!!

# open another terminal window, and SSH INTO THE VM (/home/nvidia):
brev shell evo2-h100-x-1-gpu-no-container-9587ed # <---------------- CHANGE THIS 
cd ~/.ssh
mkdir ~/.ssh # (/home/nvidia/.ssh)
chmod 700 ~/.ssh
nano authorized_keys # Paste the public key you copied earlier (cat ~/.ssh/brev_key.pub) into the authorized_keys
chmod 600 authorized_keys

# Test connection
ssh -i ~/.ssh/brev_key -v nvidia@$IP
telnet $IP 22

# 2. Open a terminal to SSH into instance:
brev ssh-config


scp -i ~/.ssh/brev_key $f1 nvidia@inst-2zI5Uf67UkfRCzqgVZvFf80WzOV-gpu-01:/home/nvidia
scp -i ~/.ssh/brev_key $f1  nvidia@$IP:/home/nvidia/

###########################################################
# Instructions:
###########################################################

# 1. Initiate NVIDIA VM (no container) by selecting launchable with H200 x 1 GPU
# 2. Install packages:
pip install matplotlib pandas seaborn scikit-learn openpyxl biopython 

# 3. Upload files from local computer to the workspace
# scp -i <PATH_TO_YOUR_PRIVATE_KEY> <LOCAL_FILE_PATH> root@<INSTANCE_IP>:<REMOTE_DIRECTORY>
# scp -i ~/.ssh/brev_key myfile.txt root@<INSTANCE_IP>:/root/verb-workspace
# scp -i <PATH_TO_YOUR_PRIVATE_KEY> file1.txt file2.txt root@<INSTANCE_IP>:/root/verb-workspace
# scp -i <PATH_TO_YOUR_PRIVATE_KEY> -r my_directory root@<INSTANCE_IP>:/root/verb-workspace
cd /root/verb-workspace
ls
brev shell


DIRL="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data"
DIRF="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/GRChr38_ref_genome"
f1="$DIRL/RovHer_chr17.txt"
f2="$DIRF/GRCh38_chr17.fasta"
f3="$DIRL/RovHer_chr19.txt"
f4="$DIRF/GRCh38_chr19.fasta"

scp -i <PATH_TO_YOUR_PRIVATE_KEY> file1.txt file2.txt nvidia@inst-2zI5Uf67UkfRCzqgVZvFf80WzOV-gpu-01:/home/nvidia

scp -i ~/.ssh/brev_key $f1 $f2 $f3 $f4
ssh -i ~/.ssh/brev_key $f1 evo2-h100-x-1-gpu-no-container-844d41


scp -i ~/.ssh/brev_key $f1 $f2 $f3 $f4 nvidia@inst-2zI5Uf67UkfRCzqgVZvFf80WzOV-gpu-01:/home/nvidia
# scp -i  SHA256:L3MYphrVRmFurFydgG4zwT8z+j07SkVkKhA3UnD3vjg $f1 $f2 $f3 $f4 nvidia@inst-2zI5Uf67UkfRCzqgVZvFf80WzOV-gpu-01:/home/nvidia


# 4. Set script parameters and run: 
#     embedding_20250625.py --> to extract embeddings from Evo2 model

##################################################################################
#  Install via helical/models/evo2
#  Helical Evo2 wrapper (e.g., helical.models.evo_2), it is not part of BioNeMo.
##################################################################################

# Via the Docker (https://github.com/helicalAI/helical/blob/release/helical/models/evo_2/README.md)
pip install helical # installs pre-req CUDA, nvidia drivers, etc...
git clone https://github.com/helicalAI/helical.git
cd helical/helical/models/evo_2
docker build -t helical_with_evo_2 .

# start new terminal (separate from root@verb-workspace:~/verb-workspace)
docker run -it --gpus all helical_with_evo_2
pip install matplotlib pandas seaborn scikit-learn openpyxl biopython


##################################################################################
# Begin analysis 
##################################################################################

INPUT_FILE="RovHer_chr17.txt"
chunk_size="60000"
python3 split_varlist.py --INPUT_FILE $INPUT_FILE --chunk_size $chunk_size

# ---- copy INPUT files from verb workspace to local environment ----
DOCKERS="06cd172a5c5d 796c2db30d74 89552f31e7f7 ba5d7344594d 0b7821cc5377 22d5430b8f9f"

DOCKERS="dfa8ed72b4d6"
verb_workspace="."
input_files=(
    "embedding_advanced_2.py"
    # "embedding_advanced.py"
    # "BRCA1_DATA.xlsx"
    # "GRCh37_chr17.fasta"
    "GRCh38_chr17.fasta"
    "GRCh38_chr19.fasta"
    "RovHer_LDLR.txt"
    "RovHer_BRCA1.txt"
)

for docker in $DOCKERS; do
for file in "${input_files[@]}"; do
    docker cp "${verb_workspace}/${file}" "${docker}:/workspace/${file}" && \
    echo "Moved ${file}." || echo "Failed to copy ${file}"
done
done
docker exec -it ${docker} ls /workspace

##################################################################################
# Embedding study: July 14 2025
##################################################################################

MODEL_SIZE="7b" # 1.5 minutes to download
WINDOW_SIZE="4096"
input_file="RovHer_BRCA1.txt" # BRCA1_DATA.xlsx, "RovHer_LDLR.txt" RovHer_BRCA1.txt
REF_CHR="17"
VAR_WIN="1" # extract ONLY the middle variant (i.e., the variant of interest at the center of the window)
# LAYERS="blocks.28.mlp.l3 blocks.27.mlp.l3 blocks.31.mlp.l3" # LAYERS="blocks.22.mlp.l3 blocks.23.mlp.l3 blocks.24.mlp.l3 blocks.26.mlp.l3 blocks.27.mlp.l3 blocks.29.mlp.l3"
LAYERS="blocks.28.mlp.l3" # LAYERS="blocks.22.mlp.l3 blocks.23.mlp.l3 blocks.24.mlp.l3 blocks.26.mlp.l3 blocks.27.mlp.l3 blocks.29.mlp.l3"

for LAYER in $LAYERS; do
  python3.11 embedding_advanced_2.py \
  --WINDOW_SIZE $WINDOW_SIZE \
  --MODEL_SIZE "7b" \
  --INPUT_FILE $input_file \
  --REF_CHR $REF_CHR \
  --LAYER $LAYER \
  --VAR_WIN $VAR_WIN
done

# ---- Move to ALL .csv from docker container to workspace ----
DOCKERS="06cd172a5c5d 796c2db30d74 89552f31e7f7 ba5d7344594d 0b7821cc5377 dfa8ed72b4d6"
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

# ---------------------------------------------------------------------- 
file="metrics_embedding.xlsx"
verb_workspace="."
docker cp ${docker}:/workspace/${file} ./${file}
docker exec -it ${docker} ls /workspace

##################################################################################
# [OPTIONAL] Install via NVIDIA BioNeMo Docker image
##################################################################################
# docker pull nvcr.io/nvidia/clara/bionemo-framework:2.3
# docker run --rm -it --gpus all \
#   nvcr.io/nvidia/clara/bionemo-framework:2.3 \
#   /bin/bash
# pip install matplotlib pandas seaborn scikit-learn openpyxl biopython
