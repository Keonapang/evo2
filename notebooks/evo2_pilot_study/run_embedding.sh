#!/bin/bash
# July 2025
# To run Evo2 pilot study to extract embeddings from 7B/40B model 
# get_embeddings(dataset, embedding_layer=None) processes the input data and extracts embeddings

# How is the model downloaded?
#    - Helical package is a Python framework designed to simplify the use of bio-foundation models like Evo2.
#    - it provides Pre-trained models for RNA and DNA tasks (e.g., Evo2).

# Steps:
# 1. Run pip install helical to install dependencies (e.g., PyTorch, NumPy, CUDA libraries) 
# 2. Set Helical framework and create a Docker container image for the Evo2 model environment.
# 3. Generate embeddings

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

##################################################################################
# Install via hugging face (new)
##################################################################################
# Evo 2 can be run using Docker (shown below), Singularity, or Apptainer.
pip install evo2
git clone https://github.com/arcinstitute/evo2
cd evo2
pip install -e . # ERROR: Project file:///home/ubuntu/evo2 has a 'pyproject.toml' and its build backend is missing the 'build_editable' hook. Since it does not have a 'setup.py' nor a 'setup.cfg', it cannot be installed in editable mode. Consider using a build backend that supports PEP 660.
docker build -t evo2 . # Build the Docker Image
docker run -it --rm --gpus '"device=0"' -v ./huggingface:/root/.cache/huggingface evo2 bash # Allocates GPU 0 to the container, ensuring that Evo2 can utilize your NVIDIA GPU for inference or training
# change dir to inside the container
python -m evo2.test.test_evo2_generation --model_name evo2_40b # model testing script 

##################################################################################
#  Install ./helical/models/evo2 via git
# added this line of code into the DOCKERFILE of /helical/evo_2:
#   RUN git reset --hard a796302
##################################################################################
sudo rm -rf ./helical

# Via the Docker (https://github.com/helicalAI/helical/blob/release/helical/models/evo_2/README.md)
pip install helical # installs pre-req CUDA, nvidia drivers, etc...
git clone https://github.com/Keonapang/helical.git #https://github.com/Keonapang/helical/blob/release/helical/models/evo_2/Dockerfile
cd helical/helical/models/evo_2
docker build -t helical_with_evo_2 . 

# start new terminal (separate from root@verb-workspace:~/verb-workspace)
docker run -it --gpus all helical_with_evo_2

##################################################################################
# Begin analysis 
##################################################################################

# OPTIONAL: if file splitting is necessary (for entire chromosomes)
INPUT_FILE="RovHer_chr17.txt"
chunk_size="60000"
python3 split_varlist.py --INPUT_FILE $INPUT_FILE --chunk_size $chunk_size

# ---- copy INPUT files from verb workspace to local environment ----
DOCKERS="a84149fd01e3 613ce836ce3c 58d8406d4461 890897b80926 c2b46a3ea16c f6bd10d68242 f2513ea2a684"
verb_workspace="."
input_files=(
    "embedding_advanced2.py"
    "RovHer_BRCA1.txt"
    "GRCh38_chr17.fasta"
    "BRCA1_DATA.xlsx"
    "GRCh37_chr17.fasta"
    # "GRCh38_chr19.fasta"
    # "RovHer_LDLR.txt"
)
for docker in $DOCKERS; do
for file in "${input_files[@]}"; do
    docker cp "${verb_workspace}/${file}" "${docker}:/workspace/${file}" && echo "Moved ${file}." || echo "Failed to copy ${file}"
done
done

# monitoring
nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv
docker stats
watch n -1 nvidia-smi

##################################################################################
# Embedding study: July 14 2025
##################################################################################
# ----------------
BLKS="12 13 14" 
BLKS="15 16 17"
BLKS="18 19 20" 
BLKS="21 22 23" 
BLKS="24 25 26" 
BLKS="27 28 29"
input_file="RovHer_BRCA1.txt" # BRCA1_DATA.xlsx, "RovHer_LDLR.txt" RovHer_BRCA1.txt
VAR_WIN=16 # extract ONLY the middle variant (i.e., the variant of interest at the center of the window)


BLKS="12 13 14" 
BLKS="15 16 17"
BLKS="18 19 20" 
BLKS="21 22 23" 
BLKS="24 25 26" 
BLKS="27 28 29"

input_file="BRCA1_DATA.xlsx" # BRCA1_DATA.xlsx, "RovHer_LDLR.txt" RovHer_BRCA1.txt

# ----------------------------------------------
VAR_WIN=64 # extract ONLY the middle variant (i.e., the variant of interest at the center of the window)
MODEL_SIZE="7b" # 1.5 minutes to download
WINDOW_SIZE="8192" 
REF_CHR="17"
REV="yes"
pip install matplotlib pandas seaborn scikit-learn openpyxl biopython
for BLK in $BLKS; do
  python3.11 embedding_advanced2.py \
    --WINDOW_SIZE $WINDOW_SIZE \
    --MODEL_SIZE $MODEL_SIZE \
    --INPUT_FILE $input_file \
    --REF_CHR $REF_CHR \
    --LAYER "blocks.${BLK}.mlp.l3" \
    --VAR_WIN $VAR_WIN \
    --REV $REV
done


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

# ---------------------------------------------------------------------- 
file="metrics_embedding.xlsx"
verb_workspace="."
docker cp ${docker}:/workspace/${file} ./${file}
docker exec -it ${docker} ls /workspace

##################################################################################
# How to monitor GPU usage
##################################################################################
nvidia-smi # monitor GPU usage
watch -n 1 nvidia-smi # monitor GPU usage every second

docker=""
docker stats $docker

##################################################################################
# [OPTIONAL] Install via NVIDIA BioNeMo Docker image
##################################################################################
# docker pull nvcr.io/nvidia/clara/bionemo-framework:2.3
# docker run --rm -it --gpus all \
#   nvcr.io/nvidia/clara/bionemo-framework:2.3 \
#   /bin/bash
# pip install matplotlib pandas seaborn scikit-learn openpyxl biopython
