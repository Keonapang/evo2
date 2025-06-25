#!/bin/bash
# To run Evo2 pilot study with different parameters
# M:\pangk\Softwares\evo2\notebooks\brca1_test\run_evo2_20250617.py

pip install matplotlib pandas seaborn scikit-learn openpyxl biopython
WINDOW_SIZE="1024"
SEQ_LENGTH="3890"
MODEL_SIZE="7b" # 1.5 minutes to download

WINDOW_SIZE="1024 2048 4096 8192"
SEQ_LENGTH="50 84 150 250 500 1000 2000 8000"
MODEL_SIZE="40b" # 33 minutes to download, model split into 4 parts (savanna_evo2_40b_base.pt.part0-3) stored in ./nemo2_evo2_40b_8k

for ws in $WINDOW_SIZE; do
  for sl in $SEQ_LENGTH; do
    echo "Running: run_evo2_20250617.py -SEQ_LENGTH $sl -WINDOW_SIZE $ws -MODEL_SIZE $MODEL_SIZE"
    python run_evo2_20250617.py --SEQ_LENGTH $sl --WINDOW_SIZE $ws --MODEL_SIZE $MODEL_SIZE
  done
done

###########################################################
# June 22-24th Exploration of Evo2 with RovHer dataset
# Scoring rare variants from the UK Biobank
###########################################################
MODEL_SIZE="7b" # 1.5 minutes to download

# Run through various combinations of parameters (WINDOW_SIZE, SEQ_LENGTH)
WINDOW_SIZE="2048 4096 8192"
SEQ_LENGTH="84 150 250 500 2000 5000 8000 10000"
input_file="RovHer_chr17.txt"
REF_CHR="17"
SUBSET_METHOD="random" # random, top, bottom, balanced

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


for ws in $WINDOW_SIZE; do
  for sl in $SEQ_LENGTH; do
    python run_evo2_rovher.py --SEQ_LENGTH $sl --WINDOW_SIZE $ws --MODEL_SIZE $MODEL_SIZE --INPUT_FILE $input_file --REF_CHR $REF_CHR --SUBSET_METHOD $SUBSET_METHOD
  done
done

cd /root/verb-workspace
#########################################
pip list | grep bionemo
bionemo-core                                 2.4.0
bionemo-esm2                                 2.4
bionemo-evo2                                 2.4
bionemo-example_model                        0.0.0
bionemo-fw                                   0.0.0
bionemo-geneformer                           2.4
bionemo-geometric                            0.0.0
bionemo-llm                                  2.4
bionemo-moco                                 0.0.2
bionemo-noodles                              0.1.0
bionemo-scdl                                 0.0.7
bionemo-size-aware-batching                  1.0.0
bionemo-testing                              2.4
bionemo-webdatamodule                        1.0.0

pip show bionemo-evo2
Name: bionemo-evo2
Version: 2.4
Summary: Library containing data preprocessing, training, and inference tooling for Evo2.
Home-page: 
Author: 
Author-email: BioNeMo Team <bionemofeedback@nvidia.com>
License: 


from bionemo.core.data.load import load
import bionemo
print(dir(bionemo))
import bionemo.evo2
print("Loading evo2...")
print(dir(bionemo.evo2))
# This indicates that bionemo.evo2 is an empty module. It does not currently expose any functionality (e.g., no classes, methods, or functions).

pip install nemo-toolkit

_bf16.tar.gz.untar# ls -lh /root/.cache/bionemo/ea4a3f5c9c26d5edc10bdc85165c090ad0ff23ac2670d4f61244f5f0d9d5e817-nemo2_evo2_1b_8k_bf16.tar.gz.untar
total 8.0K
drwxr-xr-x 2 ubuntu ubuntu 4.0K Mar 11 22:26 context
drwxr-xr-x 2 ubuntu ubuntu 4.0K Mar 11 22:26 weights

##################################################################################
# Alternative approach (June 24) Install via helical/models/evo2
#  Helical Evo2 wrapper (e.g., helical.models.evo_2), it is not part of NVIDIA BioNeMo.
##################################################################################

# instructions: https://github.com/helicalAI/helical/blob/release/helical/models/evo_2/README.md
pip install helical
pip install openpyxl
# Via the Docker image
git clone https://github.com/helicalAI/helical.git
cd helical/helical/models/evo_2
docker build -t helical_with_evo_2 .

# start new terminal (separate from root@verb-workspace:~/verb-workspace)
docker run -it --gpus all helical_with_evo_2

# In root@verb-workspace:~/verb-workspace, stop terminal
docker stop 8e26c6fd6d20
docker run -it --rm -v /root/verb-workspace:/root/verb-workspace helical_with_evo_2
ls /root/verb-workspace # Verify Access to the Host Directory:

##################################################################################

docker=8e26c6fd6d20
verb_workspace=a1af82b2c536
input_name="BRCA1"
SEQ_LENGTH="8192"
file_name="${input_name}_seq${SEQ_LENGTH}_embed.txt"

docker exec -it ${docker} ls /workspace
docker cp ${docker}:/workspace/${input_name}_seq${SEQ_LENGTH}_embed.txt /tmp/${input_name}_seq${SEQ_LENGTH}_embed.txt
docker cp /tmp/${input_name}_seq${SEQ_LENGTH}_embed.txt ${verb_workspace}:/root/verb-workspace/${input_name}_seq${SEQ_LENGTH}_embed.txt

docker exec -it ${docker} ls /workspace
file_name="${input_name}_seq${SEQ_LENGTH}_embed.npy"
docker cp ${docker}:/workspace/${file} /tmp/${file}
docker cp /tmp/${file} ${verb_workspace}:/root/verb-workspace/${file}


file="joined_embeddings.xlsx"
docker exec -it ${docker} ls /workspace
docker cp ${docker}:/workspace/${file} /tmp/${file}
docker cp /tmp/${file} ${verb_workspace}:/root/verb-workspace/${file}


# other way around 
file="joined_embeddings.xlsx"
docker exec -it ${docker} ls /workspace
docker cp ${docker}:/workspace/${file} /tmp/${file}
docker cp /tmp/${file} ${verb_workspace}:/root/verb-workspace/${file}
