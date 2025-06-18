#!/bin/bash
# To run Evo2 pilot study with different parameters
# M:\pangk\Softwares\evo2\notebooks\brca1_test\run_evo2_20250617.py

pip install matplotlib pandas seaborn scikit-learn openpyxl biopython
WINDOW_SIZE="1024 2048 4096 8192"
SEQ_LENGTH="2000 5000 8000"
MODEL_SIZE="7b" # 1.5 minutes to download

WINDOW_SIZE="1024"
SEQ_LENGTH="10"
MODEL_SIZE="40b"

WINDOW_SIZE="1024 2048 4096 8192"
SEQ_LENGTH="50 84 150 250 500 1000 2000 8000"
MODEL_SIZE="40b" # 33 minutes to download, model split into 4 parts (savanna_evo2_40b_base.pt.part0-3) stored in ./nemo2_evo2_40b_8k

for ws in $WINDOW_SIZE; do
  for sl in $SEQ_LENGTH; do
    echo "Running: run_evo2_20250617.py -SEQ_LENGTH $sl -WINDOW_SIZE $ws -MODEL_SIZE $MODEL_SIZE"
    python run_evo2_20250617.py --SEQ_LENGTH $sl --WINDOW_SIZE $ws --MODEL_SIZE $MODEL_SIZE
  done
done
