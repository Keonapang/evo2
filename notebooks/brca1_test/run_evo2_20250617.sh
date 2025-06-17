#!/bin/bash

WINDOW_SIZE="100"
SEQ_LENGTH="5"
MODEL_SIZE="7b"

WINDOW_SIZE="1024 2048 4096 8192"
SEQ_LENGTH="10 50 84 150 250 500 1000 2000"
MODEL_SIZE="7b"

pip install matplotlib pandas seaborn scikit-learn openpyxl biopython

for ws in $WINDOW_SIZE; do
  for sl in $SEQ_LENGTH; do
    echo "Running: run_evo2_20250617.py -SEQ_LENGTH $sl -WINDOW_SIZE $ws -MODEL_SIZE $MODEL_SIZE"
    python run_evo2_20250617.py --SEQ_LENGTH $sl --WINDOW_SIZE $ws --MODEL_SIZE $MODEL_SIZE
  done
done