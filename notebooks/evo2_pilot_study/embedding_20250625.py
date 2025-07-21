# To run Evo2 pilot study to extract embeddings from 7B model 
# get_embeddings(dataset, embedding_layer=None) processes the input data and extracts embeddings
# which embedding layer is best?

# June - July 2025

# MODEL_SIZE="7b" # 1.5 minutes to download
# WINDOW_SIZE="8192"
# input_file="RovHer_BRCA1.txt"
# REF_CHR="17"
# LAYERS="blocks.5.mlp.l3 blocks.10.mlp.l3 blocks.15.mlp.l3 blocks.20.mlp.l3 blocks.25.mlp.l3 blocks.28.mlp.l3 blocks.31.mlp.l3" # string
# LAYERS="blocks.5.mlp.l3 blocks.10.mlp.l3 blocks.15.mlp.l3 blocks.20.mlp.l3 blocks.25.mlp.l3 blocks.28.mlp.l3"

# for LAYER in $LAYERS; do
#   python3.11 embedding_20250625.py \
#   --WINDOW_SIZE $WINDOW_SIZE \
#   --MODEL_SIZE "7b" \
#   --INPUT_FILE $input_file \
#   --REF_CHR $REF_CHR \
#   --LAYER $LAYER
# done

# pip install matplotlib pandas seaborn scikit-learn openpyxl biopython

import argparse
import glob
import gzip
import json
import math
import os
from pathlib import Path
import time
import openpyxl
import subprocess
import warnings
warnings.simplefilter("ignore", FutureWarning)

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import torch
from Bio import SeqIO
from sklearn.metrics import auc, roc_auc_score, roc_curve
FAST_CI_MODE: bool = os.environ.get("FAST_CI_MODE", False)

########################################################################
#  Load arguments
########################################################################
parser = argparse.ArgumentParser(description="Evo2 embeddings")
parser.add_argument("--WINDOW_SIZE", type=int, required=True, default=4096, help="Window size")
parser.add_argument("--MODEL_SIZE", type=str, required=True, help="7B or 40B")
parser.add_argument("--INPUT_FILE", type=str, required=True, help="input_file")
parser.add_argument("--REF_CHR", type=str, required=True, help="ref chromosome")
parser.add_argument("--LAYER", type=str, help="embedding layer")
# parser.add_argument("--SUBSET_METHOD", type=str, help="SUBSET_METHOD")
# parser.add_argument("--SEQ_LENGTH", type=int, required=True, default=100, help="SEQ_LENGTH")

args = parser.parse_args()
WINDOW_SIZE = args.WINDOW_SIZE
MODEL_SIZE = args.MODEL_SIZE
input_file = args.INPUT_FILE
chr = args.REF_CHR
LAYER = args.LAYER # could be none 
# subset_method=args.SUBSET_METHOD # balanced or random; bottom or top; row 
# SEQ_LENGTH = args.SEQ_LENGTH # this is row_num if subset_method is row

# WINDOW_SIZE = 8192
# MODEL_SIZE = "7b"
# input_file = "RovHer_LDLR.txt"
# chr = "19"

# Initialize THREE results file
input_name = os.path.splitext(input_file)[0] # Remove the file extension
region=input_name

output_csv=Path(f"{input_name}_{LAYER}.csv")
out_metrics=Path(f"metrics_embedding.csv")
output_sum_csv = Path(f"{input_name}_{LAYER}_sum.csv")

print(f"embedding_20250625.py  -WIN {WINDOW_SIZE} -LAYER {LAYER} -MODEL_SIZE {MODEL_SIZE} -input_file {input_file}  -chr {chr}")

################################################# 
# Import model
################################################
from helical.models.evo_2 import Evo2, Evo2Config
evo2_config = Evo2Config(batch_size=1)  # Configure Evo2
evo2 = Evo2(configurer=evo2_config)

# inspect model 
# evo2.model
# print(evo2.config["model_map"])  # Check the contents of the model_map
# {'model_name': 'evo2_7b', 'model_hf_name': 'arcinstitute/evo2_7b', 'default_embedding_layer': 'blocks.31.mlp.l3', 'vocab_size': 512, 'hidden_size': 4096, 'num_filters': 4096, 'hcl_layer_idxs': [2, 6, 9, 13, 16, 20, 23, 27, 30], 'hcm_layer_idxs': [1, 5, 8, 12, 15, 19, 22, 26, 29], 'hcs_layer_idxs': [0, 4, 7, 11, 14, 18, 21, 25, 28], 'attn_layer_idxs': [3, 10, 17, 24, 31], 'hcm_filter_length': 128, 'hcl_filter_groups': 4096, 'hcm_filter_groups': 256, 'hcs_filter_groups': 256, 'hcs_filter_length': 7, 'num_layers': 32, 'short_filter_length': 3, 'num_attention_heads': 32, 'short_filter_bias': False, 'mlp_init_method': 'torch.nn.init.zeros_', 'mlp_output_init_method': 'torch.nn.init.zeros_', 'eps': 1e-06, 'state_size': 16, 'rotary_emb_base': 100000000000, 'rotary_emb_scaling_factor': 128, 'use_interpolated_rotary_pos_emb': True, 'make_vocab_size_divisible_by': 8, 'inner_size_multiple_of': 16, 'inner_mlp_size': 11264, 'log_intermediate_values': False, 'proj_groups': 1, 'hyena_filter_groups': 1, 'column_split_hyena': False, 'column_split': True, 'interleave': True, 'evo2_style_activations': True, 'model_parallel_size': 1, 'pipe_parallel_size': 1, 'tie_embeddings': True, 'mha_out_proj_bias': True, 'hyena_out_proj_bias': True, 'hyena_flip_x1x2': False, 'qkv_proj_bias': False, 'use_fp8_input_projections': True, 'max_seqlen': 1048576, 'final_norm': True, 'use_flash_attn': True, 'use_flash_rmsnorm': False, 'use_flash_depthwise': False, 'use_flashfft': False, 'use_laughing_hyena': False, 'inference_mode': True, 'prefill_style': 'fft', 'mlp_activation': 'gelu', 'print_activations': False}

# List all layers in the model
# for name, module in evo2.model.named_modules():
#     print(name)

########################################################################
# Define functions 
########################################################################
def initialize_metrics_file(metrics_xlsx):
    """
    Initializes the metrics file (.xlsx), adds the updated columns to the file.
    Parameters:
        metrics_xlsx (str or Path): The file path for the metrics Excel file.
    """
    # Ensure the input is a Path object for compatibility
    metrics_xlsx = Path(metrics_xlsx)
    columns = [
        "MODEL_SIZE", "input_file", "WINDOW_SIZE", "LAYER", "N_vars", 
        "embedding_shape", "time_embeddings",
    ]
    # Create a new metrics file if it doesn't exist
    if not metrics_xlsx.exists():
        df = pd.DataFrame(columns=columns)
        df.to_excel(metrics_xlsx, index=False)
    
    return metrics_xlsx
# Append results to both .xlsx and .csv files
def append_metrics_row(metrics_files, row_data):
    """
    Appends a row of data to the metrics files (.xlsx and .csv).
    Ensures the new columns are included in the metrics files.
    """
    metrics_xlsx = metrics_files
    # Update the .xlsx file
    df_xlsx = pd.read_excel(metrics_xlsx)
    df_xlsx = pd.concat([df_xlsx, pd.DataFrame([row_data])], ignore_index=True)
    df_xlsx.to_excel(metrics_xlsx, index=False)

def sample_data(df, sample_frac=1.0, balanced=True, disable=True, random_state=42):
    """Sample dataframe, optionally with balanced classes.
    """
    if disable:
        return df
    if balanced: # Get the number of rows in the dataframe
        num_rows_minor_class = math.ceil(len(df[df["class"] == "LOF"]) * sample_frac)
        return (
            pd.concat(
                [
                    df[df["class"] == "LOF"].sample(n=num_rows_minor_class, random_state=random_state),
                    df[df["class"] == "FUNC/INT"].sample(n=num_rows_minor_class, random_state=random_state),
                ]
            )
            .sample(frac=1.0, random_state=random_state)
            .reset_index(drop=True)
        )
    else: # Calculate the number of rows to sample
        return df.sample(frac=sample_frac, random_state=random_state).reset_index(drop=True)

def subset_dataframe(df, seq):
    """
    Randomly subsets the dataframe to SEQ_LENGTH number of rows.
    Returns: pandas.DataFrame - A subset of the dataframe with SEQ_LENGTH rows.
    """
    print("Number of rows to extract:", seq) 
    if seq > len(df):
        raise ValueError(f"SEQ_LENGTH ({seq}) is greater than the number of rows in the DataFrame ({len(df)}).")
    subset_df = df.sample(n=seq, random_state=42)
    print("New subset:", subset_df.shape) 
    return subset_df

def parse_sequences(pos, ref, alt, refseq, window_size=WINDOW_SIZE):
    """Parse reference and variant sequences from the reference genome sequence.
    Returns:  tuple
        (reference_sequence, variant_sequence)
    """
    p = pos - 1  # Convert to 0-indexed position
    full_seq = refseq
    ref_seq_start = max(0, p - window_size // 2)
    ref_seq_end = min(len(full_seq), p + window_size // 2)
    ref_seq = refseq[ref_seq_start:ref_seq_end]
    snv_pos_in_ref = min(window_size // 2, p)
    var_seq = ref_seq[:snv_pos_in_ref] + alt + ref_seq[snv_pos_in_ref + 1 :]
    # Sanity checks
    assert len(var_seq) == len(ref_seq)
    assert ref_seq[snv_pos_in_ref] == ref
    assert var_seq[snv_pos_in_ref] == alt
    return ref_seq, var_seq

################################################
# 1. Load reference file
################################################
ref_file = f"GRCh38_chr{chr}.fasta"
refseq = str(next(SeqIO.parse(open(ref_file, "rt"), "fasta")).seq)
print(len(refseq))

################################################
# 2. Load SNV dataset from local environment
################################################

data = pd.read_csv(input_file, sep="\t")  # Assuming the file is tab-delimited
N_rows = data.shape[0]
print(f"Input file: {input_file}  - Num rows: {N_rows}")

# Split `PLINK_SNP_NAME` into 4 new columns: chrom, pos, ref, alt
data[["chrom", "pos", "ref", "alt"]] = data["PLINK_SNP_NAME"].str.split(":", expand=True)
data["pos"] = data["pos"].astype(int)

# Rename 'LOF_DISRUPTIVE' to 'class'
data.rename(columns={"LOF_DISRUPTIVE": "class"}, inplace=True)
data["class"] = data["class"].replace({0: "FUNC/INT", 1: "LOF"})
print(f"{input_name} : {data.shape}")

################################################
# 3. PARSE SEQUENCE FOR EACH VARIANT
################################################
start_time = time.time()

# if output_csv exists, remove it
if output_csv.exists():
    print(f"Removing existing file: {output_csv}")
    output_csv.unlink()  # Remove the file

for idx, row in data.iterrows(): # idx is the index number of the row.
    row = data.iloc[idx]
    ref, var = parse_sequences(row['pos'], row['ref'], row['alt'],refseq,WINDOW_SIZE)
    print(f'[Position {idx}]  Length of ref seq: {len(ref)}    variant seq: {len(var)}') # 8192
    
    # variables for excel output file 
    input_length = len(var)
    var_pos = int((input_length / 2) + 1)
    
    # Generate embeddings with get_embeddings(dataset, embedding_layer=None) 
    dataset_ref = evo2.process_data([ref])
    embedding_ref = evo2.get_embeddings(dataset_ref,LAYER) # LAYER
    dataset_var = evo2.process_data([var])
    embedding_var = evo2.get_embeddings(dataset_var,LAYER)
    
    # ref and var sequence embedding at the middle "variant" position
    first_embedding1 = embedding_ref['embeddings'][0]
    selected_ref = first_embedding1[var_pos]
    first_embedding = embedding_var['embeddings'][0]
    selected_var = first_embedding[var_pos]
    print(f"embedding shape: {selected_var.shape}")

    # Calculate delta embedding
    delta_embedding = selected_ref - selected_var  # Element-wise subtraction (NumPy arrays)
    final_df = pd.DataFrame(delta_embedding, columns=["delta_embedding"])
    
    # Transpose final_df and reset index
    reshaped_df = final_df.T  # Convert from [4096 rows x 1 column] to [1 row x 4096 columns]
    reshaped_df.columns = [f"e{i}" for i in range(1, 4097)]  # Rename columns as e1, e2, ..., e4096
    reshaped_df = reshaped_df.reset_index(drop=True)  # Remove the index
    
    extra_columns_df = pd.DataFrame({"PLINK_SNP_NAME": [row['PLINK_SNP_NAME']],"input_file": [input_name],"RovHer_score": [row['yhat']], "layer": [LAYER]})
    
    # Concatenate extra_columns_df with reshaped_df horizontally
    reshaped_df_final = pd.concat([extra_columns_df, reshaped_df], axis=1)
    print(f"reshaped_df_final: {reshaped_df_final.shape}")
    
    # Check if the output CSV exists
    if not os.path.exists(output_csv): 
        reshaped_df_final.to_csv(output_csv, index=False)
        print(f"File {output_csv} created and data written.")
    else: 
        reshaped_df_final.to_csv(output_csv, mode='a', header=False, index=False)
        print(f"Data appended to existing file {output_csv}.")
    print(f"------------------------------------------")
    print(f" ")

end_time = time.time()
t2 = end_time - start_time    

################################################
# Write to metrics files
################################################
metrics_file = initialize_metrics_file(out_metrics)
new_row = {
    "MODEL_SIZE": MODEL_SIZE,
    "input_file": input_file,
    "WINDOW_SIZE": WINDOW_SIZE,
    "LAYER": LAYER,
    "N_vars": N_rows,
    "embedding_shape": str(selected_var.shape),
    "time_embeddings": t2, 
}
append_metrics_row(metrics_file, new_row)

################################################
# abs_sum_df contains the absolute sum of each embedding col
################################################

# # Select columns starting with "e" (e1, e2, ..., e4096)
# reshaped_df_final = pd.read_csv(output_csv)
# embedding_columns = [col for col in reshaped_df_final.columns if col.startswith("e")]

# # Calculate the absolute sum for each embedding column
# abs_sum_values = reshaped_df_final[embedding_columns].abs().sum(axis=0)

# abs_sum_df = pd.DataFrame({
#     "embedding_col": embedding_columns,
#     "abs_sum": abs_sum_values.values
# })
# abs_sum_df.to_csv(output_sum_csv, index=False)
