# More advanced script to extract REFERENCE, VARIANT SEQ and DELTA embeddings from Evo2 7B model 
# Keona Pang
# July -Aug 2025
# Variations:
    #   - Which layers in get_embeddings(dataset, embedding_layer=blk) to extract?
    #           - blocks.11.mlp.l3 to blocks.29.mlp.l3 for 7B model 
    #   - How to concatenate the embeddings?
    #           - delta + delta reverse complement 
    #           - ref + ref reverse complement + var + var reverse complement
    #   - how many embedding rows to extract surrounding the middle SNV itself? 
    #          - extract JUST the SNV itself (mid-point)
    #          - extract the surrounding context (+/- 128nt window), including the SNV itself
    #          - For EACH SNV, we will extract an embedding file (128 by 4096), which we will then average, to collapse into a vector (1 by 4096)
    #   - Output files:
    #       - delta embedding  and reverse    #   file="RovHer_${gene}_${LAYER}_delta.csv"
    #       - ref embedding and reverse       #   file="RovHer_${gene}_${LAYER}_ref.csv"
    #       - var embedding and reverse           file="RovHer_${gene}_${LAYER}_var.csv"
    #       - metrics_embedding.csv

# $VAR_WINVAR_WIN variable: add the embeddings surrounding the variant itself

# for LAYER in $LAYERS; do
#   python3.11 embedding_20250625.py \
#   --WINDOW_SIZE $WINDOW_SIZE \
#   --MODEL_SIZE "7b" \
#   --INPUT_FILE $input_file \
#   --REF_CHR $REF_CHR \
#   --LAYER $LAYER \
#   --VAR_WIN 1 # 1 or 0
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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import torch
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.metrics import auc, roc_auc_score, roc_curve
FAST_CI_MODE: bool = os.environ.get("FAST_CI_MODE", False)

########################################################################
#  Load arguments
########################################################################
parser = argparse.ArgumentParser(description="Evo2 embeddings")
parser.add_argument("--WINDOW_SIZE", type=int, required=True, default=4096, help="Window size")
parser.add_argument("--MODEL_SIZE", type=str, required=True, help="7B or 40B")
parser.add_argument("--INPUT_FILE", type=str, required=True, help="input_file")
parser.add_argument("--LAYER", required=True,type=str, help="embedding layer")
parser.add_argument("--VAR_WIN",required=True,type=int, help="how many variant embeddings to extract")
parser.add_argument("--REV", type=str, required=True, help="yes or no (generate reverse complement embeddings?)")
parser.add_argument("--REF_CHR", type=str, help="ref chromosome")

args = parser.parse_args()
WINDOW_SIZE = args.WINDOW_SIZE
MODEL_SIZE = args.MODEL_SIZE
input_file = args.INPUT_FILE
chr = args.REF_CHR
VAR_WIN = args.VAR_WIN # 1 or 0
LAYER = args.LAYER # could be none 
REV = args.REV # yes or no (generate reverse complement embeddings?)

# WINDOW_SIZE = 8192
# MODEL_SIZE = "7b"
# chr = "17"
# LAYER = "blocks.28.mlp.l3"
# VAR_WIN = 128 # 1 or higher
# REV = "yes" 
# input_file = "BRCA1_DATA.xlsx"
# input_file = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data/BRCA1_DATA.xlsx"

# Output files
input_name = os.path.splitext(input_file)[0] # Remove the file extension
region=input_name
out_metrics=Path(f"metrics_embedding.csv")

        # delta_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_delta_rev.csv")
        # ref_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_ref_rev.csv")
        # var_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_var_rev.csv")

        # delta_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_delta.csv")
        # ref_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_ref.csv")
        # var_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_var.csv")

print(f"\nembedding_advanced.py  -WIN {WINDOW_SIZE}  -LAYER {LAYER}  -MODEL_SIZE {MODEL_SIZE}  -input_file {input_file}  -VAR_WIN {VAR_WIN}\n")

################################################# 
# Import model
################################################
from helical.models.evo_2 import Evo2, Evo2Config
evo2_config = Evo2Config( # Default is model_name="evo2-7b"
    batch_size=1        
)
evo2 = Evo2(configurer=evo2_config)
# evo2.model
# print(evo2.config["model_map"])  # Check model_map
# {'model_name': 'evo2_7b', 'model_hf_name': 'arcinstitute/evo2_7b', 'default_embedding_layer': 'blocks.31.mlp.l3', 'vocab_size': 512, 'hidden_size': 4096, 'num_filters': 4096, 'hcl_layer_idxs': [2, 6, 9, 13, 16, 20, 23, 27, 30], 'hcm_layer_idxs': [1, 5, 8, 12, 15, 19, 22, 26, 29], 'hcs_layer_idxs': [0, 4, 7, 11, 14, 18, 21, 25, 28], 'attn_layer_idxs': [3, 10, 17, 24, 31], 'hcm_filter_length': 128, 'hcl_filter_groups': 4096, 'hcm_filter_groups': 256, 'hcs_filter_groups': 256, 'hcs_filter_length': 7, 'num_layers': 32, 'short_filter_length': 3, 'num_attention_heads': 32, 'short_filter_bias': False, 'mlp_init_method': 'torch.nn.init.zeros_', 'mlp_output_init_method': 'torch.nn.init.zeros_', 'eps': 1e-06, 'state_size': 16, 'rotary_emb_base': 100000000000, 'rotary_emb_scaling_factor': 128, 'use_interpolated_rotary_pos_emb': True, 'make_vocab_size_divisible_by': 8, 'inner_size_multiple_of': 16, 'inner_mlp_size': 11264, 'log_intermediate_values': False, 'proj_groups': 1, 'hyena_filter_groups': 1, 'column_split_hyena': False, 'column_split': True, 'interleave': True, 'evo2_style_activations': True, 'model_parallel_size': 1, 'pipe_parallel_size': 1, 'tie_embeddings': True, 'mha_out_proj_bias': True, 'hyena_out_proj_bias': True, 'hyena_flip_x1x2': False, 'qkv_proj_bias': False, 'use_fp8_input_projections': True, 'max_seqlen': 1048576, 'final_norm': True, 'use_flash_attn': True, 'use_flash_rmsnorm': False, 'use_flash_depthwise': False, 'use_flashfft': False, 'use_laughing_hyena': False, 'inference_mode': True, 'prefill_style': 'fft', 'mlp_activation': 'gelu', 'print_activations': False}
# for name, module in evo2.model.named_modules():
#     print(name)

########################################################################
# Define functions 
########################################################################
def initialize_metrics_file(metrics_xlsx):
    """
    Initializes the metrics file (.xlsx), adds the updated columns to the file.
    Parameters:
        metrics_xlsx (str or Path): path for the metrics Excel file.
    """
    metrics_xlsx = Path(metrics_xlsx)
    columns = [
        "MODEL_SIZE", "input_file", "WINDOW_SIZE", "VAR_WIN", "LAYER", "N_vars", 
        "embedding_shape", "time_embeddings",
    ]
    if not metrics_xlsx.exists():
        df = pd.DataFrame(columns=columns)
        df.to_excel(metrics_xlsx, index=False)
    return metrics_xlsx

def append_metrics_row(metrics_files, row_data):
    """
    Appends a row of data to the metrics files (.xlsx and .csv).
    Ensures the new columns are included in the metrics files.
    """
    metrics_xlsx = metrics_files
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

def get_reverse_complement(sequence):
    """Get the reverse complement of a DNA sequence."""
    seq_obj = Seq(sequence)  # Create a Seq object
    return str(seq_obj.reverse_complement())  # Convert back to string

################################################
# 1. Load SNV dataset from local environment
################################################

if input_file == "BRCA1_DATA.xlsx":
    data = pd.read_excel(input_file, header=2)
    data = pd.read_excel(os.path.join('BRCA1_DATA.xlsx'), header=2)
    data = data[[
        'chromosome', 'position (hg19)', 'reference', 'alt', 'function.score.mean', 'func.class',
    ]]
    data.rename(columns={
        'chromosome': 'chrom','position (hg19)': 'pos',
        'reference': 'ref','alt': 'alt',
        'function.score.mean': 'score','func.class': 'class',
    }, inplace=True)
    data['class'] = data['class'].replace(['FUNC', 'INT'], 'FUNC/INT')
    # Create a new column called 'PLINK_SNP_NAME' in the format "chrom:pos:ref:alt"
    data['PLINK_SNP_NAME'] = data.apply(
        lambda row: f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}", axis=1
    )
else:
    data = pd.read_csv(input_file, sep="\t")  # Assuming the file is tab-delimited
    # Split `PLINK_SNP_NAME` into 4 new columns: chrom, pos, ref, alt
    data[["chrom", "pos", "ref", "alt"]] = data["PLINK_SNP_NAME"].str.split(":", expand=True)
    data["pos"] = data["pos"].astype(int)
    data.rename(columns={"LOF_DISRUPTIVE": "class"}, inplace=True)

N_rows = data.shape[0]
print(f"Input: {input_file}  - Num rows: {N_rows}")
data["class"] = data["class"].replace({0: "FUNC/INT", 1: "LOF"})
print(f"{input_name} : {data.shape}")

################################################
# 2. Load reference file
################################################
chr = data.loc[0, 'chrom']
ref_file = f"GRCh38_chr{chr}.fasta"
# ref_file="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/GRChr38_ref_genome/GRCh37_chr17.fasta"
if input_file == "BRCA1_DATA.xlsx":
    ref_file = "GRCh37_chr17.fasta"

refseq = str(next(SeqIO.parse(open(ref_file, "rt"), "fasta")).seq)
print(len(refseq)) # 81195210

################################################
# 3. PARSE SEQUENCE FOR EACH VARIANT
################################################

# set working directory
# os.chdir("/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/July25_embedding_128nt")
# if embedding files alrady exist, remove it
# if delta_csv.exists():
#     delta_csv.unlink()
# if ref_csv.exists():
#     ref_csv.unlink()
# if var_csv.exists():
#     var_csv.unlink()

start_time = time.time()
if VAR_WIN == 1:

    for idx, row in data.iterrows(): # idx is the index number of the row.
    
        delta_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_delta_rev_{idx}.csv")
        ref_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_ref_rev_{idx}.csv")
        var_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_var_rev_{idx}.csv")

        delta_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_delta_{idx}.csv")
        ref_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_ref_{idx}.csv")
        var_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_var_{idx}.csv")

        row = data.iloc[idx]
        ref, var = parse_sequences(row['pos'], row['ref'], row['alt'],refseq,WINDOW_SIZE)
        print(f'[Position {idx}]  Length of ref seq: {len(ref)}    variant seq: {len(var)}') # 8192
        
        # variables for excel output file 
        input_length = len(var)
        var_pos = int((input_length / 2) + 1) # <------- get embedding for the MIDDLE position of the variant sequence
        
        # Generate embeddings 
        dataset_ref = evo2.process_data([ref])
        embedding_ref = evo2.get_embeddings(dataset_ref,LAYER) # LAYER
        dataset_var = evo2.process_data([var])
        embedding_var = evo2.get_embeddings(dataset_var,LAYER)

        # ref and var sequence embedding at the middle "variant" position
        first_embedding1 = embedding_ref['embeddings'][0]
        selected_ref = first_embedding1[var_pos]
        first_embedding = embedding_var['embeddings'][0]
        selected_var = first_embedding[var_pos]
        embed_shape = selected_var.shape
        print(f"embedding shape: {embed_shape}")

        # Calculate delta embedding
        delta_embedding = selected_ref - selected_var  # Element-wise subtraction (NumPy arrays)
        delta_df = pd.DataFrame(delta_embedding, columns=["delta_embedding"])
        ref_df = pd.DataFrame(selected_ref, columns=["ref_embedding"])
        var_df = pd.DataFrame(selected_var, columns=["var_embedding"])
        
        # Transpose final_df and reset index
        delta_df2 = delta_df.T  # Convert from [4096 rows x 1 column] to [1 row x 4096 columns]
        delta_df2.columns = [f"e{i}_delta" for i in range(1, 4097)]  # Rename columns as e1, e2, ..., e4096
        delta_df2 = delta_df2.reset_index(drop=True)  # Remove the index

        ref_df2 = ref_df.T  # Convert from [4096 rows x 1 column] to [1 row x 4096 columns]
        ref_df2.columns = [f"e{i}_ref" for i in range(1, 4097)]  # Rename columns as e1, e2, ..., e4096
        ref_df2 = ref_df2.reset_index(drop=True)  # Remove the index

        var_df2 = var_df.T  # Convert from [4096 rows x 1 column] to [1 row x 4096 columns]
        var_df2.columns = [f"e{i}_var" for i in range(1, 4097)]  # Rename columns as e1, e2, ..., e4096
        var_df2 = var_df2.reset_index(drop=True)  # Remove the index

        # Concatenate extra_columns_df with reshaped_df horizontally
        if input_file == "BRCA1_DATA.xlsx":
            extra_columns_df = pd.DataFrame({"PLINK_SNP_NAME": [row['PLINK_SNP_NAME']],"input_file": [input_name], "layer": [LAYER]})
        else:
            extra_columns_df = pd.DataFrame({"PLINK_SNP_NAME": [row['PLINK_SNP_NAME']],"input_file": [input_name],"RovHer_score": [row['yhat']], "layer": [LAYER]})

        delta_df3 = pd.concat([extra_columns_df, delta_df2], axis=1)
        ref_df3 = pd.concat([extra_columns_df, ref_df2], axis=1)
        var_df3 = pd.concat([extra_columns_df, var_df2], axis=1)
            
        # Check if the output CSV exists
        if not os.path.exists(delta_csv): 
            delta_df3.to_csv(delta_csv, index=False)
        else: 
            delta_df3.to_csv(delta_csv, mode='a', header=False, index=False)
            print(f"Data appended to {delta_csv} and _ref.csv and _var.csv")
        if not os.path.exists(ref_csv): 
            ref_df3.to_csv(ref_csv, index=False)
        else: 
            ref_df3.to_csv(ref_csv, mode='a', header=False, index=False)
        if not os.path.exists(var_csv): 
            var_df3.to_csv(var_csv, index=False)
        else: 
            var_df3.to_csv(var_csv, mode='a', header=False, index=False)   
        print(f"------------------------------------------")
        print(f" ")

        # If yes, generate the var + ref reverse complement embeddings
        if REV =="yes":
                ref_reverse = get_reverse_complement(ref)
                var_reverse = get_reverse_complement(var)
                print(f"ref: {ref[:5]}...  ref_reverse: {ref_reverse[:5]}...") 

                dataset_ref_rev = evo2.process_data([ref_reverse])
                embedding_ref2 = evo2.get_embeddings(dataset_ref_rev,LAYER) # LAYER
                dataset_var_rev = evo2.process_data([var_reverse])
                embedding_var2 = evo2.get_embeddings(dataset_var_rev,LAYER)

                first_embedding12 = embedding_ref2['embeddings'][0]
                selected_ref2 = first_embedding12[var_pos]
                first_embedding2 = embedding_var2['embeddings'][0]
                selected_var2 = first_embedding2[var_pos]
                delta_embedding2 = selected_ref2 - selected_var2  # Element-wise subtraction (NumPy arrays)
                delta_df2 = pd.DataFrame(delta_embedding2, columns=["delta_rev_embedding"])    
                ref_df2 = pd.DataFrame(selected_ref2, columns=["ref_rev_embedding"])
                var_df2 = pd.DataFrame(selected_var2, columns=["var_rev_embedding"])  

                delta_df3 = delta_df2.T  # Convert from [4096 rows x 1 column] to [1 row x 4096 columns]
                delta_df3.columns = [f"e{i}_delta_rev" for i in range(1, 4097)]  # Rename columns as e1, e2, ..., e4096
                delta_df3 = delta_df3.reset_index(drop=True)  # Remove the index

                ref_df3 = ref_df2.T  # Convert from [4096 rows x 1 column] to [1 row x 4096 columns]
                ref_df3.columns = [f"e{i}_ref_rev" for i in range(1, 4097)]  # Rename columns as e1, e2, ..., e4096
                ref_df3 = ref_df3.reset_index(drop=True)  # Remove the index

                var_df3 = var_df2.T  # Convert from [4096 rows x 1 column] to [1 row x 4096 columns]
                var_df3.columns = [f"e{i}_var_rev" for i in range(1, 4097)]  # Rename columns as e1, e2, ..., e4096
                var_df3 = var_df3.reset_index(drop=True)  # Remove the index

                delta_df4 = pd.concat([extra_columns_df, delta_df3], axis=1)
                ref_df4 = pd.concat([extra_columns_df, ref_df3], axis=1)
                var_df4 = pd.concat([extra_columns_df, var_df3], axis=1)

                if not os.path.exists(delta_rev_csv): 
                    delta_df4.to_csv(delta_rev_csv, index=False)
                else: 
                    delta_df4.to_csv(delta_rev_csv, mode='a', header=False, index=False)
                    print(f"Appended to {delta_rev_csv} and _ref_rev.csv and _var_rev.csv")
                if not os.path.exists(ref_rev_csv): 
                    ref_df4.to_csv(ref_rev_csv, index=False)
                else: 
                    ref_df4.to_csv(ref_rev_csv, mode='a', header=False, index=False)
                if not os.path.exists(var_rev_csv): 
                    var_df4.to_csv(var_rev_csv, index=False)
                else: 
                    var_df4.to_csv(var_rev_csv, mode='a', header=False, index=False)   

if VAR_WIN > 1:
    for idx, row in data.iterrows(): # idx is the index number of the row.

        delta_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_delta_rev_{idx}.csv")
        ref_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_ref_rev_{idx}.csv")
        var_rev_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_var_rev_{idx}.csv")

        delta_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_delta_{idx}.csv")
        ref_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_ref_{idx}.csv")
        var_csv=Path(f"{input_name}_{LAYER}_VARWIN{VAR_WIN}_var_{idx}.csv")

        row = data.iloc[idx]
        ref, var = parse_sequences(row['pos'], row['ref'], row['alt'],refseq,WINDOW_SIZE)
        print(f'[SNV row {idx}]  ref seq: {len(ref)} nt    variant seq: {len(var)} nt')
        #[SNV row 1]  Length of ref seq: 8192 nt   variant seq: 8192 nt

        # Get position of MIDDLE SNV in the sequence 
        input_length = len(var)
        var_pos = int((input_length / 2) + 1) # print(f"var_pos: {var_pos}   Nucleotide: {var[var_pos]}\n") # var_pos: 4097 Nucleotide: T

        # Extract the surrounding context (+/- 128nt window), including the SNV itself
        half_win = VAR_WIN // 2
        start = max(0, var_pos - half_win - 1)  # 4032; Python is 0-indexed
        end = min(len(var), var_pos + half_win - 1) # 4160
        
        ref_small = ref[start:end]
        var_small = var[start:end] # sequence of length 128nt
        print(f'[Var row {idx}]  {var_small}') # [SNV row 1]  ref seq: 128 nt    variant seq: 128 nt
        var_pos_small = int((len(var_small) / 2) + 1) # 65
        print(f"var_pos_small: {var_pos_small}   Nucleotide: {var_small[var_pos_small]}\n") # var_pos_small: 65   Nucleotide: T

        # Generate embeddings, given input ref + var sequences 
        dataset_ref = evo2.process_data([ref_small]) # (1, 1); Dataset({features: ['input_ids'],num_rows: 1})
        print(f"dataset_ref: {dataset_ref.shape}")
        embedding_ref = evo2.get_embeddings(dataset_ref,LAYER)
        dataset_var = evo2.process_data([var_small])
        embedding_var = evo2.get_embeddings(dataset_var,LAYER)

        # ref + var embedding at the central "SNV" position (var_pos)
        ref_embedding_small = embedding_ref['embeddings'][0] # (128, 4096)
        var_embedding_small = embedding_var['embeddings'][0] # (128, 4096)
        print(f"ref_embedding_small: {ref_embedding_small.shape}")
        embed_shape = ref_embedding_small.shape

        # Force the arrays to use np.float64 to ensure high precision
        ref_embedding_small = ref_embedding_small.astype(np.float64)
        var_embedding_small = var_embedding_small.astype(np.float64)

        # Compute the mean of each column
        ref_embedding_mean_vec = np.mean(ref_embedding_small, axis=0, dtype=np.float64)  # (4096,) 
        var_embedding_mean_vec = np.mean(var_embedding_small, axis=0, dtype=np.float64) 
        print(f"ref_embedding_mean_vec: {ref_embedding_mean_vec.shape}")
        
        # Reshape the vectors to 2D arrays of 1 row 
        ref_embedding_mean_vec = ref_embedding_mean_vec.reshape(1, -1) # (, 4096)
        var_embedding_mean_vec = var_embedding_mean_vec.reshape(1, -1) 
        
        # Calculate delta embedding
        delta_embedding_mean_vec = ref_embedding_mean_vec - var_embedding_mean_vec
        
        # Transpose final_df and reset index: Rename columns as e1, e2, ..., e4096
        delta_embedding_mean_vec = pd.DataFrame(delta_embedding_mean_vec, columns=[f"e{i}_delta" for i in range(1, 4097)])
        ref_embedding_mean_vec = pd.DataFrame(ref_embedding_mean_vec, columns=[f"e{i}_ref" for i in range(1, 4097)])
        var_embedding_mean_vec = pd.DataFrame(var_embedding_mean_vec, columns=[f"e{i}_var" for i in range(1, 4097)])
        
        # Concatenate extra_columns_df with reshaped_df horizontally
        if input_file == "BRCA1_DATA.xlsx":
            cols_df = pd.DataFrame({"PLINK_SNP_NAME": [row['PLINK_SNP_NAME']],"input_file": [input_name], "layer": [LAYER]})
        else:
            cols_df = pd.DataFrame({"PLINK_SNP_NAME": [row['PLINK_SNP_NAME']],"input_file": [input_name],"RovHer_score": [row['yhat']], "layer": [LAYER]})
        
        delta_df3 = pd.concat([cols_df, delta_embedding_mean_vec], axis=1)
        ref_df3 = pd.concat([cols_df, ref_embedding_mean_vec], axis=1)
        var_df3 = pd.concat([cols_df, var_embedding_mean_vec], axis=1)
        
        # Check if the output CSV exists
        if not os.path.exists(delta_csv): 
            delta_df3.to_csv(delta_csv, index=False)
        else: 
            delta_df3.to_csv(delta_csv, mode='a', header=False, index=False)
            print(f"Data appended to {delta_csv} and _ref.csv and _var.csv")
        if not os.path.exists(ref_csv): 
            ref_df3.to_csv(ref_csv, index=False)
        else: 
            ref_df3.to_csv(ref_csv, mode='a', header=False, index=False)
        if not os.path.exists(var_csv): 
            var_df3.to_csv(var_csv, index=False)
        else: 
            var_df3.to_csv(var_csv, mode='a', header=False, index=False)   
        print(f"------------------------------------------")
        print(f" ")

        # If yes, generate the var + ref reverse complement embeddings
        if REV =="yes":
                print(f"Generating reverse complement...") 
                ref_reverse = get_reverse_complement(ref_small)
                var_reverse = get_reverse_complement(ref_small)

                dataset_ref_rev = evo2.process_data([ref_reverse])
                embedding_ref2 = evo2.get_embeddings(dataset_ref_rev,LAYER) # LAYER
                dataset_var_rev = evo2.process_data([var_reverse])
                embedding_var2 = evo2.get_embeddings(dataset_var_rev,LAYER)

                ref_embedding_small = embedding_ref2['embeddings'][0]
                var_embedding_small = embedding_var2['embeddings'][0]
                embed_shape = ref_embedding_small.shape

                # Force the arrays to use np.float64 to ensure high precision
                ref_embedding_small = ref_embedding_small.astype(np.float64)
                var_embedding_small = var_embedding_small.astype(np.float64)

                # Compute the mean of each column
                ref_embedding_mean_vec = np.mean(ref_embedding_small, axis=0, dtype=np.float64)  # (4096,) 
                var_embedding_mean_vec = np.mean(var_embedding_small, axis=0, dtype=np.float64) 
                print(f"ref_embedding_mean_vec: {ref_embedding_mean_vec.shape}")
                
                # Reshape the vectors to 2D arrays of 1 row 
                ref_embedding_mean_vec = ref_embedding_mean_vec.reshape(1, -1) # (, 4096)
                var_embedding_mean_vec = var_embedding_mean_vec.reshape(1, -1) 
                
                # Calculate delta embedding
                delta_embedding_mean_vec = ref_embedding_mean_vec - var_embedding_mean_vec
                
                # Transpose final_df and reset index: Rename columns as e1, e2, ..., e4096
                delta_embedding_mean_vec = pd.DataFrame(delta_embedding_mean_vec, columns=[f"e{i}_delta_rev" for i in range(1, 4097)])
                ref_embedding_mean_vec = pd.DataFrame(ref_embedding_mean_vec, columns=[f"e{i}_ref_rev" for i in range(1, 4097)])
                var_embedding_mean_vec = pd.DataFrame(var_embedding_mean_vec, columns=[f"e{i}_var_rev" for i in range(1, 4097)])
                
                delta_df3 = pd.concat([cols_df, delta_embedding_mean_vec], axis=1)
                ref_df3 = pd.concat([cols_df, ref_embedding_mean_vec], axis=1)
                var_df3 = pd.concat([cols_df, var_embedding_mean_vec], axis=1)
                
                if not os.path.exists(delta_rev_csv): 
                    delta_df3.to_csv(delta_rev_csv, index=False)
                else: 
                    delta_df3.to_csv(delta_rev_csv, mode='a', header=False, index=False)
                    print(f"Data appended to {delta_rev_csv} and _ref_rev.csv and _var_rev.csv\n")
                if not os.path.exists(ref_rev_csv): 
                    ref_df3.to_csv(ref_rev_csv, index=False)
                else: 
                    ref_df3.to_csv(ref_rev_csv, mode='a', header=False, index=False)
                if not os.path.exists(var_rev_csv): 
                    var_df3.to_csv(var_rev_csv, index=False)
                else: 
                    var_df3.to_csv(var_rev_csv, mode='a', header=False, index=False)   

# else:
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
    "VAR_WIN": VAR_WIN,
    "LAYER": LAYER,
    "N_vars": N_rows,
    "embedding_shape": embed_shape,
    "time_embeddings": t2, 
}
append_metrics_row(metrics_file, new_row)

################################################
# abs_sum_df contains the absolute sum of each embedding col
################################################
# output_sum_csv = Path(f"{input_name}_{LAYER}_sum.csv") # not in use

# Select columns starting with "e" (e1, e2, ..., e4096)
# reshaped_df_final = pd.read_csv(output_csv)
# embedding_columns = [col for col in reshaped_df_final.columns if col.startswith("e")]

# # Calculate the absolute sum for each embedding column
# abs_sum_values = reshaped_df_final[embedding_columns].abs().sum(axis=0)

# abs_sum_df = pd.DataFrame({
#     "embedding_col": embedding_columns,
#     "abs_sum": abs_sum_values.values
# })
# abs_sum_df.to_csv(output_sum_csv, index=False)
