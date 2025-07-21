# Run Evo2 for pilot study
# H200 (1 and 2 GPU configurations, 144 GB each)
# H100 (2 GPU configuration, 80 GB each)

# WINDOW_SIZE="1024 2048 4096 8192"
# SEQ_LENGTH="10 50 84 150 250 500 1000 2000"
# MODEL_SIZE="7b"
# input_file="RovHer_chr17.txt"
# REF_CHR="17"

# for ws in $WINDOW_SIZE; do
#   for sl in $SEQ_LENGTH; do
#     echo "Running: run_evo2_20250617.py -SEQ_LENGTH $sl -WINDOW_SIZE $ws -MODEL_SIZE $MODEL_SIZE"
#     python run_evo2_rovher.py --SEQ_LENGTH $sl --WINDOW_SIZE $ws --MODEL_SIZE $MODEL_SIZE --INPUT_FILE $input_file --REF_CHR $REF_CHR
#   done
# done

# !pip install matplotlib pandas seaborn scikit-learn openpyxl biopython

import glob
import gzip
import json
import math
import os
from pathlib import Path
import time
import subprocess
import warnings
import argparse
warnings.simplefilter("ignore", FutureWarning)
import numpy as np
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
parser = argparse.ArgumentParser(description="Evo2 pilot study")
parser.add_argument("--SEQ_LENGTH", type=int, required=True, default=100, help="SEQ_LENGTH")
parser.add_argument("--WINDOW_SIZE", type=int, required=True, default=4096, help="Window size")
parser.add_argument("--MODEL_SIZE", type=str, required=True, default="7b", help="MODEL_SIZE")
parser.add_argument("--INPUT_FILE", type=str, required=True, help="input_file")
parser.add_argument("--REF_CHR", type=str, required=True, help="REF_CHR")
parser.add_argument("--SUBSET_METHOD", type=str, required=True, default="random", help="random, top, bottom, balanced, all")
parser.add_argument("--USE_RANDOM_SEED", type=lambda x: x.lower() == 'true', required=True, 
                    default="True", help="Set to True for reproducibility, False for randomness")
OUTPUT_DIR = "fasta_files"

args = parser.parse_args()
SEQ_LENGTH = args.SEQ_LENGTH # this is row_num if subset_method is row
WINDOW_SIZE = args.WINDOW_SIZE
MODEL_SIZE = args.MODEL_SIZE
input_file = args.INPUT_FILE
chr = args.REF_CHR
subset_method=args.SUBSET_METHOD # balanced or random; bottom or top; row 
USE_RANDOM_SEED = args.USE_RANDOM_SEED
input_name = os.path.splitext(input_file)[0] # Remove the file extension
USE_RANDOM_SEED = False

# MODEL_SIZE="7b"
# WINDOW_SIZE=2048
# SEQ_LENGTH=10
# input_file="RovHer_chr17.txt"
# chr="17"
# subset_method="random" # random, top, bottom, balanced
# OUTPUT_DIR = "fasta_files"
# USE_RANDOM_SEED= False

print(f"run_evo2_rovher.py -SEQ_LENGTH {SEQ_LENGTH}  -WINDOW_SIZE {WINDOW_SIZE}  -MODEL_SIZE {MODEL_SIZE}")
print(f"                   -input_file {input_file}  -chr {chr}  -subset_method {subset_method}")

########################################################################
# Define functions 
########################################################################
def sample_data(df, sample_frac=1.0, balanced=True, disable=True, random_state=42, use_random_state=True, SEQ_LENGTH=None):
    """
    Sample dataframe, optionally with balanced classes and fixed sequence length.
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    sample_frac : float
        Fraction of data to sample
    balanced : bool
        Whether to balance classes
    disable : bool
        Whether to disable sampling
    random_state : int
        Random seed for reproducibility
    use_random_state : bool
        When True, the random_state parameter is used for reproducibility. When False, the random_state is ignored.
    SEQ_LENGTH : int
        Total number of rows to sample; must be an even number
    Returns:pandas.DataFrame of a Sampled dataframe
    """
    if disable:
        return df
    # Ensure SEQ_LENGTH is even if provided
    if SEQ_LENGTH is not None and SEQ_LENGTH % 2 != 0:
        raise ValueError("SEQ_LENGTH must be an even number.")
    rs = random_state if use_random_state else None
    if balanced:
        # If SEQ_LENGTH is specified, derive target counts
        if SEQ_LENGTH:
            target_count = SEQ_LENGTH // 2
        else:
            target_count = math.ceil(len(df[df["class"] == "LOF"]) * sample_frac)
        # Sample LOF
        lof_samples = df[df["class"] == "LOF"].sample(
            n=min(target_count, len(df[df["class"] == "LOF"])), random_state=rs
        )
        # Sample FUNC/INT
        func_samples_needed = SEQ_LENGTH - len(lof_samples) if SEQ_LENGTH else target_count
        func_samples = df[df["class"] == "FUNC/INT"].sample(
            n=min(func_samples_needed, len(df[df["class"] == "FUNC/INT"])), random_state=rs
        )
        # Combine LOF and FUNC/INT
        sampled_df = pd.concat([lof_samples, func_samples])
        # If FUNC/INT samples are insufficient, fill remaining rows with LOF
        if SEQ_LENGTH and len(sampled_df) < SEQ_LENGTH:
            remaining_count = SEQ_LENGTH - len(sampled_df)
            extra_func_samples = df[df["class"] == "FUNC/INT"].sample(
                n=remaining_count, random_state=rs, replace=True 
            )
            sampled_df = pd.concat([sampled_df, extra_func_samples])
    else:
        # If not balanced, sample based on SEQ_LENGTH or sample_frac
        if SEQ_LENGTH:
            sampled_df = df.sample(n=SEQ_LENGTH, random_state=rs)
        else:
            sampled_df = df.sample(frac=sample_frac, random_state=rs)
    sampled_df = sampled_df.sample(frac=1.0, random_state=rs).reset_index(drop=True)
    return sampled_df

def subset_dataframe(df, seq, use_random_seed, random_seed=42):
    """
    Randomly subsets the dataframe to SEQ_LENGTH number of rows, with optional reproducibility.
    Parameters:
    -----------
    df : Input dataframe to subset.
    seq : int-  The number of rows to extract (SEQ_LENGTH).
    use_random_seed : bool, optional (default=True)
        Whether to use a fixed random seed for reproducibility.
    random_seed : int, optional (default=42)
        The random seed value to use if use_random_seed is True.
    Returns:pandas.DataFrame
        A subset of the dataframe with SEQ_LENGTH rows.
    """
    print("Number of rows to extract:", seq)
    if seq > len(df):
        raise ValueError(f"SEQ_LENGTH ({seq}) is greater than the number of rows in the DataFrame ({len(df)}).")
    # Use random seed if enabled, otherwise use no seed
    seed = random_seed if use_random_seed else None
    subset_df = df.sample(n=seq, random_state=seed)
    print(f"Random seed {'enabled' if use_random_seed else 'disabled'} (random_state={seed}).")
    print("New subset shape:", subset_df.shape)
    return subset_df

def subset_top_bottom(data, subset_method, SEQ_LENGTH):
    """
    Sorts data in descending order, tplink_snp_namehen extracts rows based on the subset_method and SEQ_LENGTH.
    Returns:
        - pd.DataFrame: subsetted data.
    """
    # Sort data by the 'yhat' column in DESCEDNING order
    sorted_data = data.sort_values(by="yhat", ascending=False)
    
    # Extract rows based on the subset_method
    if subset_method == "top":
        subset = sorted_data.head(SEQ_LENGTH)  # Get the top SEQ_LENGTH rows
    elif subset_method == "bottom":
        subset = sorted_data.tail(SEQ_LENGTH)  # Get the bottom SEQ_LENGTH rows
    else:
        raise ValueError("subset_method must be either 'top' or 'bottom'")
    
    return subset

def extract_row_by_number(df, row_num):
    """
    Extracts a specific row from the dataframe based on the provided row number (1-indexed).
    Returns:
    --------
    pandas.DataFrame
        A dataframe containing only the specified row.
    """
    if row_num < 1 or row_num > len(df):
        raise ValueError(f"ROW_NUM ({row_num}) is out of bounds for the dataframe with {len(df)} rows.")
    # Extract the specific row (convert 1-indexed to 0-indexed)
    extracted_row = df.iloc[[row_num - 1]]
    return extracted_row

def parse_sequences(pos, ref, alt, refseq, window_size=WINDOW_SIZE):
    """Parse reference and variant sequences from the reference genome sequence.
    Returns:tuple
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

def generate_fasta_files(df, refseq, output_dir="fasta_files", window_size=WINDOW_SIZE):
    """Generate FASTA files for reference and variant sequences.
    Returns:pandas.DataFrame
        Dataframe with added columns for FASTA names
    """
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Paths for output files
    ref_fasta_path = output_dir / "reference_sequences.fasta"
    var_fasta_path = output_dir / "variant_sequences.fasta"
    # Track unique sequences
    ref_sequences = set()
    var_sequences = set()
    ref_seq_to_name = {}
    # Store unique sequences with metadata for writing
    ref_entries = []
    var_entries = []
    ref_names = []
    var_names = []
    # Collect unique reference and variant sequences
    for idx, row in df.iterrows():
        ref_seq, var_seq = parse_sequences(row["pos"], row["ref"], row["alt"], refseq, window_size)
        # Add to sets to ensure uniqueness
        if ref_seq not in ref_sequences:
            ref_sequences.add(ref_seq)
            ref_name = f"ref_pos_{row['pos']}_{row['ref']}_class_{row['class']}"
            ref_entries.append(f">{ref_name}\n{ref_seq}\n")
            ref_names.append(ref_name)
            ref_seq_to_name[ref_seq] = ref_name
        else:
            ref_name = ref_seq_to_name[ref_seq]
            ref_names.append(ref_name)
        if var_seq not in var_sequences:
            var_sequences.add(var_seq)
            var_name = f"var_pos_{row['pos']}_{row['ref']}to{row['alt']}_class_{row['class']}"
            var_entries.append(f">{var_name}\n{var_seq}\n")
            var_names.append(var_name)
        else:
            assert False, "Duplicate variant sequence"
    # Write unique sequences to FASTA files
    with open(ref_fasta_path, "w") as f:
        f.writelines(ref_entries)
    with open(var_fasta_path, "w") as f:
        f.writelines(var_entries)
    # Add FASTA names to dataframe
    df_with_names = df.copy()
    df_with_names["ref_fasta_name"] = ref_names
    df_with_names["var_fasta_name"] = var_names
    print(f"Total unique reference sequences: {len(ref_sequences)}")
    print(f"Total unique variant sequences: {len(var_sequences)}")
    return df_with_names

def check_fp8_support():
    """Check if FP8 is supported on the current GPU.
    FP8 requires compute capability 8.9+ (Ada Lovelace/Hopper architecture or newer).
    """
    if not torch.cuda.is_available():
        return False, "CUDA not available"
    device_props = torch.cuda.get_device_properties(0)
    compute_capability = f"{device_props.major}.{device_props.minor}"
    device_name = device_props.name
    # FP8 is supported on compute capability 8.9+ (Ada Lovelace/Hopper architecture)
    is_supported = (device_props.major > 8) or (device_props.major == 8 and device_props.minor >= 9)
    return is_supported, f"Device: {device_name}, Compute Capability: {compute_capability}"

################################################
# 1. Load reference file
################################################

ref_file = f"GRCh38_chr{chr}.fasta"

if input_file == "BRCA1_DATA.xlsx":
    ref_file = "GRCh37_chr17.fasta"

refseq = str(next(SeqIO.parse(open(ref_file, "rt"), "fasta")).seq)
print(f"ref_file: {ref_file} with length {len(refseq)}")
print(f" ")

################################################
# 2. Load SNV dataset from local environment
################################################

# Pre-selected variants from previous run
selected_vars = pd.read_excel(f"RovHer_chr17_win4096_seq{SEQ_LENGTH}_random.xlsx")
selected_vars = selected_vars[["PLINK_SNP_NAME", "evo2_delta_score"]]
selected_vars.rename(columns={"evo2_delta_score": "evo2_score_singlebatch"}, inplace=True)
print("Dimensions of selected_vars:", selected_vars.shape)

if input_file == "BRCA1_DATA.xlsx":
    data = pd.read_excel(os.path.join('BRCA1_DATA.xlsx'),header=2,)
    data = data[[
        'chromosome', 'position (hg19)', 'reference', 'alt', 'function.score.mean', 'func.class',
    ]]
    data.head(3)
    data.rename(columns={
        'chromosome': 'chrom',
        'position (hg19)': 'pos',
        'reference': 'ref',
        'alt': 'alt',
        'function.score.mean': 'score',
        'func.class': 'class',
    }, inplace=True)
    data['class'] = data['class'].replace(['FUNC', 'INT'], 'FUNC/INT')
    data.head(10)
else:
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input {input_file} does not exist.")
    # Load file type by extension
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == ".txt":
        print(f"Detected .txt file: {input_file}")
        data = pd.read_csv(input_file, sep="\t")
    elif file_extension == ".xlsx":
        print(f"Detected .xlsx file: {input_file}")
        data = pd.read_excel(input_file)
    else:
        raise ValueError(f"Unsupported file format: {file_extension}. Only .txt and .xlsx are supported.")
    # Split the `PLINK_SNP_NAME` column into 4 new columns: chrom, pos, ref, alt
    data[["chrom", "pos", "ref", "alt"]] = data["PLINK_SNP_NAME"].str.split(":", expand=True)
    data["pos"] = data["pos"].astype(int)
    # Rename the column LOF_DISRUPTIVE to class
    data.rename(columns={"LOF_DISRUPTIVE": "class"}, inplace=True)
    data["class"] = data["class"].replace({0: "FUNC/INT", 1: "LOF"})

data = pd.merge(selected_vars, data, on="PLINK_SNP_NAME", how="left")
print("Dimensions of data:", data.shape)
print("")
##############################################################################################
# Wrap your entire code block in a for loop that runs 20 times:
# ensure the random selection of 50 sequences per iteration is done without replacement, so that all rows are selected exactly once 

# Shuffle data to randomize row order
data = data.sample(frac=1, random_state=42).reset_index(drop=True)

# Split data into 20 equal subsets
subsets = np.array_split(data, 20)

# Loop through each subset
for counter, subset in enumerate(subsets, start=1):
    print("")
    print(f"Starting iteration {counter}...")

    # Use the current subset as the input data
    data = subset
    print("Dimensions of data:", data.shape)

    ############################################################
    # 3. Write FASTA files for ref and variant sequences
    ############################################################
    data = generate_fasta_files(data, refseq, output_dir=OUTPUT_DIR, window_size=WINDOW_SIZE)

    ############################################################
    # Load Evo2 Checkpoints
    ############################################################
    if MODEL_SIZE == "1b":
        from bionemo.core.data.load import load
        checkpoint_path = load("evo2/1b-8k-bf16:1.0")
    else:
        checkpoint_path = Path(f"nemo2_evo2_{MODEL_SIZE}_8k")
        if not checkpoint_path.exists() or not any(checkpoint_path.iterdir()):
            subprocess.run(
                f"evo2_convert_to_nemo2 --model-path hf://arcinstitute/savanna_evo2_{MODEL_SIZE}_base "
                f"--model-size {MODEL_SIZE} --output-dir nemo2_evo2_{MODEL_SIZE}_8k",
                shell=True,
                check=True,
            )
        else:
            print("Model directory is not empty. Skipping command.")

    ############################################################
    # Score sequences - ref and variant sequences
    ############################################################

    start_time = time.time()
    # Define output directories for prediction results
    output_dir = Path("fasta_files")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save reference and variant sequences to FASTA
    ref_fasta_path = output_dir / "reference_sequences.fasta"
    var_fasta_path = output_dir / "variant_sequences.fasta"
    predict_ref_dir = output_dir / "reference_predictions"
    predict_var_dir = output_dir / "variant_predictions"
    predict_ref_dir.mkdir(parents=True, exist_ok=True)
    predict_var_dir.mkdir(parents=True, exist_ok=True)

    fp8_supported, gpu_info = check_fp8_support()
    print(f"FP8 Support: {fp8_supported}")
    print(gpu_info)

    # Note: If FP8 is not supported, you may want to disable it in the model config
    # The Evo2 config has 'use_fp8_input_projections: True' by default
    if FAST_CI_MODE:
        model_subset_option = "--num-layers 4 --hybrid-override-pattern SDH*"
    else:
        model_subset_option = ""
    fp8_option = "--fp8" if fp8_supported else ""

    # Update predict commands to run on the full dataset
    predict_ref_command = (
        f"predict_evo2 --fasta {ref_fasta_path} --ckpt-dir {checkpoint_path} "
        f"--output-dir {predict_ref_dir} --model-size {MODEL_SIZE} --tensor-parallel-size 1  {model_subset_option} "
        f"--pipeline-model-parallel-size 1 --context-parallel-size 1 --output-log-prob-seqs {fp8_option}"
    )
    predict_var_command = (
        f"predict_evo2 --fasta {var_fasta_path} --ckpt-dir {checkpoint_path} "
        f"--output-dir {predict_var_dir} --model-size {MODEL_SIZE} --tensor-parallel-size 1 {model_subset_option} "
        f"--pipeline-model-parallel-size 1 --context-parallel-size 1 --output-log-prob-seqs {fp8_option}"
    )
    end_time = time.time()
    t2 = end_time - start_time
    print(f"Preparation: {t2} s")  

    # Score reference
    start_time = time.time()
    print(f"Scoring reference seq...")
    subprocess.run(predict_ref_command, shell=True, check=True)
    end_time = time.time()
    t3 = end_time - start_time

    # Score variant
    start_time = time.time()
    print(f"Scoring variant seq...")
    subprocess.run(predict_var_command, shell=True, check=True)
    end_time = time.time()
    t4 = end_time - start_time

    ############################################################
    # calculate delta likelihood scores
    ############################################################

    # Find and load prediction files
    ref_pred_files = glob.glob(os.path.join(predict_ref_dir, "predictions__rank_*.pt"))
    var_pred_files = glob.glob(os.path.join(predict_var_dir, "predictions__rank_*.pt"))

    # Load sequence ID maps (maps sequence ID -> prediction index)
    with open(os.path.join(predict_ref_dir, "seq_idx_map.json"), "r") as f:
        ref_seq_idx_map = json.load(f)
    with open(os.path.join(predict_var_dir, "seq_idx_map.json"), "r") as f:
        var_seq_idx_map = json.load(f)

    # Load predictions
    ref_preds = torch.load(ref_pred_files[0])
    var_preds = torch.load(var_pred_files[0])

    ref_log_probs = []
    var_log_probs = []

    for _, row in data.iterrows():
        ref_name = row["ref_fasta_name"]
        var_name = row["var_fasta_name"]
        ref_log_probs.append(ref_preds["log_probs_seqs"][ref_seq_idx_map[ref_name]].item())
        var_log_probs.append(var_preds["log_probs_seqs"][var_seq_idx_map[var_name]].item())
    data["ref_log_probs"] = ref_log_probs
    data["var_log_probs"] = var_log_probs
    data["evo2_delta_score"] = data["var_log_probs"] - data["ref_log_probs"]
    data["run"] = counter  # Add run number to the dataframe
    data.head()

    # Remove columns ref_fasta_name and var_fasta_name
    data.drop(columns=["ref_fasta_name", "var_fasta_name", "ref_log_probs", "var_log_probs"], inplace=True)
    data.drop(columns=["chrom", "pos", "ref", "alt"], inplace=True)

    # Write or append data to the CSV file
    csv_file = f"{input_name}_win{WINDOW_SIZE}seq{SEQ_LENGTH}{subset_method}.csv"
    if not os.path.exists(csv_file):
        print(f"Creating new file: {csv_file}")
        data.to_csv(csv_file, index=False, mode="w", header=True)
    else:
        print(f"Appending to existing file: {csv_file}")
        data.to_csv(csv_file, index=False, mode="a", header=False)

    print(f"Finished iteration {counter}.\n")

print("All 20 iterations completed!")
