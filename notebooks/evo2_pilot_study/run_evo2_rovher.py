# Run Evo2 for pilot study
# H200 (1 and 2 GPU configurations, 144 GB each)
# H100 (2 GPU configuration, 80 GB each)

# INPUTS FILES: 
#   - GRCh38 human reference genome from NCBI for specific chromosome
#   - INPUT_FILE: at minimum, it contains a single column with variants in the format: PLINK_SNP_NAME:chr:pos:ref:alt
#           - can be a .txt or .xlsx file
#           - may or may not contain a header row
#           - may also contain a column named "LOF_DISRUPTIVE" or "class" with class labels, e.g. LOF, FUNC/INT     

# OUTPUT: a single .xlsx file with Evo2 delta-likelihood scores for each variant (row)

# --------------------------------------------------------------------
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
# --------------------------------------------------------------------
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
parser.add_argument("--WINDOW_SIZE", type=int, required=True, default=4096, help="Sequence window size around the variant position. ")
parser.add_argument("--MODEL_SIZE", type=str, required=True, default="7b", help="7b or 40b")
parser.add_argument("--INPUT_FILE", type=str, required=True, help="a list of variants in format chr:pos:ref:alt")
parser.add_argument("--SUBSET_METHOD", type=str, required=True, default="random", help="random, top, bottom, balanced, all") # all = no subsetting
parser.add_argument("--USE_RANDOM_SEED", type=lambda x: x.lower() == 'true', required=True, default="True", help="Set to 'True' for reproducibility, False for randomness")
parser.add_argument("--SEQ_LENGTH", type=int, default=1, help="Input variant sequence length (number of rows to extract). ")
parser.add_argument("--REF_CHR", type=str, help="reference genome chromosome")

args = parser.parse_args()
SEQ_LENGTH = args.SEQ_LENGTH # this is row_num if subset_method is row
WINDOW_SIZE = args.WINDOW_SIZE
MODEL_SIZE = args.MODEL_SIZE
input_file = args.INPUT_FILE
chr = args.REF_CHR
subset_method=args.SUBSET_METHOD # balanced or random; bottom or top; row 
USE_RANDOM_SEED = args.USE_RANDOM_SEED

# Modify input_name 
input_name = os.path.splitext(input_file)[0] 
OUTPUT_DIR = "fasta_files"
row_num=SEQ_LENGTH

print(f" ")
print(f"run_evo2_rovher.py  -input_file {input_file} -WINDOW_SIZE {WINDOW_SIZE} -SUBSET_METHOD {subset_method} -MODEL_SIZE {MODEL_SIZE} -USE_RANDOM_SEED {USE_RANDOM_SEED}\n")

# if input_name.startswith("RovHer_"):
#     print(f"Removing 'RovHer_' prefix from input_name: {input_name}")

# MODEL_SIZE="7b"
# WINDOW_SIZE=2048
# SEQ_LENGTH=10
# input_file="BRCA1_DATA.xlsx"
# chr="17"
# subset_method="balanced" # random, top, bottom, balanced
# OUTPUT_DIR = "fasta_files"

# Output file, each row (variant) has a Evo2 delta-likelihood score 
if subset_method == "row":
    excel_file = f"{input_name}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}{row_num}.xlsx"
elif subset_method == "all":
    excel_file = f"{input_name}_win{WINDOW_SIZE}_{subset_method}.xlsx"
else:
    excel_file = f"{input_name}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}.xlsx"

########################################################################
# Define functions 
########################################################################

def initialize_metrics_file():
    """
    Initializes the metrics files (.csv), adds updated columns to the metrics files.
    """
    metrics_csv = Path("metrics.csv")
    columns = [
        "MODEL_SIZE", "input_file", "subset_method", "USE_RANDOM_SEED", "SEQ_LENGTH", "WINDOW_SIZE", 
        "time_score_ref", "time_score_var", "AUROC"
    ]
    # Create new metrics files if they don't exist
    if not metrics_csv.exists():
        df = pd.DataFrame(columns=columns)
        if not metrics_csv.exists():
            df.to_csv(metrics_csv, index=False)
    return metrics_csv

def append_metrics_row(metrics_files, row_data):
    """
    Appends a row of data to the metrics files (.csv).
    Ensures the new columns are included in the metrics files.
    """
    metrics_csv = metrics_files
    # Update the .csv file
    df_csv = pd.read_csv(metrics_csv)
    df_csv = pd.concat([df_csv, pd.DataFrame([row_data])], ignore_index=True)
    df_csv.to_csv(metrics_csv, index=False)

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
            ref_name = f"ref_pos_{row['pos']}_{row['ref']}"
            ref_entries.append(f">{ref_name}\n{ref_seq}\n")
            ref_names.append(ref_name)
            ref_seq_to_name[ref_seq] = ref_name
        else:
            ref_name = ref_seq_to_name[ref_seq]
            ref_names.append(ref_name)
        if var_seq not in var_sequences:
            var_sequences.add(var_seq)
            var_name = f"var_pos_{row['pos']}_{row['ref']}to{row['alt']}"
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
# 1. Load SNV dataset from local environment
################################################

if input_file == "BRCA1_DATA.xlsx":
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
else:
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input {input_file} does not exist.\n")
    # Load file type by extension
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == ".txt":
        print(f"Detected .txt file: {input_file}\n")
        data = pd.read_csv(input_file, sep="\t", header=None)  # Load without assuming a header
        # Check if the .txt file has only one column
        if data.shape[1] == 1:
            print("Detected one-column .txt file.\n")
            # Check if the column has the header PLINK_SNP_NAME
            if data.iloc[0, 0] != "PLINK_SNP_NAME":
                data.columns = ["PLINK_SNP_NAME"]
            else:
                # If it already has a header, reload the file with header
                data = pd.read_csv(input_file, sep="\t", header=0)
    elif file_extension == ".xlsx":
        print(f"Detected .xlsx file: {input_file}\n")
        data = pd.read_excel(input_file)
    else:
        raise ValueError(f"Unsupported format: {file_extension}. Only .txt and .xlsx supported.\n")
    # Split the `PLINK_SNP_NAME` column into 4 new columns: chrom, pos, ref, alt
    if "PLINK_SNP_NAME" in data.columns:
        data[["chrom", "pos", "ref", "alt"]] = data["PLINK_SNP_NAME"].str.split(":", expand=True)
        data["pos"] = data["pos"].astype(int)
    # Rename  LOF_DISRUPTIVE to class, only if it exists
    if "LOF_DISRUPTIVE" in data.columns:
        print("Renaming 'LOF_DISRUPTIVE' to 'class'.")
        data.rename(columns={"LOF_DISRUPTIVE": "class"}, inplace=True)
    else:
        print("'LOF_DISRUPTIVE' column not found. Skipping renaming.\n")

print("Dimensions of data:", data.shape, "\n")

################################################
# 2. Load reference file
################################################

# Extract the value in the 'chrom' column from the first row
# if data['chrom'].nunique() == 1:  # Check if there's only one unique value in 'chrom'
chr = data.loc[0, 'chrom']
print("Chromosome selected:", chr)
ref_file = f"GRCh38_chr{chr}.fasta"

if input_file == "BRCA1_DATA.xlsx":
    ref_file = "GRCh37_chr17.fasta"
refseq = str(next(SeqIO.parse(open(ref_file, "rt"), "fasta")).seq)
print(f"ref_file: {ref_file} with length {len(refseq)}\n")

################################################
# 3. Filter and subset data
################################################

# remove any rows that are exact duplicates
data = data.drop_duplicates()

if subset_method == "balanced":
    data = sample_data(data, SEQ_LENGTH=SEQ_LENGTH, balanced=True, disable=False, random_state=42, use_random_state=USE_RANDOM_SEED)
if subset_method == "random":
    data = subset_dataframe(data, SEQ_LENGTH, use_random_seed=USE_RANDOM_SEED, random_seed=42)
if subset_method == "top":
    data = subset_top_bottom(data, "top", SEQ_LENGTH)
if subset_method == "bottom":
    data = subset_top_bottom(data, "bottom", SEQ_LENGTH)
if subset_method == "row":
    data = extract_row_by_number(data, row_num=SEQ_LENGTH)
if subset_method == "all":
   print("Using all data without sampling.\n") # Sequence length is not used in this case

if "class" in data.columns:
    lof_count = data[data["class"] == "LOF"].shape[0]
    print(f"Number of rows where 'class' is 'LOF': {lof_count}\n")
    other_count = data[data["class"] == "FUNC/INT"].shape[0]
    print(f"Number of rows where 'class' is 'FUNC/INT': {other_count}\n")

############################################################
# 3. Write FASTA files for ref and variant sequences
############################################################
data = generate_fasta_files(data, refseq, output_dir=OUTPUT_DIR, window_size=WINDOW_SIZE)

############################################################
# Load Evo2 Checkpoints
############################################################
start_time = time.time()
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
        print("Checkpoint directory is not empty. Skipping command.\n")
end_time = time.time()
t1 = end_time - start_time
print(f"Load model: {t1} s")  

############################################################
# Socre sequences - ref and variant sequences
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
print(f"Scoring Reference seq: {t3} s")  

# Score variant
start_time = time.time()
print(f"Scoring variant seq...")
subprocess.run(predict_var_command, shell=True, check=True)
end_time = time.time()
t4 = end_time - start_time
print(f"Scoring Variant seq: {t4} s\n")  

############################################################
# calculate delta likelihood scores
############################################################

start_time = time.time()
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

# Remove columns ref_fasta_name and var_fasta_name
data.drop(columns=["ref_fasta_name", "var_fasta_name", "ref_log_probs", "var_log_probs"], inplace=True)
data.drop(columns=["chrom", "pos", "ref", "alt"], inplace=True)
data.head()
end_time = time.time()
t5 = end_time - start_time
print(f"Delta change: {t5} s") # 0.007634162902832031 seconds


# # Calculate AUROC of zero-shot predictions
if "class" in data.columns:
    y_true = data["class"] == "LOF" # LoF should be more negative (lower Evo2 score)
    auroc = roc_auc_score(y_true, -data["evo2_delta_score"])
    print(f"Zero-shot AUROC: {auroc:.2}")
else :
    auroc="NA"

############################################################
# Append results to metrics.xlsx
############################################################

data.to_excel(excel_file, index=False)
metrics_file = initialize_metrics_file()
new_row = {
    "MODEL_SIZE": MODEL_SIZE,
    "input_file": input_file,
    "subset_method": subset_method,
    "USE_RANDOM_SEED": USE_RANDOM_SEED,
    "SEQ_LENGTH": SEQ_LENGTH,
    "WINDOW_SIZE": WINDOW_SIZE,
    "time_score_ref": t3,
    "time_score_var": t4,
    "AUROC": auroc,
}
# Append the row to the metrics file
append_metrics_row(metrics_file, new_row)
print(f"Metrics appended to {metrics_file}.\n\n")

####################################################################################
# Plot
####################################################################################

# def plot_strip_with_means(df, x_col="evo2_delta_score", class_col="class"):
#     """Creates a strip plot with jittered points and median indicators for each class using Seaborn.
#     Parameters:
#     - df (pd.DataFrame): The input DataFrame containing data.
#     - x_col (str): The column name representing the x-axis values (e.g., evo2_delta_score).
#     - class_col (str): The column name representing the class labels.
#     Returns:
#     - matplotlib Figure: Strip plot with median indicators.
#     """
#     # NVIDIA theme colors
#     NVIDIA_GREEN = "#76B900"
#     BACKGROUND_COLOR = "#F8F8F8"
#     GRID_COLOR = "#DDDDDD"
#     FONT_COLOR = "#333333"

#     # Determine order of classes (if not already specified)
#     unique_classes = sorted(df[class_col].unique())

#     # Set up the plot with NVIDIA theme
#     plt.figure(figsize=(9, 5), facecolor=BACKGROUND_COLOR)
#     plt.style.use("default")  # Reset to default to avoid any pre-existing style

#     # Create strip plot
#     p = sns.stripplot(
#         data=df,
#         x=x_col,
#         y=class_col,
#         hue=class_col,
#         order=unique_classes,
#         palette=[NVIDIA_GREEN, "red"],
#         size=6,
#         jitter=0.3,
#         alpha=0.6,
#     )

#     # Add median indicators using boxplot
#     sns.boxplot(
#         showmeans=True,
#         meanline=True,
#         meanprops={"visible": False},
#         medianprops={"color": "black", "ls": "-", "lw": 2},
#         whiskerprops={"visible": False},
#         zorder=10,
#         x=x_col,
#         y=class_col,
#         data=df,
#         order=unique_classes,
#         showfliers=False,
#         showbox=False,
#         showcaps=False,
#         ax=p,
#     )
#     # Customize plot appearance
#     plt.title(
#         "Distribution of Delta Likelihoods Scores\nComparing Evo 2 likelihood scores for different SNV classes",
#         color=FONT_COLOR,
#         fontsize=12,
#         loc="left",
#     )
#     plt.xlabel("Delta Likelihood Score, Evo 2", color=FONT_COLOR)
#     plt.ylabel("SNV Class", color=FONT_COLOR)

#     # Customize grid and tick colors
#     plt.grid(color=GRID_COLOR, axis="x", linestyle="--", linewidth=0.5)
#     plt.tick_params(colors=FONT_COLOR)

#     # Set background color
#     plt.gca().set_facecolor(BACKGROUND_COLOR)
#     plt.gcf().set_facecolor(BACKGROUND_COLOR)

#     plt.tight_layout()

#     file_name = f"{input_name}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}_dot.png"
#     if subset_method == "row":
#         file_name = f"{input_name}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}{row_num}_dot.png"
#     else:
#         file_name = f"{input_name}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}_dot.png"

#     plt.savefig(file_name, dpi=300, bbox_inches="tight")
#     print(f"Dot plot saved to {file_name}")
# plot_strip_with_means(data, x_col="evo2_delta_score", class_col="class")


# def plot_roc_curve(df):
#     """Plots an ROC curve using Seaborn with a light NVIDIA-themed design.

#     The function assumes:
#     - `class` column as the true labels (binary, 'LOF' = 1, else 0).
#     - `evo2_delta_score` as the prediction score.
#     Parameters:
#     - df (pd.DataFrame): DataFrame containing `class` and `evo2_delta_score`.
#     Returns:
#     - matplotlib Figure: ROC Curve Visualization.
#     """
#     # NVIDIA theme colors
#     NVIDIA_GREEN = "#76B900"
#     BACKGROUND_COLOR = "#F8F8F8"
#     GRID_COLOR = "#DDDDDD"
#     FONT_COLOR = "#333333"

#     # Validate required columns
#     if "class" not in df.columns or "evo2_delta_score" not in df.columns:
#         raise ValueError("DataFrame must contain 'class' and 'evo2_delta_score' columns.")

#     # Convert 'class' to binary labels: Assume 'LOF' = 1, anything else = 0
#     y_true = (df["class"] == "LOF").astype(int)

#     # Compute ROC curve
#     fpr, tpr, _ = roc_curve(y_true, -df["evo2_delta_score"])  # Negative to align with previous logic
#     roc_auc = auc(fpr, tpr)

#     # Set up the plot with NVIDIA theme
#     plt.figure(figsize=(9, 5), facecolor=BACKGROUND_COLOR)
#     plt.style.use("default")  # Reset to default to avoid any pre-existing style

#     # Plot ROC curve
#     plt.plot(fpr, tpr, color=NVIDIA_GREEN, lw=3, label=f"ROC curve (AUROC = {roc_auc:.2f})")

#     # Plot diagonal reference line for random guessing
#     plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--")

#     # Customize plot appearance
#     plt.xlim([0.0, 1.0])
#     plt.ylim([0.0, 1.05])
#     plt.xlabel("False Positive Rate", color=FONT_COLOR, fontsize=12)
#     plt.ylabel("True Positive Rate", color=FONT_COLOR, fontsize=12)
#     plt.title(
#         "Zeroshot ROC Curve\nEvaluating the discriminative performance of Evo 2 predictions",
#         color=FONT_COLOR,
#         fontsize=16,
#         loc="left",
#     )

#     # Customize grid and tick colors
#     plt.grid(color=GRID_COLOR, linestyle="--", linewidth=0.5)
#     plt.tick_params(colors=FONT_COLOR)
#     # Set background color
#     plt.gca().set_facecolor(BACKGROUND_COLOR)
#     # Add legend
#     plt.legend(loc="lower right", frameon=True, facecolor=BACKGROUND_COLOR, edgecolor=GRID_COLOR)

#     if subset_method == "row":
#         file_name = f"{input_name}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}{row_num}_ROC.png"
#     else:
#         file_name = f"{input_name}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}_ROC.png"
#     plt.savefig(file_name, dpi=300, bbox_inches="tight")

# plot_roc_curve(data)
