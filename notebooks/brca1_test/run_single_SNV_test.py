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
warnings.simplefilter("ignore", FutureWarning)

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import torch
from Bio import SeqIO
from sklearn.metrics import auc, roc_auc_score, roc_curve
FAST_CI_MODE: bool = os.environ.get("FAST_CI_MODE", False)

######################### ######################## ########################
# Load arguments
import argparse
parser = argparse.ArgumentParser(description="Evo2 pilot study")
parser.add_argument("--SEQ_LENGTH", type=int, required=True, default=100, help="SEQ_LENGTH")
parser.add_argument("--WINDOW_SIZE", type=int, required=True, default=4096, help="Window size")
parser.add_argument("--MODEL_SIZE", type=str, help="MODEL_SIZE")
parser.add_argument("--INPUT_FILE", type=str, help="input_file")
parser.add_argument("--REF_CHR", type=str, help="REF_CHR")
parser.add_argument("--SUBSET_METHOD", type=str, help="SUBSET_METHOD")

args = parser.parse_args()
OUTPUT_DIR = "fasta_files"
SEQ_LENGTH = args.SEQ_LENGTH
WINDOW_SIZE = args.WINDOW_SIZE
MODEL_SIZE = args.MODEL_SIZE
input_file = args.INPUT_FILE
chr = args.REF_CHR
subset_method=args.SUBSET_METHOD # balanced or random; bottom or top 

########################################################################
# Define functions 
########################################################################
# Initialize the metrics files if they don't already exist
def initialize_metrics_file():
    metrics_xlsx = Path("metrics.xlsx")
    metrics_csv = Path("metrics.csv")
    if not metrics_xlsx.exists() or not metrics_csv.exists():
        columns = [
            "MODEL_SIZE", "SEQ_LENGTH", "WINDOW_SIZE", 
            "load_model", "model_prep", "score_ref", "score_var", "delta", "AUROC"
        ]
        df = pd.DataFrame(columns=columns)
        
        # Write to both .xlsx and .csv files
        if not metrics_xlsx.exists():
            df.to_excel(metrics_xlsx, index=False)
        if not metrics_csv.exists():
            df.to_csv(metrics_csv, index=False)
    return metrics_xlsx, metrics_csv

# Append results to both .xlsx and .csv files
def append_metrics_row(metrics_files, row_data):
    metrics_xlsx, metrics_csv = metrics_files
    # Update the .xlsx file
    df_xlsx = pd.read_excel(metrics_xlsx)
    df_xlsx = pd.concat([df_xlsx, pd.DataFrame([row_data])], ignore_index=True)
    df_xlsx.to_excel(metrics_xlsx, index=False)
    # Update the .csv file
    df_csv = pd.read_csv(metrics_csv)
    df_csv = pd.concat([df_csv, pd.DataFrame([row_data])], ignore_index=True)
    df_csv.to_csv(metrics_csv, index=False)

def sample_data(df, sample_frac=1.0, balanced=True, disable=True, random_state=42):
    """Sample dataframe, optionally with balanced classes.
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
    Returns:
    --------
    pandas.DataFrame
        Sampled dataframe
    """
    if disable:
        return df
    if balanced:
        # Get the number of rows in the dataframe
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
    else:
        # Calculate the number of rows to sample
        return df.sample(frac=sample_frac, random_state=random_state).reset_index(drop=True)


# Write a function that takes in seq (variable storing one an integer) and randomly subsets the dataframe (df) to SEQ_LENGTH number of rows and returns the saubset of the dataframe 
def subset_dataframe(df, seq):
    """
    Randomly subsets the dataframe to SEQ_LENGTH number of rows.
    Parameters:
    -----------
    df : pandas.DataFrame
        The input dataframe to subset.
    SEQ_LENGTH : int
        The number of rows to subset the dataframe to.
    Returns:
    --------
    pandas.DataFrame
        A subset of the dataframe with SEQ_LENGTH rows.
    """
    print("Number of rows to extract:", seq) 
    if seq > len(df):
        raise ValueError(f"SEQ_LENGTH ({seq}) is greater than the number of rows in the DataFrame ({len(df)}).")
    subset_df = df.sample(n=seq, random_state=42)
    print("New subset:", subset_df.shape) 
    return subset_df

def subset_top_bottom(data, subset_method, SEQ_LENGTH):
    """
    Extracts a subset of rows from the data based on the subset_method and SEQ_LENGTH.
    
    Parameters:
    - data (pd.DataFrame): The input data table with a 'yhat' column.
    - subset_method (str): Either "top" or "bottom".
        - "top": Extract the top SEQ_LENGTH rows with the smallest yhat values.
        - "bottom": Extract the bottom SEQ_LENGTH rows with the largest yhat values.
    - SEQ_LENGTH (int): The number of rows to extract.

    Returns:
    - pd.DataFrame: The subsetted data.
    """
    # Sort data by the 'yhat' column in ascending order
    sorted_data = data.sort_values(by="yhat", ascending=True)
    
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
    Parameters:
    -----------
    df : pandas.DataFrame
        The input dataframe from which to extract the row.
    row_num : int
        The row number to extract (1-indexed).
    Returns:
    --------
    pandas.DataFrame
        A dataframe containing only the specified row.
    Raises:
    -------
    ValueError:
        If `row_num` is less than 1 or greater than the number of rows in the dataframe.
    """
    import pandas as pd
    if row_num < 1 or row_num > len(df):
        raise ValueError(f"ROW_NUM ({row_num}) is out of bounds for the dataframe with {len(df)} rows.")
    # Extract the specific row (convert 1-indexed to 0-indexed)
    extracted_row = df.iloc[[row_num - 1]]
    return extracted_row

def parse_sequences(pos, ref, alt, refseq, window_size=WINDOW_SIZE):
    """Parse reference and variant sequences from the reference genome sequence.
    Parameters:
    -----------
    pos : int
        Position (1-indexed)
    ref : str
        Reference base
    alt : str
        Alternate base
    seq_chr17 : str
        Full chromosome 17 sequence
    window_size : int
        Size of the sequence window to extract
    Returns:
    --------
    tuple
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
    Parameters:
    -----------
    df : pandas.DataFrame
        Dataframe with variant information
    refseq : str
        Chromosome reference sequence
    output_dir : str
        Output directory for FASTA files
    window_size : int
        Size of sequence window
    Returns:
    --------
    pandas.DataFrame
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
print("Target file name:", ref_file)

# Read the reference genome sequence 
if not os.path.exists(ref_file):
    print(f"File {ref_file} does not exist")
else:
    print(f"File {ref_file} exists. Opening it...")
try:
    with open(ref_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            refseq_raw = str(record.seq)  # Convert sequence to string
            break 
except FileNotFoundError:
    print(f"Error: File {ref_file} not found in directory {os.getcwd()}.")
except Exception as e:
    print(f"An error occurred while reading the file: {e}")

refseq=refseq_raw
print(len(refseq))  # 83197441

################################################
# 2. Load SNV dataset from local environment
################################################

if not os.path.exists(input_file):
    raise FileNotFoundError(f"Input {input_file} does not exist.")

data = pd.read_csv(input_file, sep="\t")  # Assuming the file is tab-delimited

# Split the `PLINK_SNP_NAME` column into 4 new columns: chrom, pos, ref, alt
data[["chrom", "pos", "ref", "alt"]] = data["PLINK_SNP_NAME"].str.split(":", expand=True)
data["pos"] = data["pos"].astype(int)

# Rename the column LOF_DISRUPTIVE to class
data.rename(columns={"LOF_DISRUPTIVE": "class"}, inplace=True)
data["class"] = data["class"].replace({0: "FUNC/INT", 1: "LOF"})

data.head(5)
print("Dimensions of data:", data.shape)

if subset_method == "balanced":
    brca1_df = sample_data(
    brca1_df,
    sample_frac=SAMPLE_CONFIG["sample_frac"],
    balanced=SAMPLE_CONFIG["balanced"],
    disable=SAMPLE_CONFIG["disable"],
    random_state=SAMPLE_CONFIG["random_state"],
) 
if subset_method == "random":
    data = subset_dataframe(data,SEQ_LENGTH)
if subset_method == "top":
    data = subset_top_bottom(data, "top", SEQ_LENGTH)
if subset_method == "bottom":
    data = subset_top_bottom(data, "bottom", SEQ_LENGTH)
if subset_method == "row":
    row_num = 10  # row number to extract (i.e. 10th row)
    data = extract_row_by_number(data, row_num)

############################################################
# 3. Write FASTA files for ref and variant sequences
############################################################
data = generate_fasta_files(data, refseq, output_dir=OUTPUT_DIR, window_size=WINDOW_SIZE)
print("OUTPUT_DIR:", OUTPUT_DIR)

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
        print("Checkpoint directory is not empty. Skipping command.")

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
print(f"Scoring Reference seq using WINDOW_SIZE = {WINDOW_SIZE} and SEQ_LENGTH = {SEQ_LENGTH}: {t3} s")  

# Score variant
start_time = time.time()
print(f"Scoring variant seq...")
subprocess.run(predict_var_command, shell=True, check=True)
end_time = time.time()
t4 = end_time - start_time
print(f"Scoring Variant seq using WINDOW_SIZE = {WINDOW_SIZE} and SEQ_LENGTH = {SEQ_LENGTH}: {t4} s")  

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

# ideally probability of a broken variant is lower than a good one. So a bad var - good ref is negative.
data["evo2_delta_score"] = data["var_log_probs"] - data["ref_log_probs"]
data.head()
print("Saving data to current directory:", data.shape)

excel_file=f"{input_file}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}.xlsx"
data.to_excel(excel_file, index=False)

end_time = time.time()
t5 = end_time - start_time
print(f"Delta change: {t5} s") # 0.007634162902832031 seconds


####################################################################################
# Plot
####################################################################################
def plot_strip_with_means(df, x_col="evo2_delta_score", class_col="class"):
    """Creates a strip plot with jittered points and median indicators for each class using Seaborn.
    Parameters:
    - df (pd.DataFrame): The input DataFrame containing data.
    - x_col (str): The column name representing the x-axis values (e.g., evo2_delta_score).
    - class_col (str): The column name representing the class labels.
    Returns:
    - matplotlib Figure: Strip plot with median indicators.
    """
    # NVIDIA theme colors
    NVIDIA_GREEN = "#76B900"
    BACKGROUND_COLOR = "#F8F8F8"
    GRID_COLOR = "#DDDDDD"
    FONT_COLOR = "#333333"

    # Determine order of classes (if not already specified)
    unique_classes = sorted(df[class_col].unique())

    # Set up the plot with NVIDIA theme
    plt.figure(figsize=(9, 5), facecolor=BACKGROUND_COLOR)
    plt.style.use("default")  # Reset to default to avoid any pre-existing style

    # Create strip plot
    p = sns.stripplot(
        data=df,
        x=x_col,
        y=class_col,
        hue=class_col,
        order=unique_classes,
        palette=[NVIDIA_GREEN, "red"],
        size=6,
        jitter=0.3,
        alpha=0.6,
    )

    # Add median indicators using boxplot
    sns.boxplot(
        showmeans=True,
        meanline=True,
        meanprops={"visible": False},
        medianprops={"color": "black", "ls": "-", "lw": 2},
        whiskerprops={"visible": False},
        zorder=10,
        x=x_col,
        y=class_col,
        data=df,
        order=unique_classes,
        showfliers=False,
        showbox=False,
        showcaps=False,
        ax=p,
    )
    # Customize plot appearance
    plt.title(
        "Distribution of Delta Likelihoods Scores\nComparing Evo 2 likelihood scores for different SNV classes",
        color=FONT_COLOR,
        fontsize=12,
        loc="left",
    )
    plt.xlabel("Delta Likelihood Score, Evo 2", color=FONT_COLOR)
    plt.ylabel("SNV Class", color=FONT_COLOR)

    # Customize grid and tick colors
    plt.grid(color=GRID_COLOR, axis="x", linestyle="--", linewidth=0.5)
    plt.tick_params(colors=FONT_COLOR)

    # Set background color
    plt.gca().set_facecolor(BACKGROUND_COLOR)
    plt.gcf().set_facecolor(BACKGROUND_COLOR)

    plt.tight_layout()
    # return plt.gcf()
    file_name = f"{input_file}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}_dot.png"
    plt.savefig(file_name, dpi=300, bbox_inches="tight")
    print(f"Dot plot saved to {file_name}")
plot_strip_with_means(data, x_col="evo2_delta_score", class_col="class")


# Calculate AUROC of zero-shot predictions
#  class 1 is LOF which is the bad thing. That means we expect this to be more negative.
y_true = data["class"] == "LOF"
auroc = roc_auc_score(y_true, -data["evo2_delta_score"])
print(f"Zero-shot prediction AUROC: {auroc:.2}")

def plot_roc_curve(df):
    """Plots an ROC curve using Seaborn with a light NVIDIA-themed design.

    The function assumes:
    - `class` column as the true labels (binary, 'LOF' = 1, else 0).
    - `evo2_delta_score` as the prediction score.
    Parameters:
    - df (pd.DataFrame): DataFrame containing `class` and `evo2_delta_score`.
    Returns:
    - matplotlib Figure: ROC Curve Visualization.
    """
    # NVIDIA theme colors
    NVIDIA_GREEN = "#76B900"
    BACKGROUND_COLOR = "#F8F8F8"
    GRID_COLOR = "#DDDDDD"
    FONT_COLOR = "#333333"

    # Validate required columns
    if "class" not in df.columns or "evo2_delta_score" not in df.columns:
        raise ValueError("DataFrame must contain 'class' and 'evo2_delta_score' columns.")

    # Convert 'class' to binary labels: Assume 'LOF' = 1, anything else = 0
    y_true = (df["class"] == "LOF").astype(int)

    # Compute ROC curve
    fpr, tpr, _ = roc_curve(y_true, -df["evo2_delta_score"])  # Negative to align with previous logic
    roc_auc = auc(fpr, tpr)

    # Set up the plot with NVIDIA theme
    plt.figure(figsize=(9, 5), facecolor=BACKGROUND_COLOR)
    plt.style.use("default")  # Reset to default to avoid any pre-existing style

    # Plot ROC curve
    plt.plot(fpr, tpr, color=NVIDIA_GREEN, lw=3, label=f"ROC curve (AUROC = {roc_auc:.2f})")

    # Plot diagonal reference line for random guessing
    plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--")

    # Customize plot appearance
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate", color=FONT_COLOR, fontsize=12)
    plt.ylabel("True Positive Rate", color=FONT_COLOR, fontsize=12)
    plt.title(
        "Zeroshot ROC Curve\nEvaluating the discriminative performance of Evo 2 predictions",
        color=FONT_COLOR,
        fontsize=16,
        loc="left",
    )

    # Customize grid and tick colors
    plt.grid(color=GRID_COLOR, linestyle="--", linewidth=0.5)
    plt.tick_params(colors=FONT_COLOR)
    # Set background color
    plt.gca().set_facecolor(BACKGROUND_COLOR)
    # Add legend
    plt.legend(loc="lower right", frameon=True, facecolor=BACKGROUND_COLOR, edgecolor=GRID_COLOR)
    file_name = f"{input_file}_win{WINDOW_SIZE}_seq{SEQ_LENGTH}_{subset_method}_ROC.png"
    plt.savefig(file_name, dpi=300, bbox_inches="tight")
    print(f"ROC curve saved to {file_name}")
plot_roc_curve(data)

############################################################
# Append results to metrics.xlsx
############################################################

# Initialize the metrics file if necessary: creates the file with the correct headers if it doesn't exist.
metrics_file = initialize_metrics_file()

# Create a new row of results
new_row = {
    "MODEL_SIZE": MODEL_SIZE,
    "SEQ_LENGTH": SEQ_LENGTH,
    "WINDOW_SIZE": WINDOW_SIZE,
    "load_model": t1,
    "model_prep": t2,
    "score_ref": t3,
    "score_var": t4,
    "delta": t5,
    "AUROC": auroc,
}
# Append the row to the metrics file
append_metrics_row(metrics_file, new_row)
print(f"Metrics appended to {metrics_file}.")

print(f"Scoring Variant seq using WINDOW_SIZE = {WINDOW_SIZE} and SEQ_LENGTH = {SEQ_LENGTH}")  
print(f"    t3: {t3} s")  
print(f"    t4: {t4} s")  
print(f"    AUROC: {auroc:.2}")

############ 6. Extract Evo2 embeddings ###############

# evo2_model = Evo2('evo2_7b')  
# sequence_wt="TGTTCCAATGAACTTTAACACATTAGAAAA"
# sequence_mut="TGTTCCAATGAACTGTAACACATTAGAAAA"

# # Tokenize wildtype sequence
# input_ids_wt = torch.tensor(
#     evo2_model.tokenizer.tokenize(sequence_wt),  # Byte-level tokenization
#     dtype=torch.int,
# ).unsqueeze(0).to('cuda:0')  # Batch dimension + GPU acceleration

# # Tokenize variant sequence
# input_ids_mut = torch.tensor(
#     evo2_model.tokenizer.tokenize(sequence_mut),
#     dtype=torch.int,
# ).unsqueeze(0).to('cuda:0')  # Batch dimension + GPU acceleration

# # # Extract embeddings from layer 28's MLP component
# layer_name = 'blocks.28.mlp.l3' # intermediate layer (part of the model's transformer structure)

# # Extract embeddings for the wildtype sequence
# outputs_wt, embeddings_wt = evo2_model(
#     input_ids_wt,
#     return_embeddings=True,
#     layer_names=[layer_name]
# )

# # Extract embeddings for the variant sequence
# outputs_mut, embeddings_mut = evo2_model(
#     input_ids_mut,
#     return_embeddings=True,
#     layer_names=[layer_name]
# )

# # Get the embedding tensors
# embedding_tensor_wt = embeddings_wt[layer_name].squeeze(0).cpu().detach()
# embedding_tensor_mut = embeddings_mut[layer_name].squeeze(0).cpu().detach()


# # Save wildtype embeddings
# np.save("embedding_wt.npy", embedding_tensor_wt.numpy())
# torch.save(embedding_tensor_wt, "embedding_wt.pt")

# # Save variant embeddings
# np.save("embedding_mut.npy", embedding_tensor_mut.numpy())
# torch.save(embedding_tensor_mut, "embedding_mut.pt")


# # Convert embeddings to NumPy arrays
# embedding_wt_np = embedding_tensor_wt.numpy()
# embedding_mut_np = embedding_tensor_mut.numpy()

# # Perform PCA on both embeddings
# pca = PCA(n_components=2)
# embedding_wt_2d = pca.fit_transform(embedding_wt_np)
# embedding_mut_2d = pca.transform(embedding_mut_np)  # Use the same PCA transformation

# # Plot both embeddings
# plt.figure(figsize=(10, 8))

# # Wildtype embeddings
# plt.scatter(embedding_wt_2d[:, 0], embedding_wt_2d[:, 1], label='Wildtype', c='blue', alpha=0.6)

# # Variant embeddings
# plt.scatter(embedding_mut_2d[:, 0], embedding_mut_2d[:, 1], label='Variant', c='red', alpha=0.6)

# plt.legend()
# plt.colorbar(label='Position in sequence')
# plt.xlabel('PCA Dimension 1')
# plt.ylabel('PCA Dimension 2')
# plt.title('Comparison of Wildtype vs Variant Embeddings')
# plt.show()



# # Tokenize input DNA sequence (handles any length ≤1M bp)
# sequence = 'ACGT' # converted into tokens suitable for input
# input_ids = torch.tensor(
#     evo2_model.tokenizer.tokenize(sequence),  # Byte-level tokenization
#     dtype=torch.int,
# ).unsqueeze(0).to('cuda:0')  # Batch dimension + GPU acceleration

# # Extract embeddings from layer 28's MLP component
# layer_name = 'blocks.28.mlp.l3' # intermediate layer (part of the model's transformer structure)
# outputs, embeddings = evo2_model(
#     input_ids, 
#     return_embeddings=True,
#     layer_names=[layer_name]  
# )

# # Embeddings shape: (batch_size=1, sequence_length=4, hidden_dim=4096)
# print('Embeddings shape: ', embeddings[layer_name].shape)

# # save as NumPy
# import numpy as np
# embedding_tensor = embeddings[layer_name].squeeze(0).cpu().detach() # Extract the embedding tensor
# np.save("embedding.npy", embedding_tensor.numpy()) # NumPy array ,binary NumPy file format.

# # Save as PyTorch
# torch.save(embedding_tensor, "embedding.pt") # binary PyTorch tensor format.
# # Load from NumPy or Tensor file
# embedding_numpy = np.load("embedding.npy")
# embedding_tensor = torch.load("embedding.pt")
# print("Head of PyTorch embedding:\n", embedding_tensor[:5])  # First 5 rows
# print("Tail of PyTorch embedding:\n", embedding_tensor[-5:])  # Last 5 rows

# # Example: Using embeddings for variant effect prediction
refseq="TGTTCCAATGAACTTTAACACATTAGAAAA"
varseq="TGTTCCAATGAACTGTAACACATTAGAAAA"
# ref_emb = evo2_model(refseq).embeddings[layer_name]
# var_emb = evo2_model(varseq).embeddings[layer_name]

# # Calculate functional impact score
# impact_score = cosine_similarity(ref_emb, var_emb)

# # Calculate functional impact score
# impact_score = cosine_similarity(wildtype_emb, mutant_emb)

# ########### 7. Visualize embeddings (plot) ###########
# import matplotlib.pyplot as plt
# from sklearn.decomposition import PCA

# # Extract the embedding tensor (batch, sequence length, embedding dim)
# embedding_tensor = embeddings[layer_name].squeeze(0).cpu().detach().numpy()

# # Perform PCA to reduce to 2 dimensions
# pca = PCA(n_components=2)
# embedding_2d = pca.fit_transform(embedding_tensor)

# # Plot the reduced embedding
# plt.figure(figsize=(8, 8))
# plt.scatter(embedding_2d[:, 0], embedding_2d[:, 1], c=range(len(embedding_2d)), cmap='viridis')
# plt.colorbar(label='Position in sequence')
# plt.xlabel('PCA Dimension 1')
# plt.ylabel('PCA Dimension 2')
# plt.title('Visualization of Evo 2 Embedding')
# plt.show()

 # Time: 44.48s | 24.92043113708496 seconds
 # Time: 44.484 s | 25.011 seconds