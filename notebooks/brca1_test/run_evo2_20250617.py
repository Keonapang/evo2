# Run Evo2 for pilot study
# H200 (1 and 2 GPU configurations, 144 GB each)
# H100 (2 GPU configuration, 80 GB each)

# WINDOW_SIZE="1024 2048 4096 8192"
# SEQ_LENGTH="10 50 84 150 250 500 1000 2000"
# MODEL_SIZE="7b"

# for ws in $WINDOW_SIZE; do
#   for sl in $SEQ_LENGTH; do
#     echo "Running: 2_evo2_study_20250617.py -SEQ_LENGTH $sl -WINDOW_SIZE $ws -MODEL_SIZE $MODEL_SIZE"
#     python 2_evo2_study_20250617.py --SEQ_LENGTH $sl --WINDOW_SIZE $ws --MODEL_SIZE $MODEL_SIZE
#   done
# done

# !pip install matplotlib pandas seaborn scikit-learn openpyxl biopython
subprocess.run("pip install matplotlib pandas seaborn scikit-learn openpyxl biopython", shell=True, check=True)
FAST_CI_MODE: bool = os.environ.get("FAST_CI_MODE", False)

import glob
import gzip
import json
import math
import os
from pathlib import Path
import time
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import torch
from Bio import SeqIO
from sklearn.metrics import auc, roc_auc_score, roc_curve

######################### ######################## ########################
# Load arguments
parser = argparse.ArgumentParser(description="Evo2 pilot study")
parser.add_argument("--SEQ_LENGTH", type=int, required=True, default=100, help="SEQ_LENGTH")
parser.add_argument("--WINDOW_SIZE", type=int, required=True, default=4096, help="Window size")
parser.add_argument("--MODEL_SIZE", type=str, help="MODEL_SIZE")
args = parser.parse_args()
OUTPUT_DIR = "brca1_fasta_files"
DATA_DIR = "brca1"

########################################################################
# Define functions 
########################################################################
# Initialize the metrics file if it doesn't already exist
def initialize_metrics_file():
    metrics_file = Path("metrics.xlsx")
    if not metrics_file.exists():
        # Create a blank DataFrame with the specified columns
        columns = [
            "MODEL_SIZE", "SEQ_LENGTH", "WINDOW_SIZE", 
            "load_model", "model_prep", "score_ref", "score_var", "delta", "AUROC"
        ]
        df = pd.DataFrame(columns=columns)
        df.to_excel(metrics_file, index=False)
    return metrics_file

# append_metrics_row() function appends a new row of results to the file after each run.
def append_metrics_row(metrics_file, row_data):
    # Load the existing Excel file
    df = pd.read_excel(metrics_file)
    # Append the new row
    df = pd.concat([df, pd.DataFrame([row_data])], ignore_index=True)
    # Save back to the Excel file
    df.to_excel(metrics_file, index=False)

def download_data(data_dir="brca1", commit_hash="3819474bee6c24938016614411f1fa025e542bbe"):
    """Download required data files if they don't exist locally.
    Parameters:
    -----------
    data_dir : str
        Directory to store downloaded files
    commit_hash : str
        GitHub commit hash for data version
    """
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    excel_path = os.path.join(data_dir, "41586_2018_461_MOESM3_ESM.xlsx")
    genome_path = os.path.join(data_dir, "GRCh37.p13_chr17.fna.gz")

    if not os.path.exists(excel_path):
        os.system(
            f"wget https://github.com/ArcInstitute/evo2/raw/{commit_hash}/notebooks/brca1/41586_2018_461_MOESM3_ESM.xlsx -O {excel_path}"
        )
    if not os.path.exists(genome_path):
        os.system(
            f"wget https://github.com/ArcInstitute/evo2/raw/{commit_hash}/notebooks/brca1/GRCh37.p13_chr17.fna.gz -O {genome_path}"
        )
    return excel_path, genome_path

def load_genome_sequence(genome_path):
    """Load genome sequence from FASTA file.
    Parameters:
    -----------
    genome_path : str
        Path to the genome FASTA file
    Returns:
    --------
    str
        Genome sequence string
    """
    with gzip.open(genome_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return str(record.seq)
    raise ValueError("Failed to parse genome sequence")

def load_brca1_data(excel_path):
    """Load and preprocess BRCA1 data from Excel file.
    Parameters:
    -----------
    excel_path : str
        Path to the Excel file
    Returns:
    --------
    pandas.DataFrame
        Processed BRCA1 dataframe
    """
    brca1_df = pd.read_excel(excel_path, header=2)
    # Select and rename columns
    brca1_df = brca1_df[
        [
            "chromosome",
            "position (hg19)",
            "reference",
            "alt",
            "function.score.mean",
            "func.class",
        ]
    ]
    brca1_df.rename(
        columns={
            "chromosome": "chrom",
            "position (hg19)": "pos",
            "reference": "ref",
            "alt": "alt",
            "function.score.mean": "score",
            "func.class": "class",
        },
        inplace=True,
    ) # Convert to two-class system
    brca1_df["class"] = brca1_df["class"].replace(["FUNC", "INT"], "FUNC/INT")
    return brca1_df

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

# Write a function that takes in SEQ_LENGTH (variable storing one an integer) and randomly subsets the dataframe (df) to SEQ_LENGTH number of rows and returns the saubset of the dataframe 
def subset_dataframe(df, SEQ_LENGTH):
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
    print("Entering subset_dataframe function...")
    print("Original dataframe:", df.shape) 
    print("Number of rows to extract:", SEQ_LENGTH) 
    if SEQ_LENGTH > len(df):
        raise ValueError(f"SEQ_LENGTH ({SEQ_LENGTH}) is greater than the number of rows in the DataFrame ({len(df)}).")
    subset_df = df.sample(n=SEQ_LENGTH, random_state=42)
    print("New subset:", subset_df.shape) 
    return subset_df

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

############################################################
# load input data
############################################################
SAMPLE_CONFIG = {"sample_frac": 0.05, "balanced": True, "disable": False, "random_state": 42}

# 1. Download the necessary data files if not present
excel_path, genome_path = download_data(DATA_DIR)
seq_chr17 = load_genome_sequence(genome_path)

# 2. Load and preprocess variant data
brca1_df = load_brca1_data(excel_path)
print("Dimensions of brca1_df:", brca1_df.shape) # Dimensions of brca1_df: (3893, 6)

# 3. Subset data using subset_dataframe()
brca1_df2 = subset_dataframe(brca1_df,SEQ_LENGTH)
brca1_df.head(2)
print("Loaded df:", brca1_df.shape)

# 4. Write FASTA files for ref and variant sequences
brca1_df = generate_fasta_files(brca1_df, seq_chr17, output_dir=OUTPUT_DIR)
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
# if MODEL_SIZE == "1b":
#     from bionemo.core.data.load import load
#     #  This line will download the checkpoint from NGC to your $HOME/.cache/bionemo directory and return the path.
#     #  To do the same from the command line, use `CHECKPOINT_PATH=$(download_bionemo_data evo2/1b-8k-bf16:1.0)`
#     checkpoint_path = load("evo2/1b-8k-bf16:1.0")
# else:
#     checkpoint_path = Path(f"nemo2_evo2_{MODEL_SIZE}_8k")
#     # Check if the directory does not exist or is empty
#     if not checkpoint_path.exists() or not any(checkpoint_path.iterdir()):
#         !evo2_convert_to_nemo2 --model-path hf://arcinstitute/savanna_evo2_{MODEL_SIZE}_base --model-size {MODEL_SIZE} --output-dir nemo2_evo2_{MODEL_SIZE}_8k
#     else:
#         print("Checkpoint directory is not empty. Skipping command.")
        
end_time = time.time()
t1 = end_time - start_time
print(f"Load model: {t1} s")  

############################################################
# Socre sequences - ref and variant sequences
############################################################

start_time = time.time()

# Define output directories for prediction results
output_dir = Path("brca1_fasta_files")
output_dir.mkdir(parents=True, exist_ok=True)

# Save reference and variant sequences to FASTA
ref_fasta_path = output_dir / "brca1_reference_sequences.fasta"
var_fasta_path = output_dir / "brca1_variant_sequences.fasta"
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
print(f"Running command: {predict_ref_command}")
subprocess.run(predict_ref_command, shell=True, check=True)
end_time = time.time()
t3 = end_time - start_time
print(f"Scoring Reference seq using WINDOW_SIZE = {WINDOW_SIZE} and SEQ_LENGTH = {SEQ_LENGTH}: {t3} s")  

# Score variant
start_time = time.time()
print(f"Running command: {predict_var_command}")
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
print("Dimensions of brca1_df:", brca1_df.shape)

for _, row in brca1_df.iterrows():
    ref_name = row["ref_fasta_name"]
    var_name = row["var_fasta_name"]
    ref_log_probs.append(ref_preds["log_probs_seqs"][ref_seq_idx_map[ref_name]].item())
    var_log_probs.append(var_preds["log_probs_seqs"][var_seq_idx_map[var_name]].item())
brca1_df["ref_log_probs"] = ref_log_probs
brca1_df["var_log_probs"] = var_log_probs

# ideally probability of a broken variant is lower than a good one. So a bad var - good ref is negative.
brca1_df["evo2_delta_score"] = brca1_df["var_log_probs"] - brca1_df["ref_log_probs"]
brca1_df.head()
print("Saving brca1_df to current directory:", brca1_df.shape)
brca1_df.to_excel("./brca1_df.xlsx", index=False)

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
        "Distribution of Delta Likelihoods Scores\nComparing Evo 2 likelihood scores for different BRCA1 SNV classes",
        color=FONT_COLOR,
        fontsize=12,
        loc="left",
    )
    plt.xlabel("Delta Likelihood Score, Evo 2", color=FONT_COLOR)
    plt.ylabel("BRCA1 SNV Class", color=FONT_COLOR)

    # Customize grid and tick colors
    plt.grid(color=GRID_COLOR, axis="x", linestyle="--", linewidth=0.5)
    plt.tick_params(colors=FONT_COLOR)

    # Set background color
    plt.gca().set_facecolor(BACKGROUND_COLOR)
    plt.gcf().set_facecolor(BACKGROUND_COLOR)

    plt.tight_layout()
    # return plt.gcf()
    file_name = f"win{WINDOW_SIZE}_seq{SEQ_LENGTH}_dot.png"
    plt.savefig(file_name, dpi=300, bbox_inches="tight")
    print(f"Dot plot saved to {file_name}")
plot_strip_with_means(brca1_df, x_col="evo2_delta_score", class_col="class")


# Calculate AUROC of zero-shot predictions
#  class 1 is LOF which is the bad thing. That means we expect this to be more negative.
y_true = brca1_df["class"] == "LOF"
auroc = roc_auc_score(y_true, -brca1_df["evo2_delta_score"])
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
    file_name = f"win{WINDOW_SIZE}_seq{SEQ_LENGTH}_ROC.png"
    plt.savefig(file_name, dpi=300, bbox_inches="tight")
    print(f"ROC curve saved to {file_name}")
plot_roc_curve(brca1_df)



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


############ 6. Extract Evo2 embeddings ###############

# # Load 7B parameter model (40B also available)
# evo2_model = Evo2('evo2_7b')  

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
# wildtype_emb = evo2_model(sequence_wt).embeddings[layer_name]
# mutant_emb = evo2_model(sequence_mut).embeddings[layer_name]

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