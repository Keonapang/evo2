# Train an artifical neural network for supervised classification of BRCA1 variants using Evo2 embedding
# July 24, 2025
# Purpose: 
        # 1. Can we replicate RovHer framework but using neural network? Does adding embeddings help?
        # 2. Compare bewteen delta and refvar embedding concatenation methods
        # 3. Compare between using clinvar labels and class labels (LOF vs FUNC/INT)
##############################################################################################

# EMBED_COLS=("refvar")         # delta, refvar,no
# EMBED_COLS="no"       # yes or no
# ANNO_COLS="yes"        # yes or no

# CORES=5
# LAYER="blocks.28.mlp.l3"  # fixed
# REGION="RovHer_chr17"     # BRCA1_DATA, RovHer_BRCA1 or RovHer_LDLR, RovHer_chr17
# y_labels=("height_FDR")   # "clinvar" "class" "height_FDR"
#   for y_label in "${y_labels[@]}"; do
#     for EMBED_COL in "${EMBED_COLS[@]}"; do
#       python3.11 "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/NN/rovher_NN_score.py" \
#       --REGION $REGION \
#       --LAYER $LAYER \
#       --EMBED_COLS $EMBED_COL \
#       --Y_LABEL $y_label \
#       --ANNO_COLS $ANNO_COLS \
#       --EPOCH $EPOCH \
#       --i $i
#     done
#   done


# model="NN" # "NN" or "MARS"
# REGION="BRCA1_DATA" # "RovHer_BRCA1"  or "BRCA1_DATA"
# y_layer="class" # "class" or "clinvar"
# EMBED_COL="refvar" # "delta", "refvar", or "no"
# folder_name="July25_embedding" # MODIFY
# anno="no"
# VAR_WIN="256"
# reverse="yes"
# Rscript /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/NN/plot_AUC_perlayer.r \
# $model $REGION $y_layer $EMBED_COL $folder_name $anno $VAR_WIN $reverse

##############################################################################################
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import sys
import argparse
import glob
import gzip
import json
import math
import argparse
import numpy as np
import pandas as pd
import random
import copy
import time
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, BatchNormalization
from tensorflow.keras.layers import Dropout
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras import layers, models
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import precision_recall_curve
from keras.models import Model
from Bio.Seq import Seq
from Bio import SeqIO
import multiprocessing

parser = argparse.ArgumentParser(description="Evo2 embeddings")
parser.add_argument("--REGION", type=str, required=True, help="BRCA1_DATA, RovHer_BRCA1 or RovHer_LDLR, both (BRCA1 + LDLR RVs)")
parser.add_argument("--LAYER", required=True,type=str, help="embedding layer")
parser.add_argument("--Y_LABEL", type=str, required=True, help="clinvar (0, 0.25, 0.5, 0.75,1); class (LOF, FUNC/INT)")
parser.add_argument("--EMBED_COLS", required=True, type=str, help="yes or no")
parser.add_argument("--ANNO_COLS", required=True, type=str, help="yes or no")
parser.add_argument("--MODEL_SIZE", type=str, help="7B or 40B")
parser.add_argument("--i", type=int, help="number of iterations")
parser.add_argument("--EPOCH", type=int, help="number of EPOCHs")
parser.add_argument("--VAR_WIN", type=int, help="number of EPOCHs")

args = parser.parse_args()
ANNO_COLS = args.ANNO_COLS
EMBED_COLS = args.EMBED_COLS
REGION = args.REGION
LAYER = args.LAYER
y_label = args.Y_LABEL
i=args.i
MODEL_SIZE = args.MODEL_SIZE
epoch=args.EPOCH
VAR_WIN = args.VAR_WIN

# ANNO_COLS = "yes"
# EMBED_COLS="no" # delta, refvar, no
# REGION = "RovHer_chr17" # "RovHer_chr17", BRCA1_DATA, RovHer_BRCA1 RovHer_LDLR, "both" (BRCA1 + LDLR RVs)
# LAYER="blocks.28.mlp.l3"
# y_label="height_FDR" # "height_FDR", "clinvar" (0, 0.25, 0.5, 0.75,1); "class" (LOF, FUNC/INT)
# i=1
# VAR_WIN="128"

folder_name="July25_embedding"  # MODIFY

# Input Directories
server="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk"
INPUT_DIR = Path(f"{server}/evo2/{folder_name}")
INPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input training variants + labels  
if REGION == "BRCA1_DATA":
    file = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data/BRCA1_DATA.xlsx" # training variants + labels
elif REGION == "RovHer_BRCA1" or REGION == "RovHer_LDLR":
    label_file1 = f"{server}/RARity_monogenic_benchmark/BRCAexchange/BRCA1_clinvar_cleaned.txt" 
    label_file2 = f"{server}/RARity_monogenic_benchmark/LOVD_LDLR/LDLR_clinvar_curated.txt" # British heart foundation-classified variants on LOVD
elif "FDR" in y_label:
    DIR_LABEL=f"{server}/Training/{y_label}" # old folder (no cross validation): allvar_predict_saved_models_NO_GENESETS 
    INPUT_DIR=f"{server}/evo2/July20_embedding"
    label_file=f"{DIR_LABEL}/genebass_AF001_chrall_vars_gene_final.txt"
    rovher_file=f"{DIR_LABEL}/GENEBASS_RV_train_RV_pred_Oct10/cv_MARS_d3_np15_nv0.1_lowerAF3e-05_upperAF1e-02_PLINK_AF_LOF_yhat.txt"
else:
    raise ValueError("Invalid REGION specified. Choose from: BRCA1_DATA, RovHer_BRCA1, RovHer_LDLR, or a chr region (e.g., chr1).")

# Input embedding files: 
delta_file = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_delta.csv"
delta_rev_file = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_delta_rev.csv"
ref_file = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_ref.csv"
var_file = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_var.csv"
ref_rev_file = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_ref_rev.csv"
var_rev_file = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_var_rev.csv"

# Error check
if not os.path.exists(var_file):
    print(f"\nMISSING data from layer {LAYER} and VAR_WIN {VAR_WIN}\n")
    print(f"var_file: {var_file}\n")
    print(f"\nCheck dir: {INPUT_DIR}\n")
    print(f"\n================ Terminated ===============\n")
    sys.exit(1)

# Input Embedding files
if REGION == "both":
    REGION = "RovHer_BRCA1" 
    ref_file1 = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_ref.csv"
    var_file1 = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_var.csv"
    ref_rev_file1 = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_ref_rev.csv"
    var_rev_file1 = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_var_rev.csv"
    REGION = "RovHer_LDLR" 
    ref_file2 = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_ref.csv"
    var_file2 = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_var.csv"
    ref_rev_file2 = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_ref_rev.csv"
    var_rev_file2 = f"{INPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_var_rev.csv"
    
# Output Directories
if "FDR" in y_label:
    ROOTDIR = Path(f"{server}/evo2/NN/RovHer")
else:
    ROOTDIR = Path(f"{server}/evo2/NN/{folder_name}")
OUTPUT_DIR1 = Path(f"{ROOTDIR}/{y_label}") # height_FDR, clinvar, class
OUTPUT_DIR2 = Path(f"{OUTPUT_DIR1}/{EMBED_COLS}")  # delta, refvar, no
OUTPUT_DIR = Path(f"{OUTPUT_DIR2}/anno{ANNO_COLS}") # yes or no

ROOTDIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR1.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR2.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# OUTPUT files
plot1 = f"{OUTPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_AUC_loss.png"
if i > 1:
    plot2 = f"{OUTPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_ROC_{i}.png"
else:
    plot2 = f"{OUTPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_ROC.png"
plot3 = f"{OUTPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_scatter.png"
score_file = f"{OUTPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}.txt"
README_file = f"{OUTPUT_DIR}/{REGION}_{LAYER}_VARWIN{VAR_WIN}_README.txt"

# if score_file already exists,
if os.path.exists(score_file):
    print(f"\n{score_file} already exists. Exiting...\n")
    print(f"\n================ Terminated ===============\n")
    sys.exit(1)

metrics_file = f"{ROOTDIR}/{REGION}_VARWIN{VAR_WIN}_{y_label}_ANNO{ANNO_COLS}_AUC.csv"

# Plot colors (AUROC, loss, scatter)
FONT_COLOR = "#333333"
NVIDIA_GREEN = "#76B900"
GRID_COLOR = "#DDDDDD"
BACKGROUND_COLOR = "#F8F8F8"

#######################################################
# Define function
#######################################################
from keras.losses import binary_crossentropy

# for binary classification (output layer with sigmoid activation and binary_crossentropy loss)
def build_model(input_dim):
    model = Sequential()
    model.add(Dense(512, activation='relu',input_dim=input_dim))
    model.add(BatchNormalization())
    model.add(Dropout(0.3))
    model.add(Dense(128, activation='relu'))
    model.add(BatchNormalization())
    model.add(Dropout(0.3))
    model.add(Dense(32, activation='relu'))
    model.add(BatchNormalization())
    model.add(Dense(1, activation='sigmoid'))
    model.compile(optimizer='adam',loss='binary_crossentropy',metrics=['AUC'])
    return model

# for continous data, single dense neuron with a linear activation function
def build_cont_model(input_dim):
    model = Sequential()
    model.add(Dense(512, activation='relu', input_dim=input_dim))
    model.add(BatchNormalization())
    model.add(Dropout(0.3))
    model.add(Dense(128, activation='relu'))
    model.add(BatchNormalization())
    model.add(Dropout(0.3))
    model.add(Dense(32, activation='relu'))
    model.add(BatchNormalization())
    model.add(Dense(1, activation='linear'))  # Linear activation function for cont variable regression
    model.compile(optimizer='adam', loss='mean_squared_error', metrics=['MAE'])  # Regression loss and metrics
    return model

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

# Compute Binary Cross-Entropy using NumPy
def binary_cross_entropy_np(y_true, y_pred):
    """
    Calculates Binary Cross-Entropy loss for multiple samples using NumPy.
    y_true: NumPy array of actual labels (0s and 1s)
    y_pred: NumPy array of predicted probabilities (between 0 and 1)
    """
    epsilon = 1e-15  # Small value to prevent log(0)
    y_pred = np.clip(y_pred, epsilon, 1 - epsilon)  # Clip probabilities
    loss = -np.mean(y_true * np.log(y_pred) + (1 - y_true) * np.log(1 - y_pred))
    return loss

def recode_clinvar(value):
    mapping = {
        "P": 1,
        "B": 0,
        "LB": 0.25,
        "LP": 0.75,
        "LP,P": 0.75,
        "B/LB": 0.25
    }
    return mapping.get(value, 0.5)

#######################################################
# Load training labels
#######################################################

recode_map = {
    "Pathogenic": 1,
    'Pathogenic/Likely pathogenic': 1,
    'Likely pathogenic': 0.75,
    "Uncertain significance": 0.5,
    "Likely benign": 0.25,
    "Benign": 0,
    "absent": "NA",
    'Conflicting interpretations of pathogenicity': "NA",
}

# 1. Variant data + ClinVar labels 
print(f"Loading training labels (", y_label, ")...\n")
if REGION == "BRCA1_DATA":
    data = pd.read_excel(file, header=2)
    data = data[['chromosome', 'position (hg19)', 'reference', 'alt','func.class', 'clinvar',]]
    data.rename(columns={
            'chromosome': 'chrom','position (hg19)': 'pos',
            'reference': 'ref','alt': 'alt',
            'func.class': 'class', 'clinvar': 'clinvar',
        }, inplace=True)
    # Re-code values 
    data['class'] = data['class'].replace(['FUNC', 'INT'], 'FUNC/INT')
    data['class'] = data['class'].replace({'FUNC/INT': 0, 'LOF': 1})
    # print number of values that are 1
    print("Counts of FUNC/INT (0) and LoF (1) values:")   
    print(data['class'].value_counts())
    print("")
    # Create new column 
    data['PLINK_SNP_NAME'] = data.apply(
            lambda row: f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}", axis=1
    )
    # Remove columns chrom, pos, ref and alt
    data = data.drop(columns=['chrom', 'pos', 'ref', 'alt'])
    # Recode `clinvar` column
    data['clinvar'] = data['clinvar'].replace(recode_map)
    unique_clinvar_values = data['clinvar'].unique()
    print("Clinvar values: ", unique_clinvar_values)
    print("")
    NROWS=data.shape[0]
    if y_label == "clinvar":
        data = data[data['clinvar'] != "NA"]
        print("\nAfter removing NA in clinvar:", data.shape)

if REGION == "RovHer_BRCA1" or REGION == "RovHer_LDLR":
    # BRCA1
    ACMG_col1 = pd.read_csv(label_file1, sep="\t", usecols=["PLINK_SNP_NAME", "ACMG_final"])
    ACMG_col1 = ACMG_col1.rename(columns={"ACMG_final": "clinvar"})
    gene = pd.read_csv("/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data/RovHer_BRCA1.txt", sep="\t", usecols=["PLINK_SNP_NAME", "LOF_DISRUPTIVE"])
    gene = gene.rename(columns={"LOF_DISRUPTIVE": "class"})
    ACMG_col1 = pd.merge(ACMG_col1, gene, on="PLINK_SNP_NAME", how="left")

    # LDLR
    ACMG_col2 = pd.read_csv(label_file2, sep="\t", usecols=["PLINK_SNP_NAME", "clinvar_clnsig"])
    ACMG_col2 = ACMG_col2.rename(columns={"clinvar_clnsig": "clinvar"})

    gene = pd.read_csv("/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data/RovHer_LDLR.txt", sep="\t", usecols=["PLINK_SNP_NAME", "LOF_DISRUPTIVE"])
    gene = gene.rename(columns={"LOF_DISRUPTIVE": "class"})
    ACMG_col2 = pd.merge(ACMG_col2, gene, on="PLINK_SNP_NAME", how="left")

    # Combine 
    data = pd.concat([ACMG_col1, ACMG_col2], ignore_index=True)
    print(f"BRCA1 and LDLR merged: {data.shape}")

    # If ClinVar is chosen, remove NA rows
    NROWS=data.shape[0]
    if y_label == "clinvar":
        data = data[~data["clinvar"].isin(["", "NA", "CCP"])]
        data["clinvar"] = data["clinvar"].apply(recode_clinvar)
        print(data["clinvar"].value_counts(dropna=False))

if "FDR" in y_label:
    print(f"\nLoading {y_label} and 75 functional annotations...\n")
    data = pd.read_table(label_file, sep="\s+")

    print(f"Loading RovHer scores...\n")
    score_df = pd.read_csv(rovher_file, sep="\t", usecols=["PLINK_SNP_NAME", "yhat"])
    score_df.rename(columns={"yhat": "rovher_score"}, inplace=True)
    data = pd.merge(data, score_df, on="PLINK_SNP_NAME", how="left")
    print("")
    # Subset to chr17 rows
    data = data[data['PLINK_SNP_NAME'].str.startswith("17:")]
    NROWS=data.shape[0]


# Check for duplicates
data = data[~data['PLINK_SNP_NAME'].duplicated(keep='first')]

#######################################################
# Load embeddings (+ reverse complement if available)
#######################################################

if EMBED_COLS == "delta":
        print("\nLoading delta embeddings...")
        delta = pd.read_csv(delta_file)
        if os.path.exists(delta_rev_file):
            print("Loading delta reverse complement embeddings...\n")
            delta_reverse = pd.read_csv(delta_rev_file)

if EMBED_COLS == "refvar":
        if REGION == "both":
            var1 = pd.read_csv(var_file1)
            var_reverse1 = pd.read_csv(var_rev_file1)
            ref1 = pd.read_csv(ref_file1)
            ref_reverse1 = pd.read_csv(ref_rev_file1)
            var2 = pd.read_csv(var_file2)
            var_reverse2 = pd.read_csv(var_rev_file2)
            ref2 = pd.read_csv(ref_file2)
            ref_reverse2 = pd.read_csv(ref_rev_file2)
            var = pd.concat([var1, var2], ignore_index=True)
            var_reverse = pd.concat([var_reverse1, var_reverse2], ignore_index=True)
            ref = pd.concat([ref1, ref2], ignore_index=True)
            ref_reverse = pd.concat([ref_reverse1, ref_reverse2], ignore_index=True)
        else:
            print("Loading ref + var embeddings...\n")
            chunk_size = 10000
            chunks = []
            for chunk in tqdm(pd.read_csv(var_file, chunksize=chunk_size), desc="Progress", unit="rows"):
                chunks.append(chunk)
            var = pd.concat(chunks, ignore_index=True)
            chunks = []
            for chunk in tqdm(pd.read_csv(ref_file, chunksize=chunk_size), desc="Progress", unit="rows"):
                chunks.append(chunk)
            ref = pd.concat(chunks, ignore_index=True)
            if os.path.exists(var_rev_file) and os.path.exists(ref_rev_file):
                print("\nLoading reverse complement embeddings...\n")
                var_reverse = pd.read_csv(var_rev_file)
                ref_reverse = pd.read_csv(ref_rev_file)

#######################################################
#  Align embeddings with training variant labels; subsetting rows
#######################################################

if EMBED_COLS == "delta":
        if "delta_reverse" not in locals():
            final_common_snp_names = list(
                set(data['PLINK_SNP_NAME'])
                .intersection(delta['PLINK_SNP_NAME'])
            )
            data = data[data['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            delta = delta[delta['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            if not (delta.shape[0] == data.shape[0]):
                raise ValueError("Number of rows in embeddings do not match number of rows in data.")
        else:
            # Step 1: Compute the strict intersection of PLINK_SNP_NAME across all dfs
            final_common_snp_names = list(
                set(data['PLINK_SNP_NAME'])
                .intersection(delta['PLINK_SNP_NAME'])
                .intersection(delta_reverse['PLINK_SNP_NAME'])
            )
            # Step 2: Filter all dfs simultaneously based on the common SNP names
            data = data[data['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            delta = delta[delta['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            delta_reverse = delta_reverse[delta_reverse['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            print("Filtered labels file (data):", data.shape)
            print("delta:", delta.shape, "delta_reverse:", delta_reverse.shape)
            # Check if the number of rows match
            if not (delta.shape[0] == data.shape[0] and
                    delta_reverse.shape[0] == data.shape[0]):
                raise ValueError("Number of rows in embeddings do not match number of rows in data.")

if EMBED_COLS == "refvar":
        if "var_reverse" not in locals():
            final_common_snp_names = list(
                set(data['PLINK_SNP_NAME'])
                .intersection(var['PLINK_SNP_NAME'])
                .intersection(ref['PLINK_SNP_NAME'])
            )
            data = data[data['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            var = var[var['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            ref = ref[ref['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            if not (var.shape[0] == data.shape[0] and
                    ref.shape[0] == data.shape[0]):
                raise ValueError("Number of rows in embeddings do not match number of rows in data.")
            print("var:", var.shape, "ref:", ref.shape)
        
        else:
            final_common_snp_names = list(
                set(data['PLINK_SNP_NAME'])
                .intersection(var['PLINK_SNP_NAME'])
                .intersection(var_reverse['PLINK_SNP_NAME'])
                .intersection(ref['PLINK_SNP_NAME'])
                .intersection(ref_reverse['PLINK_SNP_NAME'])
            )
            # Step 2: Filter all dfs based on the common SNP names
            data = data[data['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            var = var[var['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            var_reverse = var_reverse[var_reverse['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            ref = ref[ref['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            ref_reverse = ref_reverse[ref_reverse['PLINK_SNP_NAME'].isin(final_common_snp_names)].reset_index(drop=True)
            print("Filtered labels:", data.shape)
            print("var:", var.shape, "var_rev:", var_reverse.shape, "ref:", ref.shape, "ref_rev:", ref_reverse.shape)
            
            # 3. Check if the number of rows match
            if not (var.shape[0] == data.shape[0] and
                    var_reverse.shape[0] == data.shape[0] and
                    ref.shape[0] == data.shape[0] and
                    ref_reverse.shape[0] == data.shape[0]):
                raise ValueError("Number of rows in embeddings do not match number of rows in data.")

#######################################################
# Subset columns
#   - Drop the 'input_file' and 'layer' columns
#######################################################
plink_snp_names = data['PLINK_SNP_NAME']  # Assuming all DataFrames have the same PLINK_SNP_NAME column
N_vars = data.shape[0]

if "FDR" in y_label:
    rovher_score = data['rovher_score']  # Assuming all DataFrames have the same PLINK_SNP_NAME column

if EMBED_COLS == "delta":
        print(f"-------------------------------------\n")
        if "delta_reverse" not in locals():
            delta = delta.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
        else:
            delta = delta.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
            delta_reverse = delta_reverse.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
            print(f"Reverse complement embeddings: {delta_reverse.shape}") 
        print(f"Variant embeddings: {delta.shape}")

if EMBED_COLS == "refvar":
        if "ref_reverse" not in locals():
            var = var.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
            ref = ref.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
        else:
            var = var.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
            var_reverse = var_reverse.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
            ref = ref.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
            ref_reverse = ref_reverse.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
            print(f"Reverse complement embeddings: {var_reverse.shape}") 
        print(f"Variant embeddings: {var.shape}")
print(f"-------------------------------------\n")

#######################################################
# Clean up 75 functional annotations; removing cols 
#######################################################
if y_label == "class" or "FDR" in y_label:
    data = data.drop(columns=['clinvar'], errors='ignore')
    print(f"Data shape after dropping class: {data.shape}")

print(f"data: {data.shape}\n")
print(data.iloc[:2, :8], "\n")  # Use .iloc for positional indexing

if ANNO_COLS == "yes":
    if data.isnull().values.any():
        print("There are NA values in the data.")
        sys.exit(1)
    columns_to_remove = [
    "PLINK_SNP_NAME", "Gene.refGene", "height_FDR", "rovher_score",
    "VEST4_score", "ClinPred_score", "gMVP_score", "fathmm-XF_coding_score",
    "REVEL_score", "M-CAP_score", "MetaLR_score", "MetaRNN_score", 
    "BayesDel_noAF_score", "BayesDel_addAF_rankscore"
    ]
    data2 = data.drop(columns=columns_to_remove, errors='ignore')
    anno_vec = data2.values    

#######################################################
# Build feature vector by concatenation (embeddings)
#######################################################

if ANNO_COLS == "yes" and EMBED_COLS == "no":
    feature_vec = data2.values
    print(f"Feature_vec (annotations only): {feature_vec.shape}")
    del data2

if EMBED_COLS == "delta":
        
        if "delta_reverse" not in locals():  # Check if delta_reverse doesn't exist
            feature_vec = np.hstack([
                delta.values,  # delta embeddings only
            ])
        else:
            feature_vec = np.hstack([
                delta.values,         # delta embeddings
                delta_reverse.values, # Reverse complement
            ])
elif EMBED_COLS == "refvar":
        if "ref_reverse" not in locals():
            feature_vec = np.hstack([
                ref.values,         # Reference embeddings
                var.values,         # Variant embeddings
            ])
        else:
            feature_vec = np.hstack([
            ref.values,         # Reference embeddings
            ref_reverse.values, # Reverse complement of reference
            var.values,         # Variant embeddings
            var_reverse.values  # Reverse complement of variant
        ])
else:
        raise ValueError("Invalid EMBED_COLS specified. Choose from: delta, refvar.")
    
if ANNO_COLS == "yes":
    print(f"Concatenating {feature_vec.shape} embeddings with {anno_vec.shape} annos...\n")
    feature_vec = np.hstack([
            anno_vec,  # Embeddings
            feature_vec   # Annotations
    ])
print(f"Concatenated x-feature vector: {feature_vec.shape}") #  (812, 8194)
print(f"feature_vec:\n{feature_vec[:4, :5]}\n")  # Use [:4, :5] to show first 4 rows and 5 columns


#######################################################
# Extract y-label vectors
#######################################################
if y_label == "class":
    train_y = data["class"].values
    lof_count = data[data["class"] == 1].shape[0]
    other_count = data[data["class"] == 0].shape[0]
    print(f"Number of 'LOF': {lof_count}   'FUNC/INT': {other_count}\n")

if y_label == "clinvar":
    train_y = data["clinvar"].values
    P_count = data[data["clinvar"] == 1].shape[0]
    LP_count = data[data["clinvar"] == 0.75].shape[0]
    VUS_count = data[data["clinvar"] == 0.5].shape[0]
    LB_count = data[data["clinvar"] == 0.25].shape[0]
    B_count = data[data["clinvar"] == 0].shape[0]
    print(f"'P': {P_count}  'LP': {LP_count} 'VUS': {VUS_count} 'LB': {LB_count} 'B': {B_count} \n")

if "FDR" in y_label:
    train_y = data[y_label].values
    mean_value = np.mean(train_y)
    print(f"Mean of FDRs: {mean_value}")
    lower_value = np.min(train_y)
    upper_value = np.max(train_y)
    print(f"Lower: {lower_value}  Upper: {upper_value}\n")

del data
print(f"-------------------------------------\n")

###############################################################################
# Split dataset (if you are using class or clinvar as labels)
#   - shuffles data before splitting (default)
#   - stratify ensures that the class distribution of 
#     train_y is preserved in both the train and test sets (ONLY used when y contains categorical labels)

# Train/Validation (80%): Used for training and tuning the model.
# Test (20%): Held out for final evaluation.
###############################################################################

internal_validation_split="yes"

if y_label == "class" or y_label == "clinvar":

    if internal_validation_split == "yes":
        # Split dataset into test (20%) and remaining training set (80%)
        if "FDR" in y_label:
            X_train, X_test, y_train, y_test, snp_train, snp_test, rovher_train, rovher_test = train_test_split(
                feature_vec, train_y, plink_snp_names, rovher_score, test_size=0.2, random_state=42
            )
        else:
            X_train, X_test, y_train, y_test, snp_train, snp_test = train_test_split(
                feature_vec, train_y, plink_snp_names, test_size=0.2, random_state=42, stratify=train_y
            )
    else: # MANUALLY split into test (20%) and remaining training set (80%)
        if "FDR" in y_label:
            X_train_val, X_test, y_train_val, y_test, snp_train_val, snp_test, rovher_train_val, rovher_test_val = train_test_split(
            feature_vec, train_y, plink_snp_names, rovher_score, test_size=0.2, random_state=42
            )
            X_train, X_val, y_train, y_val, snp_train, snp_val, rovher_train, rovher_val = train_test_split(
                X_train_val, y_train_val, snp_train_val, rovher_train_val, test_size=0.2, random_state=42
            )    
        else: 
            X_train_val, X_test, y_train_val, y_test, snp_train_val, snp_test = train_test_split(
            feature_vec, train_y, plink_snp_names, test_size=0.2, random_state=42, stratify=train_y
            )
            X_train, X_val, y_train, y_val, snp_train, snp_val = train_test_split(
                X_train_val, y_train_val, snp_train_val, test_size=0.2, random_state=42, stratify=y_train_val
            )
        X_val = X_val.astype('float64')
        y_val = y_val.astype('float64')
        X_train_val = X_train_val.astype('float64')
        y_train_val = y_train_val.astype('float64')
        print(f"Validation set size: {X_val.shape}")
        print(f"snp_test: {snp_test.shape}")
        print(f"snp_train: {snp_train.shape}")

    print(f"X_train:\n{X_train[:4, :5]}\n") 

    X_train = X_train.astype('float64')
    y_train = y_train.astype('float64')
    X_test = X_test.astype('float64')
    y_test = y_test.astype('float64')

    print(f"\nTraining set size: {X_train.shape}") # (403, 16384)
    print(f"Test set size: {X_test.shape}\n") # (127, 16384)

#######################################################
# Train ANN - no cross validation
#######################################################

if y_label == "class" or y_label == "clinvar":

    start_time = time.time()
    input_dim = feature_vec.shape[1]
    
    # Create ANN model, with output layer for binary classification
    ANN_model = build_model(input_dim) # ANN_model.summary()

    if internal_validation_split == "yes":
        history = ANN_model.fit(
        X_train, y_train, 
        epochs=epoch, 
        batch_size=64, 
        validation_split=0.15, 
        )
        y_test_pred = ANN_model.predict(X_test)

        results = pd.DataFrame({
        'PLINK_SNP_NAME': snp_test.values,
        # 'rovher_score': rovher_test.flatten(),
        'rovher_NN_score': y_test_pred.flatten() 
        })
        
    else:
        history = ANN_model.fit(
        X_train, y_train, 
        epochs=epoch, 
        batch_size=64, 
        verbose=2
        )
        y_test_pred = ANN_model.predict(X_val)

        results = pd.DataFrame({
        'PLINK_SNP_NAME': snp_test.values,
        # 'rovher_score': rovher_test_val.flatten(),
        'rovher_NN_score': y_test_pred.flatten() 
        })
    print(f"Predictions shape: {y_test_pred.shape}") 
    results.to_csv(score_file, sep='\t', index=False)

    end_time = time.time()
    t = end_time - start_time
    t_minutes = t / 60
    print("Time (mins): ", t)
    print(f"\nResults: {score_file}\n")

#######################################################
# Train ANN - cross validation
#######################################################

if "chr" in REGION:
    start_time = time.time()
    input_dim = feature_vec.shape[1]
    ANN_model = build_cont_model(input_dim)
    # Initialize 5 folds
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    np.random.seed(42)
    final_predictions = pd.DataFrame(columns=["PLINK_SNP_NAME", "rovher_NN_score"])
    # final_predictions = pd.DataFrame(columns=["PLINK_SNP_NAME", "rovher_score", "rovher_NN_score"])
    fold_no = 1
    # Iterate over each fold
    for train_index, test_index in kf.split(feature_vec):
        print(f"\n======== Neural network fold {fold_no} ======== \n")
        # Split into training and test sets for the current fold
        X_train, X_test = feature_vec[train_index], feature_vec[test_index]
        y_train, y_test = train_y[train_index], train_y[test_index]
        snp_train, snp_test = np.array(plink_snp_names)[train_index], np.array(plink_snp_names)[test_index]
        # rovher_train, rovher_test = np.array(rovher_score)[train_index], np.array(rovher_score)[test_index]
        print(X_train[:5, :2])
        print(f"X_train: {X_train.shape}   y_train: {y_train.shape}")

        X_train = X_train.astype('float64')
        y_train = y_train.astype('float64')
        X_test = X_test.astype('float64')
        y_test = y_test.astype('float64')
        
        ANN_model.fit(X_train, y_train, epochs=epoch, batch_size=64, verbose=2)
        predictions = ANN_model.predict(X_test)
        fold_predictions = pd.DataFrame({
            "PLINK_SNP_NAME": snp_test,
            # "rovher_score": rovher_test.flatten(),
           "rovher_NN_score": predictions.flatten()
        })
        final_predictions = pd.concat([final_predictions, fold_predictions])
        print(f"Final predictions: {final_predictions.shape}")
        if final_predictions.isna().any().any():
            print("WARNING : Theres NA in final_predictions!")
            print(final_predictions[final_predictions.isna().any(axis=1)])
    fold_no += 1
    final_predictions.to_csv(score_file, index=False, sep='\t')
    print(f"Final: {final_predictions.shape}")
    end_time = time.time()
    t = end_time - start_time
    t_minutes = t / 60  # Convert seconds to minutes
    print("Training time (mins): ", t)
    print(f"\nResults: {score_file}\n")

#######################################################
# PLOT: Scatter plot of rovher_NN_score vs rovher_score
# (for continous)
#######################################################

if "FDR" in y_label:
    df = pd.read_csv(score_file, sep='\t')
    plt.figure(figsize=(8, 6), facecolor=BACKGROUND_COLOR)
    plt.style.use("default")
    plt.scatter(
        df["rovher_NN_score"],  # x-axis
        df["rovher_score"],     # y-axis
        color=NVIDIA_GREEN,     # Points color
        alpha=0.8,              # Transparency
        edgecolor="black",      # Border color for points
        linewidth=0.5
    )
    # Customize axes labels
    plt.xlabel("RovHer scores (NN)", color=FONT_COLOR, fontsize=12)
    plt.ylabel("RovHer scores (MARS)", color=FONT_COLOR, fontsize=12)
    plt.title("Chr17 RV RovHer scores generated by neural network vs MARS\n y-variable: height FDR", color=FONT_COLOR, fontsize=14)
    plt.grid(color=GRID_COLOR, linestyle="--", linewidth=0.7)
    plt.gca().set_facecolor(BACKGROUND_COLOR)
    plt.tick_params(colors=FONT_COLOR)
    plt.tight_layout()
    plt.savefig(plot3)
    plt.show()

#######################################################
# Evlauate on test/validation set (binary classification)
#######################################################

if y_label == "clinvar":
    print(f"\n------------- AUC (P/LP vs B/LB) ------------\n")
    
    # Remove VUS; which are rows where y_test == 0.5
    test_mask = y_test != 0.5
    X_test_filtered = X_test[test_mask]
    y_test_filtered = y_test[test_mask]
    if internal_validation_split == "no":
        val_mask = y_val != 0.5
        X_val_filtered = X_val[val_mask]
        y_val_filtered = y_val[val_mask]
    
    # Recode labels: 0.75 (LP) is recoded to 1, 0.25 (LB) is recoded to 0 
    y_test_filtered = np.where(y_test_filtered == 0.75, 1, y_test_filtered) 
    y_test_filtered = np.where(y_test_filtered == 0.25, 0, y_test_filtered)

    if internal_validation_split == "no":
        y_val_filtered = np.where(y_val_filtered == 0.75, 1, y_val_filtered)
        y_val_filtered = np.where(y_val_filtered == 0.25, 0, y_val_filtered)
        
    # Predict probabilities on the test and validation sets
    y_test_pred_prob = ANN_model.predict(X_test_filtered).ravel()
    if internal_validation_split == "no":
        y_val_pred_prob = ANN_model.predict(X_val_filtered).ravel()
    
    # Calculate AUROC for the test/validation set
    auc_test = roc_auc_score(y_test_filtered, y_test_pred_prob)
    if internal_validation_split == "no":  
        auc_val = roc_auc_score(y_val_filtered, y_val_pred_prob)
    
    # For plotting ROC curve
    fpr_test, tpr_test, thresholds_test = roc_curve(y_test_filtered, y_test_pred_prob)
    if internal_validation_split == "no":
        fpr_val, tpr_val, thresholds_val = roc_curve(y_val_filtered, y_val_pred_prob)
    print(f"AUC (Test): {auc_test:.4f}")

if y_label == "class":

    print(f"\n------------- AUC (LOF vs FUNC/INT) ------------\n")
    # Calculate AUROC for the test/validation set

    y_test_pred_prob = ANN_model.predict(X_test).ravel()
    auc_test = roc_auc_score(y_test, y_test_pred_prob)
    fpr_test, tpr_test, thresholds_test = roc_curve(y_test, y_test_pred_prob)
    
    if internal_validation_split == "no":
        y_val_pred_prob = ANN_model.predict(X_val).ravel()
        auc_val  = roc_auc_score(y_val, y_val_pred_prob)
        fpr_val, tpr_val, thresholds_val = roc_curve(y_val, y_val_pred_prob)

    print(f"AUC (Test): {auc_test:.4f}")

#######################################################
# PLOT: training loss/AUC (binary classification)
#######################################################

if y_label == "class" or y_label == "clinvar":
    plt.figure(figsize=(8, 5))

    # Plot Training Loss
    plt.plot(history.history['loss'], label='Training Loss', color='blue')

    # Plot Training AUC
    plt.plot(history.history['auc'], label='Training AUC', color='orange')

    # Add title, labels, and legend
    plt.title(f'Neural Network training on {y_label}')
    plt.xlabel('Epochs')
    plt.ylabel('Value')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig(plot1)

#######################################################
# Plot: AUROC Curve (for binary classification)
#######################################################

if y_label == "class" or y_label == "clinvar":

    plt.figure(figsize=(8, 6),facecolor=BACKGROUND_COLOR)
    plt.style.use("default")
    # Test set
    plt.plot(fpr_test, tpr_test, label=f"Test AUC = {auc_test:.3f}", color=NVIDIA_GREEN, linewidth=2.5)

    # Validation set
    if internal_validation_split == "no":  
        plt.plot(fpr_val, tpr_val, label=f"Validation AUC = {auc_val:.3f}", color=NVIDIA_GREEN, linewidth=2.5)

    # Random guess line
    plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--")

    plt.xlabel("False Positive Rate", color=FONT_COLOR, fontsize=12)
    plt.ylabel("True Positive Rate", color=FONT_COLOR, fontsize=12)
    plt.title(f"ROC Curve for {y_label}", color=FONT_COLOR, fontsize=16)
    plt.legend(loc="lower right", fontsize=10)
    plt.legend(loc="lower right", frameon=True, fontsize=10, facecolor=BACKGROUND_COLOR, edgecolor=GRID_COLOR)

    # Grid and layout adjustments
    plt.grid(alpha=0.3, color=GRID_COLOR, linestyle="--", linewidth=0.5)
    plt.tick_params(colors=FONT_COLOR)
    plt.tight_layout()
    plt.savefig(plot2)
    plt.show()

#######################################################
# Append AUROC to file  (binary classification)
#######################################################

if y_label == "class" or y_label == "clinvar":
    if not os.path.exists(metrics_file):
        with open(metrics_file, 'w') as f:
            f.write("REGION,y_label,LAYER,VAR_WIN,EMBED_COLS,ANNO_COLS,i,N_vars,AUC_test,Train_mins\n")

    # Append data 
    if os.path.exists(metrics_file):
        with open(metrics_file, 'a') as f:
            f.write(f"{REGION},{y_label},{LAYER},{VAR_WIN},{EMBED_COLS},{ANNO_COLS},{i},{N_vars},{auc_test:.4f},{t_minutes:.2f}\n")

print(f"\n --------- Neural network with {y_label} labels ---------\n")
print(f"RESULTS : {OUTPUT_DIR}\n")