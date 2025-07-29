# Train an artifical neural network for supervised classification of BRCA1 variants using Evo2 embedding
# July 25, 2025
# Purpose: 
        # 1. Can we replicate the evo2 analysis, achieving the same AUC
        # 2. Compare bewteen delta and refvar embedding concatenation methods
        # 3. Compare between using clinvar labels and class labels (LOF vs FUNC/INT)
##############################################################################################

# REGION="RovHer_BRCA1"             # BRCA1_DATA, RovHer_BRCA1 or RovHer_LDLR, "both" (BRCA1 + LDLR RVs)
# y_labels=("class" "clinvar") # "clinvar" "class" "height_FDR"
# COMBOS=("delta")                  # delta, refvar
# LAYER="blocks.28.mlp.l3"          # fixed

# for i in {1..10}; do
#   for y_label in "${y_labels[@]}"; do
#     for COMBO in "${COMBOS[@]}"; do
#       python3.11 "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/NN/BRCA1-Simulation6.py" \
#       --REGION $REGION \
#       --LAYER $LAYER \
#       --COMBO $COMBO \
#       --Y_LABEL $y_label \
#       --i $i
#     done
#   done
# done

##############################################################################################
import sys
import argparse
import glob
import gzip
import json
import math
import os
import argparse
import numpy as np
import pandas as pd
import random
import copy
import time
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, BatchNormalization
from tensorflow.keras.layers import Dropout
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras import layers, models

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import precision_recall_curve
from keras.models import Model

from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser(description="Evo2 embeddings")
parser.add_argument("--REGION", type=str, required=True, help="BRCA1_DATA, RovHer_BRCA1 or RovHer_LDLR, both (BRCA1 + LDLR RVs)")
parser.add_argument("--LAYER", required=True,type=str, help="embedding layer")
parser.add_argument("--COMBO", required=True,type=str, help="delta, refvar")
parser.add_argument("--Y_LABEL", type=str, help="clinvar (0, 0.25, 0.5, 0.75,1); class (LOF, FUNC/INT)")
parser.add_argument("--SUBSET_METHOD", type=str, help="random, top, bottom, balanced, all")
parser.add_argument("--MODEL_SIZE", type=str, help="7B or 40B")
parser.add_argument("--i", type=int, help="number of iterations")

args = parser.parse_args()
MODEL_SIZE = args.MODEL_SIZE
SUBSET_METHOD = args.SUBSET_METHOD
REGION = args.REGION
LAYER = args.LAYER
COMBO = args.COMBO
y_label = args.Y_LABEL
i=args.i
# REGION = "BRCA1_DATA" # BRCA1_DATA, RovHer_BRCA1 or RovHer_LDLR, "both" (BRCA1 + LDLR RVs)
# LAYER="blocks.28.mlp.l3"
# COMBO="refvar" # delta, refvar
# y_label="clinvar" # clinvar (0, 0.25, 0.5, 0.75,1); class (LOF, FUNC/INT)

# Input Directories
INPUT_DIR = Path("/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/July1_embedding")
INPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input data: 
delta_file = f"{INPUT_DIR}/{REGION}_{LAYER}_delta.csv"
delta_rev_file = f"{INPUT_DIR}/{REGION}_{LAYER}_delta_rev.csv"
ref_file = f"{INPUT_DIR}/{REGION}_{LAYER}_ref.csv"
var_file = f"{INPUT_DIR}/{REGION}_{LAYER}_var.csv"
ref_rev_file = f"{INPUT_DIR}/{REGION}_{LAYER}_ref_rev.csv"
var_rev_file = f"{INPUT_DIR}/{REGION}_{LAYER}_var_rev.csv"

# Embedding input files
if REGION == "BRCA1_DATA":
    file = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data/BRCA1_DATA.xlsx" # training variants + labels
else:
    DIR="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk"
    label_file1 = f"{DIR}/RARity_monogenic_benchmark/BRCAexchange/BRCA1_clinvar_cleaned.txt" 
    label_file2 = f"{DIR}/RARity_monogenic_benchmark/LOVD_LDLR/LDLR_clinvar_curated.txt" # British heart foundation-classified variants on LOVD

if REGION == "both":
    REGION = "RovHer_BRCA1" 
    ref_file1 = f"{INPUT_DIR}/{REGION}_{LAYER}_ref.csv"
    var_file1 = f"{INPUT_DIR}/{REGION}_{LAYER}_var.csv"
    ref_rev_file1 = f"{INPUT_DIR}/{REGION}_{LAYER}_ref_rev.csv"
    var_rev_file1 = f"{INPUT_DIR}/{REGION}_{LAYER}_var_rev.csv"
    REGION = "RovHer_LDLR" 
    ref_file2 = f"{INPUT_DIR}/{REGION}_{LAYER}_ref.csv"
    var_file2 = f"{INPUT_DIR}/{REGION}_{LAYER}_var.csv"
    ref_rev_file2 = f"{INPUT_DIR}/{REGION}_{LAYER}_ref_rev.csv"
    var_rev_file2 = f"{INPUT_DIR}/{REGION}_{LAYER}_var_rev.csv"
    

# Output Directories
OUTPUT_DIR = Path(f"/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/NN/BRCA1_LDLR_{COMBO}/{y_label}")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# OUTPUT FILES
plot1 = f"{OUTPUT_DIR}/{REGION}_{LAYER}_AUC_loss.png"
if i is not None:
    plot2 = f"{OUTPUT_DIR}/{REGION}_{LAYER}_ROC_{i}.png"
else:
    plot2 = f"{OUTPUT_DIR}/{REGION}_{LAYER}_ROC.png"
output_file = f"{OUTPUT_DIR}/{REGION}_{LAYER}_AUC.txt"

#######################################################
# Define function
#######################################################
from keras.losses import binary_crossentropy

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
    data = data[['chromosome', 'position (hg19)', 'reference', 'alt', 'function.score.mean', 'func.class', 'clinvar',]]
    data.rename(columns={
            'chromosome': 'chrom','position (hg19)': 'pos',
            'reference': 'ref','alt': 'alt',
            'function.score.mean': 'score','func.class': 'class', 'clinvar': 'clinvar',
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
    # Recode `clinvar` column
    data['clinvar'] = data['clinvar'].replace(recode_map)
    unique_clinvar_values = data['clinvar'].unique()
    print("Clinvar values: ", unique_clinvar_values)
    print("")
    NROWS=data.shape[0]
    if y_label == "clinvar":
        data = data[data['clinvar'] != "NA"]
        print("After removing NA in clinvar:", data.shape)

else:
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


#######################################################
# Load embeddings
#######################################################

if COMBO == "delta":
    # Variant + reverse complement embeddings 
    delta = pd.read_csv(delta_file)
    delta_reverse = pd.read_csv(delta_rev_file)

if COMBO == "refvar":
    if REGION == "BRCA1_DATA" or REGION == "RovHer_BRCA1" or REGION == "RovHer_LDLR":
        # 1. Variant + reverse complement 
        var = pd.read_csv(var_file)
        var_reverse = pd.read_csv(var_rev_file)
        # 2. Reference + reverse complement
        ref = pd.read_csv(ref_file)
        ref_reverse = pd.read_csv(ref_rev_file)
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

#######################################################
# Subset rows
#######################################################

# Check for duplicates
data = data[~data['PLINK_SNP_NAME'].duplicated(keep='first')]

if COMBO == "delta":
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
    # Tallies
    print("Filtered labels file (data):", data.shape)
    print("delta:", delta.shape, "delta_reverse:", delta_reverse.shape)
    # Check if the number of rows match
    if not (delta.shape[0] == data.shape[0] and
            delta_reverse.shape[0] == data.shape[0]):
        raise ValueError("Number of rows in embeddings do not match number of rows in data.")

if COMBO == "refvar":
    # Step 1: Compute the strict intersection of PLINK_SNP_NAME across all dfs
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
    # Tallies
    print("Filtered labels:", data.shape)
    print("var:", var.shape, "var_rev:", var_reverse.shape, "ref:", ref.shape, "ref_rev:", ref_reverse.shape)
    # Check if the number of rows match
    if not (var.shape[0] == data.shape[0] and
            var_reverse.shape[0] == data.shape[0] and
            ref.shape[0] == data.shape[0] and
            ref_reverse.shape[0] == data.shape[0]):
        raise ValueError("Number of rows in embeddings do not match number of rows in data.")

#######################################################
# Subset columns
#######################################################

print(f"-------------------------------------\n")
# Drop the 'input_file' and 'layer' columns
if COMBO == "delta":
    delta = delta.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
    delta_reverse = delta_reverse.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
    print(f"Variant embeddings: {delta.shape}")
    print(f"Variant reverse comp. embeddings: {delta_reverse.shape}\n")

if COMBO == "refvar":
    var = var.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
    var_reverse = var_reverse.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
    ref = ref.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
    ref_reverse = ref_reverse.drop(columns=['PLINK_SNP_NAME','input_file', 'layer'])
    print(f"Variant embeddings: {var.shape}") # (631, 4096)
    print(f"Variant reverse comp. embeddings: {var_reverse.shape}") # (631, 4096)
    print(f"Reference embeddings: {ref.shape}") # (631, 4096)
    print(f"Reference reverse comp. embeddings: {ref_reverse.shape}\n") # (631, 4096)

print(f"-------------------------------------\n")

#######################################################
# Build x-feature vector by concatenation
#######################################################

if COMBO == "delta":
    # feature vector for each SNV  (631, 8192)
    feature_vec = np.hstack([
        delta.values,         # delta embeddings
        delta_reverse.values, # Reverse complement
    ])
    
if COMBO == "refvar":
    # feature vector for each SNV (3893, 16384) | 16384 features per SNV
    feature_vec = np.hstack([
        ref.values,         # Reference embeddings
        ref_reverse.values, # Reverse complement of reference
        var.values,         # Variant embeddings
        var_reverse.values  # Reverse complement of variant
    ])

print(f"Concatenated x-feature vector: {feature_vec.shape}") #  (812, 8194)

# Extract y-label vectors
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

print(f"-------------------------------------\n")

#######################################################
# Split dataset
#######################################################

internal_validation_split="no"

if internal_validation_split == "yes":
    # Split dataset into test (20%) and remaining training data (80%)
    X_train, X_test, y_train, y_test = train_test_split(
        feature_vec, train_y, test_size=0.2, random_state=42, stratify=train_y
    )
    print(f"Training set size: {X_train.shape}") 
else: 
    # MANUALLY split into test (20%) and remaining training data (80%)
    X_train_val, X_test, y_train_val, y_test = train_test_split(
        feature_vec, train_y, test_size=0.2, random_state=42, stratify=train_y
    )
    # Split remaining training data into train (80%) and validation (20%)
    X_train, X_val, y_train, y_val = train_test_split(
        X_train_val, y_train_val, test_size=0.2, random_state=42, stratify=y_train_val
    )
    X_val = X_val.astype('float32')
    y_val = y_val.astype('float32')
    X_train_val = X_train_val.astype('float32')
    y_train_val = y_train_val.astype('float32')
    print(f"Validation set size: {X_val.shape}")

X_train = X_train.astype('float32')
y_train = y_train.astype('float32')
X_test = X_test.astype('float32')
y_test = y_test.astype('float32')

print(f"Training set size: {X_train.shape}") # (403, 16384)
print(f"Test set size: {X_test.shape}\n") # (127, 16384)

#######################################################
# Train ANN
#######################################################
# Create ANN model, with output layer for binary classification
input_dim = feature_vec.shape[1]
def build_model():
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

ANN_model = build_model()
# ANN_model.summary()

start_time = time.time()
if internal_validation_split == "yes":
    # reserves 15% of the training data (X_train, y_train) for validation during training
    history = ANN_model.fit(
        X_train, y_train, 
        epochs=100, 
        batch_size=64, 
        validation_split=0.15, 
    )
    print("Internal validation split enabled during training\n")
else:
    history = ANN_model.fit(
        X_train, y_train, 
        epochs=100, 
        batch_size=64, 
        verbose=2
    )
    print("No internal validation split occured.\n")

end_time = time.time()
print("Training time: ", end_time - start_time)

#######################################################
# Evlauate on test/validation set
#######################################################

if y_label == "clinvar":
    print(f"\n------------- AUC (P/LP vs B/LB) ------------\n")
    
    # Remove VUS; which are rows where y_test == 0.5
    test_mask = y_test != 0.5
    X_test_filtered = X_test[test_mask]
    y_test_filtered = y_test[test_mask]
    val_mask = y_val != 0.5
    X_val_filtered = X_val[val_mask]
    y_val_filtered = y_val[val_mask]
    
    # Recode labels: 0.75 (LP) is recoded to 1, 0.25 (LB) is recoded to 0 
    y_test_filtered = np.where(y_test_filtered == 0.75, 1, y_test_filtered) 
    y_test_filtered = np.where(y_test_filtered == 0.25, 0, y_test_filtered)
    print(f"y_test: {y_test.shape}")
    print(f"y_test_filtered: {y_test_filtered.shape}\n")
    y_val_filtered = np.where(y_val_filtered == 0.75, 1, y_val_filtered)
    y_val_filtered = np.where(y_val_filtered == 0.25, 0, y_val_filtered)
    
    # Predict probabilities on the test and validation sets
    y_test_pred_prob = ANN_model.predict(X_test_filtered).ravel()
    y_val_pred_prob = ANN_model.predict(X_val_filtered).ravel()
    
    # Calculate AUROC for the test/validation set
    auc_test = roc_auc_score(y_test_filtered, y_test_pred_prob)
    if internal_validation_split == "no":  
        auc_val = roc_auc_score(y_val_filtered, y_val_pred_prob)
    
    # For plotting ROC curve
    fpr_test, tpr_test, thresholds_test = roc_curve(y_test_filtered, y_test_pred_prob)
    fpr_val, tpr_val, thresholds_val = roc_curve(y_val_filtered, y_val_pred_prob)

else:
    print(f"\n------------- AUC (LOF vs FUNC/INT) ------------\n")

    y_test_pred_prob = ANN_model.predict(X_test).ravel()
    y_val_pred_prob = ANN_model.predict(X_val).ravel()

    # Calculate AUROC for the test/validation set
    auc_test = roc_auc_score(y_test, y_test_pred_prob)
    if internal_validation_split == "no":  
        auc_val  = roc_auc_score(y_val, y_val_pred_prob)

    # For plotting ROC curve
    fpr_test, tpr_test, thresholds_test = roc_curve(y_test, y_test_pred_prob)
    fpr_val, tpr_val, thresholds_val = roc_curve(y_val, y_val_pred_prob)

print(f"AUC (Test): {auc_test:.4f}  (Validation): {auc_val:.4f}")
print("Results in:", f"{OUTPUT_DIR}\n")

#######################################################
# Plotting: train loss/AUC
#######################################################
# Plot training and validation loss
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(history.history['loss'], label='Training Loss')
plt.title('Training Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()

# Plot training and validation AUC
plt.subplot(1, 2, 2)
plt.plot(history.history['auc'], label='Training AUC')
plt.title('Training AUC')
plt.xlabel('Epochs')
plt.ylabel('AUC')
plt.legend()
plt.tight_layout()
plt.show()
plt.savefig(plot1)


#######################################################
# ROC
#######################################################
FONT_COLOR = "#333333"
NVIDIA_GREEN = "#76B900"
GRID_COLOR = "#DDDDDD"
BACKGROUND_COLOR = "#F8F8F8"

plt.figure(figsize=(8, 6),facecolor=BACKGROUND_COLOR)
plt.style.use("default")
# Test set
plt.plot(fpr_test, tpr_test, label=f"Test AUC = {auc_test:.3f}", color="blue", linewidth=2.5)

# Validation set
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
# Appending to results tab-separated file 
#######################################################
if i is None:
    i = 1  # Default iteration if not specified
if not os.path.exists(output_file):
    with open(output_file, 'w') as f:
        f.write("iteration\ty_label\tREGION\tCOMBO\tAUC_test\tAUC_val\n")

if os.path.exists(output_file):
    with open(output_file, 'a') as f:
        f.write(f"{i}\t{y_label}\t{REGION}\t{COMBO}\t{auc_test:.4f}\t{auc_val:.4f}\n")
