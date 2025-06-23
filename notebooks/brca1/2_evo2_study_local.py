# Run Evo2 for pilot study
# -------- Variables ----------
# gene
# chromosome
# input varaint file
# output directory 
# how many variants
# window size

!pip install matplotlib pandas seaborn scikit-learn openpyxl
from Bio import SeqIO
import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import requests
import seaborn as sns
from sklearn.metrics import roc_auc_score

# Load arguments
parser = argparse.ArgumentParser(description="Evo2 pilot study")
parser.add_argument("--numvar", type=int, required=True, help="Number of variants to annotate per run")
parser.add_argument("--win", type=int, required=True, default=500, help="Window size")
parser.add_argument("--input", type=str, required=True,help="Input file with variants")
parser.add_argument("--chr", type=int, required=True, help="Chromosome reference file")
parser.add_argument("--gene", type=str, help="Gene name")
args = parser.parse_args()

##########################################################################
numvar=1000
win=8192
input="RovHer_chr17.txt"
chr = "17"
gene = "BRCA1"
##########################################################################
refdir = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/GRChr38_ref_genome"
dir = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2"
dir2 = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/data"

# set output path
outdir = f"./output"
os.makedirs(outdir, exist_ok=True)
print(f"Output dir : {outdir}")

# input files
input_file = f"{dir2}/{input}"
ref_file = f"{refdir}/GRCh38_chr{chr}.fasta"

################################################
# Load SNV dataset from local environment
################################################

# Check if file exists
if not os.path.exists(input_file):
    raise FileNotFoundError(f"Input {input_file} does not exist.")

data = pd.read_csv(input_file, sep="\t")  # Assuming the file is tab-delimited

# Split the `PLINK_SNP_NAME` column into 4 new columns: chrom, pos, ref, alt
data[["chrom", "pos", "ref", "alt"]] = data["PLINK_SNP_NAME"].str.split(":", expand=True)

# Convert `pos` to integer for correct data type
data["pos"] = data["pos"].astype(int)

# Print the dimensions of the resulting DataFrame
print("Final data dimensions:", data.shape)
data.head(10)

################################################
# REFERENCE GENOME
################################################

os.chdir(refdir)
print("Target file name:", ref_file)

# Read the reference genome sequence 
if not os.path.exists(ref_file):
    print(f"File {ref_file} does not exist in {refdir}.")
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

print(len(refseq_raw))  # 81195210 | 83257441

# refseq = refseq_raw.lstrip("N") # clean up: Remove leading 'N's
refseq=refseq_raw
print(f"First bases of the sequence: {refseq[:20]}") 
print(len(refseq))  # 83197441

################################################

################################################

# Extracts A window of 8,192 bases surrounding the variant position.
def parse_sequences(pos, ref, alt):
    """
    Parse reference and variant sequences from the reference genome sequence.
    """
    p = pos - 1  # Convert to 0-indexed position
    full_seq = refseq
    ref_seq_start = max(0, p - win // 2) # defines a window of sequence context
    ref_seq_end = min(len(full_seq), p + win // 2)
    ref_seq = refseq[ref_seq_start:ref_seq_end] # extracts a sequence (window) from the full reference seq 
    snv_pos_in_ref = min(win // 2, p) # 0-base position of the variant (ref or alt) within the extracted sequence window.
    var_seq = ref_seq[:snv_pos_in_ref] + alt + ref_seq[snv_pos_in_ref + 1:] # Create the Variant Sequence
    
    print(f"Reference base at pos {row['pos']} in refseq: {refseq[pos - 1]}") # T
    print(f"Expected reference base from brca1_df: {row['ref']}") # T
    print(f"ref_seq_start: {ref_seq_start}") # 41272038
    print(f"ref_seq_end: {ref_seq_end}") # 41280230
    print(f"snv_pos_in_ref: {snv_pos_in_ref}") # 4096
    print(f"len(var_seq): {len(var_seq)}") # 8192
    print(f"len(ref_seq): {len(ref_seq)}") # 8192
    # Sanity checks
    assert len(var_seq) == len(ref_seq)
    assert ref_seq[snv_pos_in_ref] == ref # checks that the ref base in brca1_df matches the base in refseq at position pos
    assert var_seq[snv_pos_in_ref] == alt
    return ref_seq, var_seq

# Parse sequences for the first variant
row = data.iloc[1]
ref_seq, var_seq = parse_sequences(row['pos'], row['ref'], row['alt']) # 41276135, T. G

print(f'refseq length: {len(refseq)}') # 81,195,210 | 83197441
print(f'Reference SNV: {len(ref_seq)}') # 8192
print(f'Variant SNV: {len(var_seq)}') # 8192

print(row)
print('--')
print(f'Reference, SNV 0: ...{ref_seq[4082:4112]}...') # ...TGTTCCAATGAACTTTAACACATTAGAAAA...
print(f'Variant, SNV 0:   ...{var_seq[4082:4112]}...') # ...TGTTCCAATGAACTGTAACACATTAGAAAA...
print(f'Reference sequence for chr17: {refseq[4082:4112]}') # CACGCACCTGCTACACTCCTTCTTAGGGCT



