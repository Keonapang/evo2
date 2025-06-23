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

# clean up
refseq = refseq_raw.lstrip("N")  # Remove leading 'N's
print(f"First bases of the sequence: {refseq[:20]}") # AAGCTTCTCACCCTGTTCCT
print(len(refseq))  # 83197441

# Either convert to uppercase or REMOVE repetitive sequences
refseq = refseq.upper()  # Convert all bases to uppercase
# def remove_soft_masked(seq):
#     return ''.join([base for base in seq if base.isupper()])

# remove_soft_masked(refseq)

print(refseq_raw[:1000])  # 
print(refseq[:1000])  # AAGCTTCTCACCCTGTTC | GATCATGCAGCTCTTCC
refseq = refseq_raw
print(len(refseq))  # 81,195,210 | 83197441

# Extract the first 1,000 nucleotides from refseq_orig
substring_to_check = refseq_orig[:100]

# Check if the substring exists in refseq
if substring_to_check in refseq:
    print("The first 1,000 nucleotides of refseq_orig are found in refseq.")
else:
    print("The first 1,000 nucleotides of refseq_orig are NOT found in refseq.")

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
row = brca1_df.iloc[1]
ref_seq, var_seq = parse_sequences(row['pos'], row['ref'], row['alt']) # 41276135, T. G

print(f'refseq length: {len(refseq)}') # 81,195,210 | 83197441
print(f'Reference SNV: {len(ref_seq)}') # 8192
print(f'Variant SNV: {len(var_seq)}') # 8192

print(row)
print('--')
print(f'Reference, SNV 0: ...{ref_seq[4082:4112]}...') # ...TGTTCCAATGAACTTTAACACATTAGAAAA...
print(f'Variant, SNV 0:   ...{var_seq[4082:4112]}...') # ...TGTTCCAATGAACTGTAACACATTAGAAAA...
print(f'Reference sequence for chr17: {refseq[4082:4112]}') # CACGCACCTGCTACACTCCTTCTTAGGGCT

#########################
from evo2.models import Evo2
model = Evo2('evo2_1b_base')# Load model

# Build mappings of unique reference sequences
ref_seqs = []
ref_seq_to_index = {}

# Parse sequences and store indexes
ref_seq_indexes = []
var_seqs = []

for _, row in brca1_df.iterrows():
    ref_seq, var_seq = parse_sequences(row['pos'], row['ref'], row['alt'])

    # Get or create index for reference sequence
    if ref_seq not in ref_seq_to_index:
        ref_seq_to_index[ref_seq] = len(ref_seqs)
        ref_seqs.append(ref_seq)
    
    ref_seq_indexes.append(ref_seq_to_index[ref_seq])
    var_seqs.append(var_seq)

ref_seq_indexes = np.array(ref_seq_indexes)

print(f'Scoring likelihoods of {len(ref_seqs)} reference sequences with Evo 2...')
ref_scores = model.score_sequences(ref_seqs)
print(f'Scoring likelihoods of {len(var_seqs)} variant sequences with Evo 2...')
var_scores = model.score_sequences(var_seqs)

# Subtract score of corresponding reference sequences from scores of variant sequences
delta_scores = np.array(var_scores) - np.array(ref_scores)[ref_seq_indexes]

# Add delta scores to dataframe
brca1_df[f'evo2_delta_score'] = delta_scores
brca1_df.head(10)

plt.figure(figsize=(4, 2))

# Plot stripplot of distributions
p = sns.stripplot(
    data=brca1_df,
    x='evo2_delta_score',
    y='class',
    hue='class',
    order=['FUNC/INT', 'LOF'],
    palette=['#777777', 'C3'],
    size=2,
    jitter=0.3,
)
# Mark medians from each distribution
sns.boxplot(showmeans=True,
            meanline=True,
            meanprops={'visible': False},
            medianprops={'color': 'k', 'ls': '-', 'lw': 2},
            whiskerprops={'visible': False},
            zorder=10,
            x="evo2_delta_score",
            y="class",
            data=brca1_df,
            showfliers=False,
            showbox=False,
            showcaps=False,
            ax=p)
plt.xlabel('Delta likelihood score, Evo 2')
plt.ylabel('BRCA1 SNV class')
plt.tight_layout()
plt.show()

# Calculate AUROC of zero-shot predictions
y_true = (brca1_df['class'] == 'LOF')
auroc = roc_auc_score(y_true, -brca1_df['evo2_delta_score'])

print(f'Zero-shot prediction AUROC: {auroc:.2}') # Zero-shot prediction AUROC: 0.73

############ 6. Extract Evo2 embeddings ###############

# Load 7B parameter model (40B also available)
evo2_model = Evo2('evo2_7b')  

# Tokenize input DNA sequence (handles any length ≤1M bp)
sequence = 'ACGT' # converted into tokens suitable for input
input_ids = torch.tensor(
    evo2_model.tokenizer.tokenize(sequence),  # Byte-level tokenization
    dtype=torch.int,
).unsqueeze(0).to('cuda:0')  # Batch dimension + GPU acceleration

# Extract embeddings from layer 28's MLP component
layer_name = 'blocks.28.mlp.l3' # intermediate layer (part of the model's transformer structure)
outputs, embeddings = evo2_model(
    input_ids, 
    return_embeddings=True,
    layer_names=[layer_name]  
)

# Embeddings shape: (batch_size=1, sequence_length=4, hidden_dim=4096)
print('Embeddings shape: ', embeddings[layer_name].shape)

# save as NumPy
import numpy as np
embedding_tensor = embeddings[layer_name].squeeze(0).cpu().detach() # Extract the embedding tensor
np.save("embedding.npy", embedding_tensor.numpy()) # NumPy array ,binary NumPy file format.

# Save as PyTorch
torch.save(embedding_tensor, "embedding.pt") # binary PyTorch tensor format.
# Load from NumPy or Tensor file
embedding_numpy = np.load("embedding.npy")
embedding_tensor = torch.load("embedding.pt")
print("Head of PyTorch embedding:\n", embedding_tensor[:5])  # First 5 rows
print("Tail of PyTorch embedding:\n", embedding_tensor[-5:])  # Last 5 rows

# Example: Using embeddings for variant effect prediction
wildtype_emb = evo2_model(sequence_wt).embeddings[layer_name]
mutant_emb = evo2_model(sequence_mut).embeddings[layer_name]

# Calculate functional impact score
impact_score = cosine_similarity(wildtype_emb, mutant_emb)


########### 7. Visualize embeddings (plot) ###########
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Extract the embedding tensor (batch, sequence length, embedding dim)
embedding_tensor = embeddings[layer_name].squeeze(0).cpu().detach().numpy()

# Perform PCA to reduce to 2 dimensions
pca = PCA(n_components=2)
embedding_2d = pca.fit_transform(embedding_tensor)

# Plot the reduced embedding
plt.figure(figsize=(8, 8))
plt.scatter(embedding_2d[:, 0], embedding_2d[:, 1], c=range(len(embedding_2d)), cmap='viridis')
plt.colorbar(label='Position in sequence')
plt.xlabel('PCA Dimension 1')
plt.ylabel('PCA Dimension 2')
plt.title('Visualization of Evo 2 Embedding')
plt.show()