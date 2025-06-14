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
parser.add_argument("--chr_ref", type=int, required=True, help="Chromosome reference file")
parser.add_argument("--gene", type=str, help="Gene name")
args = parser.parse_args()
##########################################################################
# # Define the local directory where files will be saved
# dir <- "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/GRChr38_ref_genome"
# if (!dir.exists(dir)) {
#   dir.create(dir, recursive = TRUE)
# }
# base_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/"

# for (chr in 1:22) {
#   # Construct the file name and URL
#   file_name <- paste0("chr", chr, ".fa.gz")
#   file_url <- paste0(base_url, file_name)
#   local_file_path <- file.path(dir, file_name)
#   message(paste("Downloading", file_name, "from", file_url, "..."))
  #   tryCatch({
#     download.file(file_url, destfile = local_file_path, mode = "wb")
#     message(paste("Saved", file_name, "to", local_file_path))
#   }, error = function(e) {
#     message(paste("Failed to download", file_name, ":", e$message))
#   })
# }

# set output path
root = "/home/ubuntu/nvidia-workbench"
os.makedirs(root, exist_ok=True)
outdir = f"{root}/output"
os.makedirs(outdir, exist_ok=True)
print(f"Output dir : {outdir}")

# Check if input files exists
if not os.path.exists(input):
    raise FileNotFoundError(f"Input variants {input} does not exist.")
file_name = f"chr{chr}.fa.gz"
if not os.path.exists(file_name):
    raise FileNotFoundError(f"Reference {file_name} does not exist.")

####### 2. Load SNV dataset from local environment ############

brca1_df = pd.read_excel(
    os.path.join('41586_2018_461_MOESM3_ESM.xlsx'),
    header=2, # Skip the first 2 rows (0-based index) and use the 3rd row as the header
)
brca1_df = brca1_df[[
    'chromosome', 'position (hg19)', 'reference', 'alt', 'function.score.mean', 'func.class',
]]
brca1_df.head(3)

# Rename columns
brca1_df.rename(columns={
    'chromosome': 'chrom',
    'position (hg19)': 'pos',
    'reference': 'ref',
    'alt': 'alt',
    'function.score.mean': 'score',
    'func.class': 'class',
}, inplace=True)

# Convert to two-class system
brca1_df['class'] = brca1_df['class'].replace(['FUNC', 'INT'], 'FUNC/INT')
brca1_df.head(10)

# Read the reference genome sequence of chromosome 17 f"chr{chr}.fa.gz"
with gzip.open(file_name, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        refseq = str(record.seq)  # Convert sequence to string
        break  # Only process the first record in the file
print(refseq[:20])

# Extracts A window of 8,192 bases surrounding the variant position.
def parse_sequences(pos, ref, alt):
    """
    Parse reference and variant sequences from the reference genome sequence.
    """
    p = pos - 1  # Convert to 0-indexed position
    full_seq = refseq
    ref_seq_start = max(0, p - win // 2) # defines a window of sequence context
    ref_seq_end = min(len(full_seq), p + win // 2)
    ref_seq = refseq[ref_seq_start:ref_seq_end] # extracts a subsequence (window) from the full chromosome 17 reference sequence 
    snv_pos_in_ref = min(win // 2, p) # 0-base position of the variant (ref or alt) within the extracted sequence window.
    var_seq = ref_seq[:snv_pos_in_ref] + alt + ref_seq[snv_pos_in_ref + 1:] # Create the Variant Sequence
    # Sanity checks
    assert len(var_seq) == len(ref_seq)
    assert ref_seq[snv_pos_in_ref] == ref # checks that the ref base in brca1_df matches the base in refseq at position pos
    assert var_seq[snv_pos_in_ref] == alt
    return ref_seq, var_seq

# Parse sequences for the first variant
row = brca1_df.iloc[1]
ref_seq, var_seq = parse_sequences(row['pos'], row['ref'], row['alt']) # 41276135, T. G

print(f'refseq length: {len(refseq)}') # 81,195,210
print(f'Reference SNV: {len(ref_seq)}') # 8192
print(f'Variant SNV: {len(var_seq)}') # 8192

print(row)
print('--')
print(f'Reference, SNV 0: ...{ref_seq[4082:4112]}...') # ...TGTTCCAATGAACTTTAACACATTAGAAAA...
print(f'Variant, SNV 0:   ...{var_seq[4082:4112]}...') # ...TGTTCCAATGAACTGTAACACATTAGAAAA...
print(f'Reference sequence for chr17: {refseq[4082:4112]}') # CACGCACCTGCTACACTCCTTCTTAGGGCT

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