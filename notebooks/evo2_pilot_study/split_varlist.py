# Run Evo2 for pilot study
# H200 (1 and 2 GPU configurations, 144 GB each)
# H100 (2 GPU configuration, 80 GB each)

# WINDOW_SIZE="1024 2048 4096 8192"
# SEQ_LENGTH="10 50 84 150 250 500 1000 2000"
# MODEL_SIZE="7b"
# input_file="RovHer_chr17.txt"
# REF_CHR="17"

# !pip install matplotlib pandas seaborn scikit-learn openpyxl biopython

import argparse
import pandas as pd
import os
import numpy as np
########################################################################
#  Load arguments
########################################################################
parser = argparse.ArgumentParser(description="Evo2 pilot study")
parser.add_argument("--INPUT_FILE", type=str, required=True, help="input_file")
parser.add_argument("--chunk_size", type=int, help="chunk_size (default: 3000)", default=3000)

args = parser.parse_args()
INPUT_FILE = args.INPUT_FILE
chunk_size = args.chunk_size

print(f"split_varlist.py -input_file {INPUT_FILE}")

########################################################################
# Define functions 
########################################################################

def split_data_into_chunks(input_file, chunk_size):
    """
    Splits a large table into smaller chunks and saves them as .txt files.

    Args:
        input_file (str): Path to the input file (e.g., 'chr17.txt').
        chunk_size (int): Number of rows per chunk. Default is 3000.

    Returns:
        None
    """
    # Load the data
    data = pd.read_csv(input_file, sep="\t")
    
    # Get the total number of rows
    total_rows = len(data)
    print(f"Total number of rows in the dataset: {total_rows}")
    
    # Calculate the number of chunks needed
    num_chunks = (total_rows // chunk_size) + (1 if total_rows % chunk_size != 0 else 0)
    print(f"Splitting data into {num_chunks} chunks, each with up to {chunk_size} rows.")
    
    # Loop through the data in chunks and save each chunk to a file
    base_name = os.path.splitext(input_file)[0]  # Get the base file name without extension
    for i in range(num_chunks):
        # Get the rows for the current chunk
        start_row = i * chunk_size
        end_row = min((i + 1) * chunk_size, total_rows)
        chunk = data.iloc[start_row:end_row]
        
        # Define the output file name
        output_file = f"{base_name}_{i + 1}.txt"
        
        # Save the chunk to a file
        chunk.to_csv(output_file, sep="\t", index=False)
        print(f"Saved chunk {i + 1} to {output_file}")
    
    print("All chunks have been saved.")

################################################
# 2. Load SNV dataset from local environment
################################################

if INPUT_FILE == "BRCA1_DATA.xlsx":
    data = pd.read_excel(os.path.join('BRCA1_DATA.xlsx'),header=2,)
else:
    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(f"Input {INPUT_FILE} does not exist.")
    # Load file type by extension
    file_extension = os.path.splitext(INPUT_FILE)[1].lower()
    if file_extension == ".txt":
        print(f"Detected .txt file: {INPUT_FILE}")
        data = pd.read_csv(INPUT_FILE, sep="\t")
    elif file_extension == ".xlsx":
        print(f"Detected .xlsx file: {INPUT_FILE}")
        data = pd.read_excel(INPUT_FILE)
    else:
        raise ValueError(f"Unsupported file format: {file_extension}. Only .txt and .xlsx are supported.")

print("Dimensions of INPUT_FILE:", data.shape)

################################################
# 3. filter: Remove any duplicate rows
################################################
data = data.drop_duplicates()

################################################
# 4. Split data into chunks
################################################
split_data_into_chunks(INPUT_FILE, chunk_size)


