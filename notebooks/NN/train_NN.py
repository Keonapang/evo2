# Train a neural network for supervised classification of BRCA1 variants using Evo2 embedding

import os
import numpy as np
import pandas as pd
DIR = "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/evo2/neural_net"
input_path = os.path.join(DIR, "input_vectors.csv")
print(input_path)
#######################################################

# Step 1: create samples of input ref/var vector pairs

# Ref_vector 1
ref_vector1 = np.random.uniform(low=0.00001, high=0.01, size=8192)
# Variant_vector 1 
var_vector1 = np.random.uniform(low=0.00001, high=0.01, size=8192)

# Step 2: Create label_vector 
label_vector = np.random.choice([0, 1], size=8192, p=[0.8, 0.2])

# Step 3: Create a single dataframe by combining the three vectors
df_combined = pd.DataFrame({
    "Ref_vector": ref_vector1,
    "Variant_vector": var_vector1,
    "Label": label_vector
})

# Step 4: save to CSV 
print(df_combined.head())
print(df_combined.shape)
df_combined.to_csv(input_path, index=False)

