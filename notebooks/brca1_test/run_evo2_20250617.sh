#!/bin/bash
# To run Evo2 pilot study with different parameters
# M:\pangk\Softwares\evo2\notebooks\brca1_test\run_evo2_20250617.py

pip install matplotlib pandas seaborn scikit-learn openpyxl biopython
WINDOW_SIZE="1024"
SEQ_LENGTH="3890"
MODEL_SIZE="7b" # 1.5 minutes to download

WINDOW_SIZE="1024"
SEQ_LENGTH="10"
MODEL_SIZE="40b"

WINDOW_SIZE="1024 2048 4096 8192"
SEQ_LENGTH="50 84 150 250 500 1000 2000 8000"
MODEL_SIZE="40b" # 33 minutes to download, model split into 4 parts (savanna_evo2_40b_base.pt.part0-3) stored in ./nemo2_evo2_40b_8k

for ws in $WINDOW_SIZE; do
  for sl in $SEQ_LENGTH; do
    echo "Running: run_evo2_20250617.py -SEQ_LENGTH $sl -WINDOW_SIZE $ws -MODEL_SIZE $MODEL_SIZE"
    python run_evo2_20250617.py --SEQ_LENGTH $sl --WINDOW_SIZE $ws --MODEL_SIZE $MODEL_SIZE
  done
done

docker run -it --gpus all -p 8888:8888 helical_with_evo_2 jupyter notebook --ip=0.0.0.0 --allow-root

#########################################
# Install via helical/models/evo2
# instructions: https://github.com/helicalAI/helical/blob/release/helical/models/evo_2/README.md

# Via the Docker image
git clone https://github.com/helicalAI/helical.git
cd helical/helical/models/evo_2
docker build -t helical_with_evo_2 .
docker run -it --gpus all helical_with_evo_2

# via conda

conda create -n helical-env-with-evo-2 python=3.11
conda activate helical-env-with-evo-2

conda install cuda-toolkit=12.4 -c nvidia

export CUDNN_PATH=$CONDA_PREFIX/lib/python3.11/site-packages/nvidia/cudnn
export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/lib/python3.11/site-packages/nvidia/nvtx/include

pip install torch==2.6.0
pip install "helical[evo-2]@git+https://github.com/helicalAI/helical.git"

git clone git@github.com:Zymrael/vortex.git
cd vortex
git checkout f243e8e
sed -i 's/torch==2.5.1/torch==2.6.0/g' pyproject.toml
make setup-full
cd ..

pip install torch==2.6.0 torchvision

# Example Usage For Getting Embeddings

from helical.models.evo_2 import Evo2, Evo2Config

evo2_config = Evo2Config(batch_size=1)

evo2 = Evo2(configurer=evo2_config)

sequences = ["ACGT" * 1000]

dataset = evo2.process_data(data)

embeddings = evo2.get_embeddings(dataset)
# Get the last embedding of each sequence
print(embeddings["embeddings"][0][embeddings["original_lengths"][0]-1])
print(embeddings["embeddings"][1][embeddings["original_lengths"][1]-1])
print(embeddings["original_lengths"])


# Example Usage For Sequence Generation
from helical.models.evo_2 import Evo2, Evo2Config

evo2_config = Evo2Config(batch_size=1)

evo2 = Evo2(configurer=evo2_config)

sequences = ["ACGT" * 1000]

dataset = evo2.process_data(data)

generate = evo2.generate(dataset)

# Print the generated sequences
print(generate)