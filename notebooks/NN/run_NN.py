# Train an artifical neural network for supervised classification of BRCA1 variants using Evo2 embedding
# July 24, 2025
# Purpose: 
        # 1. Can we replicate RovHer framework but using neural network? Does adding embeddings help?
        # 2. Compare bewteen delta and refvar embedding concatenation methods
        # 3. Compare between using clinvar labels and class labels (LOF vs FUNC/INT)
##############################################################################################

ANNO_COLS = "yes"
EMBED_COLS="no" # delta, refvar, no
REGION = "RovHer_chr17" # "RovHer_chr17", BRCA1_DATA, RovHer_BRCA1 RovHer_LDLR, "both" (BRCA1 + LDLR RVs)
LAYER="blocks.28.mlp.l3"
y_label="height_FDR" # "height_FDR", "clinvar" (0, 0.25, 0.5, 0.75,1); "class" (LOF, FUNC/INT)
i=1
VAR_WIN="128"
# CORES=5
# LAYER="blocks.28.mlp.l3"  # fixed
# REGION="RovHer_chr17"     # BRCA1_DATA, RovHer_BRCA1 or RovHer_LDLR, RovHer_chr17
# y_label="height_FDR"   # "clinvar" "class" "height_FDR"
# python3.11 /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/NN/run_NN.py \
#       --REGION $REGION \
#       --LAYER $LAYER \
#       --COMBO $COMBO \
#       --Y_LABEL $y_label \
#       --EMBED_COLS $EMBED_COLS \
#       --ANNO_COLS $ANNO_COLS \
#       --CORES $CORES \
#       --i $i


##############################################################################################

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import sys
import argparse
import multiprocessing
import subprocess

parser = argparse.ArgumentParser(description="Evo2 embeddings")
parser.add_argument("--REGION", type=str, required=True, help="BRCA1_DATA, RovHer_BRCA1 or RovHer_LDLR, both (BRCA1 + LDLR RVs)")
parser.add_argument("--LAYER", required=True,type=str, help="embedding layer")
parser.add_argument("--Y_LABEL", type=str, required=True, help="clinvar (0, 0.25, 0.5, 0.75,1); class (LOF, FUNC/INT)")
parser.add_argument("--EMBED_COLS", required=True, type=str, help="yes or no")
parser.add_argument("--ANNO_COLS", required=True, type=str, help="yes or no")
parser.add_argument("--MODEL_SIZE", type=str, help="7B or 40B")
parser.add_argument("--i", type=int, help="number of iterations")
parser.add_argument("--EPOCH", type=int, help="number of EPOCHs")
parser.add_argument("--VAR_WIN", type=int, help="number of surrounding variant embeddings to average")
parser.add_argument("--CORES", type=int, default=2, help="number of cores")

args = parser.parse_args()
ANNO_COLS = args.ANNO_COLS
EMBED_COLS = args.EMBED_COLS
REGION = args.REGION
LAYER = args.LAYER
y_label = args.Y_LABEL
i=args.i
epoch=args.EPOCH
VAR_WIN = args.VAR_WIN
MODEL_SIZE = args.MODEL_SIZE
CORES = args.CORES

# =================== Run sequentially ===================
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import sys
import argparse
import multiprocessing
import subprocess

BLKS=["12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]        # delta, refvar 
Y_LABEL="class" # "clinvar" "class" "height_FDR"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="no"        # yes or no
REGION="BRCA1_DATA"

BLKS=["12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]        # delta, refvar 
Y_LABEL="clinvar" # "clinvar" "class" "height_FDR"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="yes"        # yes or no
REGION="BRCA1_DATA"

BLKS=["12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]        # delta, refvar 
Y_LABEL="clinvar" # "clinvar" "class" "height_FDR"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="no"        # yes or no
REGION="BRCA1_DATA"

BLKS=["12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29"]        # delta, refvar 
Y_LABEL="clinvar" # "clinvar" "class" "height_FDR"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="yes"        # yes or no
REGION="BRCA1_DATA"
#------------------------------- July 27

BLKS=["6", "7", "8", "9", "10", "11"]  

Y_LABEL="clinvar" # "clinvar" "class" "height_FDR"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="yes"
REGION="BRCA1_DATA" # BRCA1_DATA, RovHer_BRCA

Y_LABEL="clinvar" # "clinvar" "class" "height_FDR"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="no"
REGION="BRCA1_DATA" # BRCA1_DATA, RovHer_BRCA1

Y_LABEL="clinvar" # "clinvar" "class" "height_FDR"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="no"
REGION="RovHer_BRCA1" # BRCA1_DATA, RovHer_BRCA1


BLKS=["17","18","19", "20","21", "22", "24", "24"]
Y_LABEL="clinvar" # "clinvar" "class" "height_FDR"
EMBED_COLS="refvar" # delta, refvar, no
ANNO_COLS="no"
REGION="BRCA1_DATA" # BRCA1_DATA, RovHer_BRCA1


#########  No embeddings, just functional annotations #########
BLKS=["28"]
VAR_WINS = ["1"]

REGION="BRCA1_DATA"
Y_LABEL="class" # "clinvar"
EMBED_COLS="delta"   # delta, refvar, no
ANNO_COLS="yes"

REGION="BRCA1_DATA"
Y_LABEL="class" # "clinvar" "class"
EMBED_COLS="delta"   # delta, refvar, no
ANNO_COLS="yes"

REGION="BRCA1_DATA"
Y_LABEL="class" # "clinvar"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="yes"

REGION="BRCA1_DATA"
Y_LABEL="clinvar" # "clinvar"
EMBED_COLS="refvar"   # delta, refvar, no
ANNO_COLS="yes"


REGION="RovHer_BRCA1"
Y_LABEL="clinvar" # "clinvar"
EMBED_COLS="no"   # delta, refvar, no
ANNO_COLS="yes"

#-------------------------------
CORES=12
EPOCH=80
iterations=6
def run_script(blk, var_win, cores, i):
    LAYER = f"blocks.{blk}.mlp.l3"
    cmd = (
        f"python3.11 /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/NN/rovher_NN2.py "
        f"--REGION {REGION} "
        f"--LAYER {LAYER} "
        f"--Y_LABEL {Y_LABEL} "
        f"--EMBED_COLS {EMBED_COLS} "
        f"--ANNO_COLS {ANNO_COLS} "
        f"--EPOCH {EPOCH} "
        f"--i {i} "
        f"--VAR_WIN {var_win}"
    )
    print(f"\nRUN : {cmd}\n")
    os.system(cmd)

if __name__ == "__main__":
    for var_win in VAR_WINS:  # Iterate over different VAR_WIN values
        for blk in BLKS:  # Iterate over block layers
            for i in range(0, iterations + 1):
                run_script(blk, var_win, CORES, i)










#------------------------------- July 30 #-------------------------------

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import sys
import argparse
import multiprocessing
import subprocess

BLKS=["28"]
Y_LABEL="height_FDR" # "clinvar" "class" "height_FDR"
EMBED_COLS="no"   # delta, refvar, no
ANNO_COLS="yes"
REGION="RovHer_chr17"
VAR_WINS=["1"]

#-------------------------------
CORES=12
EPOCH=80
iterations=6
def run_script(blk, var_win, cores, i):
    LAYER = f"blocks.{blk}.mlp.l3"
    cmd = (
        f"python3.11 /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/NN/rovher_NN2.py "
        f"--REGION {REGION} "
        f"--LAYER {LAYER} "
        f"--Y_LABEL {Y_LABEL} "
        f"--EMBED_COLS {EMBED_COLS} "
        f"--ANNO_COLS {ANNO_COLS} "
        f"--EPOCH {EPOCH} "
        f"--i {i} "
        f"--VAR_WIN {var_win}"
    )
    print(f"\nRUN : {cmd}\n")
    os.system(cmd)

if __name__ == "__main__":
    for var_win in VAR_WINS:  # Iterate over different VAR_WIN values
        for blk in BLKS:  # Iterate over block layers
            for i in range(0, iterations + 1):
                run_script(blk, var_win, CORES, i)


# note: can remove rovher_NN.py 
# =================== Run in parallel (all runs) ===================


BLKS=["12", "13", "14", "15", "16", "17", "18"]        # delta, refvar 
BLKS=["13"]        # delta, refvar 
EMBED_COLS="refvar"       # yes or no
ANNO_COLS="no"        # yes or no
VAR_WIN="128"
REGION="BRCA1_DATA"
CORES = 15

EPOCH=100
Y_LABEL="class" # "clinvar" "class" "height_FDR"
i=1
def run_script(args):
    """Execute the external script with different parameters."""
    blk, combo = args  # Unpack the arguments
    LAYER = f"blocks.{blk}.mlp.l3"
    cmd = (
        f"python3.11 /mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/NN/rovher_NN2.py "
        f"--REGION {REGION} "
        f"--LAYER {LAYER} "
        f"--Y_LABEL {Y_LABEL} "
        f"--EMBED_COLS {EMBED_COLS} "
        f"--ANNO_COLS {ANNO_COLS} "
        f"--EPOCH {EPOCH} "
        f"--i {i} "
        f"--VAR_WIN {VAR_WIN}"
    )
    print(f"\n\n{cmd}\n\n")
    os.system(cmd)

if __name__ == "__main__":
    args_list = [(blk, blk) for blk in BLKS]
    with multiprocessing.Pool(processes=len(BLKS)) as pool:
        pool.map(run_script, args_list)

