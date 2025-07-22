# Train MARS using 5-fold cross validation 
# Adding Evo2 7B model embeddings
# June 26, 2025

# Trial name scenarios:
# A      Train on variants of both genes
# B      Train on variants of BRCA1, test on LDLR 
# C      Train on variants of LDLR, test on BRCA1 
# D      Train on BRCA1, test on held-out set of BRCA1 
# E      Train on LDLR, test on held-out set of LDLR 


####################### TRAINING NEW MARS #######################
trial_name="A B C D E" #  B C D E
d=3
np=2000 # Reducing nprune reduces exhaustive search time; maximum number of permissible terms in the final pruned model.
nv=0.1
AF=0
cores=30
anno_cols="yes"
embedding_col="no" # "delta" (delta embeddings), "no" (no embeddings), "refvar" (ref + var embeddings)
blks="28" # 22 23 24 26 27 29 | 10 15 20 25 28 31 (or not accounted for if embedding_col="no")
# nk (maximum number of  forward pass terms)
script_path="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/MARS"

for t in $trial_name; do
for blk in $blks; do
  Rscript ${script_path}/train_MARS.r $AF $d $np $nv $t $anno_cols $embedding_col $blk $cores 
done
done

####################### TRAINING ROVHER #######################
trial_name="Chr17" #  B C D E
d=3
np=15 # Reducing nprune reduces exhaustive search time; maximum number of permissible terms in the final pruned model.
nv=0.1
AF=0.00003
cores=10
anno_cols="yes"
embedding_col="delta" # FIXED
blks="28" # (not accounted for if embedding_col="no")
script_path="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/MARS"

for t in $trial_name; do
for blk in $blks; do
  Rscript ${script_path}/train_cv_MARS.r $AF $d $np $nv $t $anno_cols $embedding_col $blk $cores 
done
done

####################### PLOTTING ##########################

trial_name="B" # 
d=3
np=2000 # Reducing nprune reduces exhaustive search time; maximum number of permissible terms in the final pruned model.
nv=0.1
AF=0
anno_cols="yes"
cores=2
for t in $trial_name; do
Rscript ${script_path}/plot_MARS.r $AF $d $np $nv $t $anno_cols $embedding_cols
done
