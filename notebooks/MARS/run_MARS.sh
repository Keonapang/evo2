# Train MARS using 5-fold cross validation 
# Adding Evo2 7B model embeddings
# June 26, 2025

# Trial name scenarios:
# A      Train on variants of both genes
# B      Train on variants of BRCA1, test on LDLR 
# C      Train on variants of LDLR, test on BRCA1 
# D      Train on BRCA1, test on held-out set of BRCA1 
# E      Train on LDLR, test on held-out set of LDLR 


####################### TRAINING MARS for CLINVAR #######################
anno_cols_list="no yes"
embed_cols="refvar" # "delta" (delta embeddings), "no" (no embeddings), "refvar" (ref + var embeddings)
# blks="20"
blks="20 19 21 18 17"
VAR_WIN="128"

anno_cols_list="yes"
embed_cols="delta" # "delta" (delta embeddings), "no" (no embeddings), "refvar" (ref + var embeddings)
blks="20 19 21 18 17"

anno_cols_list="no yes"
embed_cols="refvar"
VAR_WIN="128"
blks="16 25"
blks="14"

blks="12 13 14 15 16 17 18 19 20" # 22 23 24 26 27 29 | 10 15 20 25 28 31 (or not accounted for if embedding_col="no")
blks="21 22 23 24 25 26 27 28 29"

trial_name="D" #  B C D E (uses RovHer_BRCA1.txt)
d=3
np=200 # Reducing nprune reduces exhaustive search time; maximum number of permissible terms in the final pruned model.
nv=0.1
AF=0
reverse="no" # "yes or "no" (to include reverse complement embeddings)
cores=12
for t in $trial_name; do
for anno_cols in $anno_cols_list; do
  for blk in $blks; do
    script_path="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/MARS"
    Rscript ${script_path}/train_MARS_clinvar.r $AF $d $np $nv $t $anno_cols $embed_cols $blk $cores $reverse $VAR_WIN
  done
done
done



####################### TRAINING ROVHER on Chr17 RVs #######################
# ONLY embedding layer available is 28

# y variable is "height_FDR"
trial_name="chr17" #  B C D E
d=3
np=150 # 50, 20000, Reducing nprune reduces exhaustive search time; maximum number of permissible terms in the final pruned model.
nv=0.1
AF=0 # 0.00003 or 0
anno_cols="yes" # no or yes (to include annotation cols)
embedding_col="refvar" # "delta" , "no" (no embeddings), "refvar"
pval_thresh=0.005 # numeric (linreg pre-filtering to keep only sig. predictors to pval < pval_thresh)
reverse="yes" # "yes or "no" (to include reverse complement embeddings)
blks="28" # (not accounted for if embedding_col="no")
cores=15
VARWIN="1" # "128" or "256"
script_path="/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/notebooks/MARS"

for t in $trial_name; do
for blk in $blks; do
  Rscript ${script_path}/train_cv_MARS_new.r $AF $d $np $nv $t $anno_cols $embedding_col $blk $cores $reverse $pval_thresh $VARWIN
done
done

for t in $trial_name; do
for blk in $blks; do
  Rscript ${script_path}/train_cv_MARS_new2.r $AF $d $np $nv $t $anno_cols $embedding_col $blk $cores $reverse $pval_thresh $VARWIN
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
