from pathlib import Path
import os
from pyscenic.prune import df2regulons
from pyscenic.aucell import aucell
import pandas as pd 
pd.options.mode.chained_assignment = None  # default='warn'
import anndata as ad
import pickle

print('Started')

# read in expression matrix 
# set a working directory
wdir = "/lustre/groups/ml01/workspace/samantha.bening/"
os.chdir( wdir )

adata = ad.read_h5ad('Bachelor/data2/veo_ibd_balanced.h5ad')
# make expression matrix 
ex_matrix = adata.to_df()
print("adata df read in")

df = pd.read_csv("Bachelor/src/SCENICfiles/new/reg_full10k.csv")
#print(df.head())

regulons = df2regulons(df)
#print(regulons[:2])

#with open("Bachelor/src/SCENICfiles/new/reg_full_df2regulons.p", "wb") as f:
#    regulons = pickle.dump(regulons, f)

#print("regulon pickle dumped")

auc_mtx = aucell(ex_matrix, regulons, num_workers=20)
print("AUCell method done! ")
auc_mtx.to_csv("Bachelor/src/SCENICfiles/new/auc_mtx.csv")