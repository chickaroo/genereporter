import os
from pyscenic.prune import df2regulons
import pandas as pd 
import pickle

print('Started')

# read in expression matrix 
# set a working directory
wdir = "/lustre/groups/ml01/workspace/samantha.bening/"
os.chdir( wdir )

df = pd.read_csv("Bachelor/src/SCENICfiles/new/reg_full10k.csv")
#print(df.head())

regulons = df2regulons(df)
#print(regulons[:2])

with open("Bachelor/src/SCENICfiles/new/reg_full_df2regulons.p", "wb") as f:
    regulons = pickle.dump(regulons, f)

print("regulon pickle dumped")

#auc_mtx = aucell(ex_matrix, regulons, num_workers=30)
#rint("AUCell method done! ")
#auc_mtx.to_csv("Bachelor/src/SCENICfiles/auc_mtx.csv")