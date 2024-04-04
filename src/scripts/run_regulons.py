import os
import glob
import pickle
import pandas as pd
import anndata as ad

from dask.diagnostics import ProgressBar

#from distributed import Client, LocalCluster

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons


if __name__ == '__main__':
    # read in expression matrix 
    # set a working directory
    wdir = "/lustre/groups/ml01/workspace/samantha.bening/"
    os.chdir( wdir )

    adata = ad.read_h5ad('Bachelor/data2/veo_ibd_balanced.h5ad')
    # make expression matrix 
    ex_matrix = adata.to_df()

    # load ranking databases
    db_fnames = glob.glob("data/scenic_dbs/hg38_*.genes_vs_motifs.rankings.feather")
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    # load adjacencies
    adjacencies = pd.read_csv("Bachelor/src/SCENICfiles/TFtg_adj.csv", names=['TF', 'target', 'importance'])

    # load adjacencies
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    print(modules[:20])

    # set up Dask cluster
    # create custom Dask client
    #local_cluster = LocalCluster(n_workers=31, # put in one less than the number of cores you gave the job
    #                            threads_per_worker=5, 
    #                            processes=True,
    #                            memory_limit="10GiB") 
    #custom_client = Client(local_cluster)

    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(dbs, modules,
                    "data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", 
                    num_workers=30)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)

    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv("Bachelor/src/SCENICfiles/motifs.csv")
    with open("Bachelor/src/SCENICfiles/regulons.p", "wb") as f:
        pickle.dump(regulons, f)

    #custom_client.shutdown()
    #local_cluster.close()