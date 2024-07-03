import os
import numpy as np
import scanpy as sc
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
from argparse import ArgumentParser

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("--celltype", type=str, nargs='+')
    args = parser.parse_args()
    celltype = " ".join(args.celltype)

    # read in expression matrix 
    # set a working directory
    wdir = "/lustre/groups/ml01/workspace/samantha.bening/Bachelor/"
    os.chdir( wdir )

    adata = sc.read_h5ad('data2/veo_ibd_balanced.h5ad')

    # create custom Dask client
    local_cluster = LocalCluster(n_workers=24, # put in one less than the number of cores you gave the job
                                threads_per_worker=2, 
                                processes=True,
                                memory_limit="9GiB") 
    custom_client = Client(local_cluster)


    # subset for a very broad cell type
    adata = adata[adata.obs['celltype_l1'] == celltype]
    # filter for genes not expressed in e.g. 30 or more cells
    sc.pp.filter_genes(adata, min_cells=40)
    # randomly sample 5k cells from this subset (compute limit)
    a = np.zeros(int(adata.obs['celltype_l1'].value_counts()[celltype]), dtype=int)
    a[:10000] = 1
    np.random.shuffle(a)
    a = a.astype(bool)
    adata = adata[a, :]
    # make expression matrix 
    ex_matrix = adata.to_df(layer='raw') # use raw layer here

    # do dask distributed and have multiple jobs?? 


    # run GRNBoost2
    network = grnboost2(expression_data=ex_matrix,
                        tf_names='all', # gene-gene adjacencies
                        client_or_address=custom_client)

    # filter for only importance >= 0.001 
    network = network[network['importance'] >= 0.001]

    network.to_csv(f'src/SCENICfiles/new/gg_adj/gg_adj_{celltype.replace(" ", "_")}.csv',  header=False, index=False)