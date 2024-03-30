import os
import pandas as pd
import scanpy as sc

from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names
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
    local_cluster = LocalCluster(n_workers=30, # put in one less than the number of cores you gave the job
                                threads_per_worker=10) 
    custom_client = Client(local_cluster)


    # subset for a very broad cell type
    adata = adata[adata.obs['celltype_l1'] == celltype]
    # filter for genes not expressed in e.g. 30 or more cells
    sc.pp.filter_genes(adata, min_cells=30)
    # make expression matrix 
    ex_matrix = adata.to_df()
    # make 'tf_list' which is just the gene list here (gene-gene adjacencies)
    tf_names = list(adata.var_names)


    # run GRNBoost2
    network = grnboost2(expression_data=ex_matrix,
                        tf_names=tf_names,
                        client_or_address=custom_client)

    # filter for only importance >= 0.001 
    network = network[network['importance'] >= 0.001]

    network.to_csv(f'src/SCENICfiles/gg_adj_{celltype.replace(" ", "_")}.csv',  header=False, index=False)