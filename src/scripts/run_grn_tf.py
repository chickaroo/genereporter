import os
import numpy as np
import scanpy as sc
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from distributed import Client, LocalCluster
from argparse import ArgumentParser

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("--data", type=str, default = 'veo_ibd_balanced.h5ad')
    parser.add_argument("--output", type=str, default = 'src/SCENICfiles/new')
    parser.add_argument("--subset", type=int, default = 20000)
    args = parser.parse_args()
    data_file = f'data2/{args.data}'
    output_dir = args.output
    subset_size = args.subset
    print("Args read in")

    # read in expression matrix 
    # set a working directory
    wdir = "/lustre/groups/ml01/workspace/christopher.lance/genereporter/"
    os.chdir( wdir )

    adata = sc.read_h5ad(data_file)
    print("adata read in")

    # create custom Dask client
    local_cluster = LocalCluster(n_workers=20, # put in one less than the number of cores you gave the job
                                threads_per_worker=3, 
                                processes=True,
                                memory_limit="13GiB") 
    custom_client = Client(local_cluster)

    # filter for genes not expressed in e.g. 30 or more cells
    sc.pp.filter_genes(adata, min_cells=30)
    # randomly sample {subset_size} cells from this subset (compute limit)
    a = np.zeros(len(adata), dtype=int)
    a[:subset_size] = 1
    np.random.shuffle(a)
    a = a.astype(bool)
    adata = adata[a, :]
    # make expression matrix 
    ex_matrix = adata.to_df(layer='raw') # use raw layer here

    # load tf_names
    tf_names = load_tf_names("data/allTFs_human.txt")
    print("Preprocessing done")


    # run GRNBoost2
    network = grnboost2(expression_data=ex_matrix,
                        tf_names=tf_names, # TF-target gene adjacencies
                        client_or_address=custom_client, 
                        verbose=True)

    print("GRNBoost2 done")
    # filter for only importance >= 0.001 
    network = network[network['importance'] >= 0.001]

    network.to_csv(f'{output_dir}/TFtg_adj.csv',  header=True, index=False)
    print("TF-target gene adjacencies written to file")