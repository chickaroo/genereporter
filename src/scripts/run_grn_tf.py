import os
import numpy as np
import scanpy as sc
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from distributed import Client, LocalCluster
from dask.diagnostics import ProgressBar
from datetime import datetime
from argparse import ArgumentParser

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("--data", type=str, default = 'data2/veo_ibd_balanced.h5ad')
    parser.add_argument("--output", type=str, default = 'src/SCENICfiles/rerun_combined_data')
    parser.add_argument("--subset", type=int, default = 0) # default 0 = no subsetting, pass e.g. 10000 to subset to 10k cells
    args = parser.parse_args()
    data_file = args.data
    output_dir = args.output
    subset = args.subset
    print("\tArgs read in")

    # read in expression matrix 
    # set a working directory
    wdir = "cd"
    os.chdir( wdir )

    adata = sc.read_h5ad(data_file)
    print("\tAdata read in")

    # create custom Dask client
    local_cluster = LocalCluster(n_workers=8, # less workers to avoid memory issues
                                threads_per_worker=2, # threads for faster computation
                                processes=True, 
                                memory_limit="48GiB", # memory limit per worker (total 384GB on cluster of total 400GB)
                                silence_logs="debug") 
    custom_client = Client(local_cluster)

    # filter for genes not expressed in e.g. 30 or more cells
    sc.pp.filter_genes(adata, min_cells=30)

    if subset != 0:
        # randomly sample {subset} cells from this subset (compute limit) if parameter is not 0
        a = np.zeros(len(adata), dtype=int)
        a[:subset] = 1
        np.random.shuffle(a)
        a = a.astype(bool)
        adata = adata[a, :]
    # make expression matrix 
    ex_matrix = adata.to_df(layer='raw') # use raw layer here, no normalization

    # load tf_names
    tf_names = load_tf_names("data/allTFs_human.txt")
    print("\tPreprocessing done")


    # run GRNBoost2
    with ProgressBar():
        start_time = datetime.now()
        network = grnboost2(expression_data=ex_matrix,
                            tf_names=tf_names, # TF-target gene adjacencies
                            client_or_address=custom_client, 
                            verbose=True)
        end_time = datetime.now()

    print(f"GRN inference completed in {end_time - start_time}")
    # filter for only importance >= 0.001 
    network = network[network['importance'] >= 0.001]
    print(f"Generated network with {len(network)} regulatory links")

    network.to_csv(f'{output_dir}/TFtg_adj.csv',  header=True, index=False)
    print("\tTF-target gene adjacencies written to file")