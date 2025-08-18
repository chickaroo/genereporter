import os
import numpy as np
import scanpy as sc
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
import dask.array as da
from dask_jobqueue import SLURMCluster
from argparse import ArgumentParser

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("--celltype", type=str, nargs='+')
    parser.add_argument("--data", type=str, default = 'veo_ibd_balanced.h5ad')
    parser.add_argument("--output", type=str, default = 'src/SCENICfiles/new')
    parser.add_argument("--subset", type=int, default = 20000)
    parser.add_argument("--cluster", type=str, default = "distributed")
    args = parser.parse_args()
    data_file = f'data2/{args.data}'
    output_dir = args.output
    subset_size = args.subset
    celltype = " ".join(args.celltype)
    print("\tArgs read in")

    # read in expression matrix 
    # set a working directory
    wdir = "/lustre/groups/ml01/workspace/christopher.lance/genereporter/"
    os.chdir( wdir )

    adata = sc.read_h5ad(data_file)
    print("\tAdata read in")

    # alternative: create distributed dask client over multiple parallel jobs

    if args.cluster == "distributed":
        print("Using distributed cluster")
        # create a SLURM cluster
        cluster = SLURMCluster(
            queue='cpu_p',                    # --partition=cpu_p
            cores=32,                          # --cpus-per-task=8
            memory="100GB",                    # --mem=20G
            walltime="6:00:00",              # --time=12:00:00
            job_name="scenicdist",            # --job-name=scenicplus
            job_extra_directives=[
                "--nodes=1",                  # --nodes=1
                "--qos=cpu_normal",           # --qos=cpu_normal
                "--nice=10000",               # --nice=10000
                "-o scenicdist_%j.log",       # stdout log file (overrides default)
                "-e scenicdist_%j.err"        # stderr log file (overrides default)
            ]
        )
        cluster.scale(5)  # activate 10 workers
    else:
        # create local Dask client
        print("Using local cluster")
        cluster = LocalCluster(n_workers=30, # put in one less than the number of cores you gave the job
                                    threads_per_worker=1) 

        
    custom_client = Client(cluster) # create distributed client


    # subset for a very broad cell type
    adata = adata[adata.obs['celltype_l1'] == celltype]
    # filter for genes not expressed in e.g. 30 or more cells
    sc.pp.filter_genes(adata, min_cells=30)
    # randomly sample {subset_size} number of cells from this cell type (compute limit) 
    # NOT for myeloid: only has total of 8k cells 
    if celltype != 'Myeloid':
        a = np.zeros(int(adata.obs['celltype_l1'].value_counts()[celltype]), dtype=int)
        a[:subset_size] = 1
        np.random.shuffle(a)
        a = a.astype(bool)
        adata = adata[a, :]
    # make expression matrix 
    ex_matrix = adata.to_df(layer='raw') # use raw layer here
    print(f"Expression matrix size: {ex_matrix.shape}")
    # try to convert to dask array
    ex_matrix_da = da.from_array(ex_matrix.values, chunks=(2000, ex_matrix.shape[1]))
    ex_matrix_da = ex_matrix_da.persist()  # Persist the array in memory

    print("\tPreprocessing done")


    # run GRNBoost2
    network = grnboost2(expression_data=ex_matrix_da,
                        gene_names=ex_matrix.columns,
                        tf_names='all', # gene-gene adjacencies
                        client_or_address=custom_client, 
                        verbose=True)

    print("\tGRNBoost2 done")
    # filter for only importance >= 0.001 
    network = network[network['importance'] >= 0.001]

    network.to_csv(f'{output_dir}/gg_adj/gg_adj_{celltype.replace(" ", "_")}.csv',  header=True, index=False)
    print(f"\tGene-gene adjacencies for {celltype} written to file")