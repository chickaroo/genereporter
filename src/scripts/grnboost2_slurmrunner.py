#!/usr/bin/env python3

import argparse
import scanpy as sc
import logging
import os
from dask_jobqueue.slurm import SLURMRunner
from dask.distributed import Client, get_worker
from arboreto.algo import grnboost2

def preprocess(data_file, celltype):
    adata = sc.read_h5ad(data_file)
    adata = adata[adata.obs['celltype_l1'] == celltype]
    sc.pp.filter_genes(adata, min_cells=30)
    ex_matrix = adata.to_df(layer='raw')
    return ex_matrix

def worker_debug_fn(x):
    logger = logging.getLogger("distributed.worker")
    worker = get_worker()
    logger.info(f"{worker.name}: processing task with input {x}")
    return x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_file", required=True, default = 'data2/veo_ibd_balanced.h5ad')
    parser.add_argument("--celltype", required=True, default='Myeloid')
    parser.add_argument("--output", required=True, default='src/SCENICfiles/tester')
    args = parser.parse_args()

    wdir = "/lustre/groups/ml01/workspace/christopher.lance/genereporter/"
    os.chdir( wdir )

    # Launch all processes (scheduler, client, workers)
    runner = SLURMRunner(scheduler_file="scheduler-{job_id}.json")

    with runner:
        with Client(runner) as client:
            # Only the client process continues past this point.
            client.wait_for_workers(runner.n_workers)
            print("Dask client ready with workers")

            # Preprocessing happens only once, in the client process:
            ex_matrix = preprocess(args.data_file, args.celltype)
            print("Expression matrix preprocessing done")

            # Run a dummy test submission to trace workers
            futures = client.map(worker_debug_fn, range(5))
            client.gather(futures)

            print("Launching GRNBoost2...")
            network = grnboost2(
                expression_data=ex_matrix,
                tf_names='all',
                client_or_address=client, 
                verbose=True
            )

            print("GRNBoost2 done")
            # filter for only importance >= 0.001 
            network = network[network['importance'] >= 0.001]

            network.to_csv(f'{args.output}/gg_adj/gg_adj_{args.celltype.replace(" ", "_")}.csv',  header=True, index=False)
            print(f"\tGene-gene adjacencies for {args.celltype} written to file")
