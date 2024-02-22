from dask.distributed import Client
import dask.dataframe as dd
import dask.multiprocessing
import subprocess

subprocess.call("pyscenic grn SCENICfiles/data_filtered_scenic.loom SCENICfiles/genes.txt -o SCENICfiles/gene_gene_adj.csv --num_workers 20", shell=True)