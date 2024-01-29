# import dependencies
import os
import numpy as np
import pyscenic
import pandas as pd
import scanpy as sc
import loompy as lp
import glob
from MulticoreTSNE import MulticoreTSNE as TSNE
import argparse

# set variables for file paths to read from and write to:

parser = argparse.ArgumentParser(
    description='Run SCENIC pipeline on a dataset.'
)
parser.add_argument('-wdir', type=str, help='working directory of "Bachelor" repository')
parser.add_argument('-f_loom_path_scenic', type=str, default="SCENICfiles/data_filtered_scenic.loom",
                    help='path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.')
parser.add_argument('-f_anndata_path', type=str, default="data/output/adata.h5ad",
                    help='path to anndata object, which will be updated to store Scanpy results as they are generated below')
parser.add_argument('-f_adjacencies', type=str, default="SCENICfiles/adj.csv",)

parser.add_argument('-dir_db_ranking', type=str, default='SCENICfiles/mc_v10_clust/gene_based/')
parser.add_argument('-f_motif_path', type=str, default='SCENICfiles/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tb')

parser.add_argument('-f_pyscenic_output', type=str, default="SCENICfiles/pyscenic_output.loom",
                    help='path to pyscenic output')
parser.add_argument('-f_final_loom', type=str, default="SCENICfiles/scenic_integrated_output.loom",
                    help='loom output, generated from a combination of Scanpy and pySCENIC results')

# set a working directory
#wdir = "C:/Users/saman/00_Bachelorarbeit/Bachelor/"
wdir = parser.parse_args().wdir
os.chdir( wdir )

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = os.path.join(os.getcwd(), parser.parse_args().f_loom_path_scenic)

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = os.path.join(os.getcwd(), parser.parse_args().f_anndata_path)

# path to adjacency matrix
f_adjacencies = os.path.join(os.getcwd(), parser.parse_args().f_adjacencies)

# path to ranking databases
dir_db_ranking = os.path.join(os.getcwd(), parser.parse_args().dir_db_ranking)
f_db_names = ' '.join( glob.glob(dir_db_ranking) )

# path to motif databases
f_motif_path = os.path.join(os.getcwd(), parser.parse_args().f_motif_path)

# path to pyscenic output
f_pyscenic_output = os.path.join(os.getcwd(), parser.parse_args().f_pyscenic_output)

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = os.path.join(os.getcwd(), parser.parse_args().f_final_loom)

sc.settings.verbosity = 0 # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_versions()
sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = -1

# read in expression data 
adata = sc.read_h5ad(f_anndata_path)

# STEP 1: GRN inference, generation of co-expression modules

# read in adjacency matrix
adjacencies = pd.read_csv(f_adjacencies, index_col=False)

# STEP 2-3: Regulon prediction aka cisTarget from CLI
