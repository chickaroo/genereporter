# import dependencies
import os
import numpy as np
import scanpy as sc
import loompy as lp
from argparse import ArgumentParser

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("--data", type=str, default = 'veo_ibd_balanced.h5ad')
    parser.add_argument("--output", type=str, default = 'src/SCENICfiles/new')
    args = parser.parse_args()
    data_file = f'data2/{args.data}'
    output_dir = args.output
    print("args read in")

    # set a working directory
    wdir = "/lustre/groups/ml01/workspace/christopher.lance/genereporter"
    os.chdir( wdir )

    # # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
    f_loom_path_scenic = f"{output_dir}/data_filtered_scenic.loom"

    adata = sc.read_h5ad(data_file)
    print("adata read in")

    # create basic row and column attributes for the loom file:
    row_attrs = {
        "Gene": np.array(adata.var_names) ,
    }
    col_attrs = {
        "CellID": np.array(adata.obs_names) ,
        "nGene": np.array( np.sum(adata.layers["raw"].transpose()>0 , axis=0)).flatten() ,
        "nUMI": np.array( np.sum(adata.layers['raw'].transpose() , axis=0)).flatten() ,
    }
    lp.create( f_loom_path_scenic, adata.layers['raw'].transpose(), row_attrs, col_attrs)
    print("loom file created")