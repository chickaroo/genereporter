# import dependencies
import os
import numpy as np
import scanpy as sc
import loompy as lp

if __name__ == '__main__':

    # set a working directory
    wdir = "/lustre/groups/ml01/workspace/samantha.bening/Bachelor/"
    os.chdir( wdir )

    # # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
    f_loom_path_scenic = "src/SCENICfiles/new/data_filtered_scenic.loom"

    adata = sc.read_h5ad('data2/veo_ibd_balanced.h5ad')

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