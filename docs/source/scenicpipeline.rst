SCENIC Pipeline
--------------------------------

This is a detailed guide on how to use the provided scripts to run the SCENIC pipeline on a SLURM cluster.

The ``scenic_pipeline.sh`` bash script runs through the three step SCENIC pipeline automatically, based on these parameters. The parameters are not optional.

Please note that currently, the working directory is hard-coded for our specific SLURM cluster and environment. You may need to modify this to fit your own cluster's specifications.

**PARAMETERS:**

**-a (adata file)**
   Name of adata file in the working directory. 
   
   Example: ``data2/veo_ibd_balanced.h5ad``

**-o (output directory)** 
   Path to existing output directory in the working directory.
   
   Example: ``src/SCENICfiles/new``

**-s (subset size)**
   Number of cells in subset to calculate. 
   
   Example: ``20000``
   
   **Note:** Pass ``0`` to NOT subset and keep entire data.

Output Files
------------

The output files are as follows:

* ``TFtg_adj.csv.gz``: Adjacencies file
* ``regulons_output.csv.gz``: Regulons file
* ``veo_ibd_balanced_aucell_final.h5ad``: Adata file with AUCell values
* ``metadata.txt``: Metadata file with details about run