#!/bin/bash
#SBATCH --job-name=submit_scenic
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --output=submit_scenic_output.txt
#SBATCH --error=submit_scenic_error.txt
#SBATCH --nice=10000
 
# Run your job commands
 
echo "Starting SCENIC Pipeline"

# Load modules
conda activate bachelor_env

# test first 

pyscenic -h

conda deactivate

# # make f_db_names list
# pyscenic ctx adj.tsv \ 
#     $(ls SCENICfiles/mc_v10_clust/gene_based/*.feather) \
#     --annotations_fname SCENICfiles/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tb \
#     --expression_mtx_fname SCENICfiles/data_filtered_scenic.loom \
#     --output SCENICfiles/reg.csv \
#     --mask_dropouts \
#     --num_workers 20

 
# echo "Done with Steps 2-3"

# pyscenic aucell \ 
#     SCENICfiles/data_filtered_scenic.loom \
#     reg.csv \
#     --output SCENICfiles/pyscenic_output.loom \
#     --num_workers 20


# echo "Done with Step 4"

# TODO: Vizualize results in python script!! 