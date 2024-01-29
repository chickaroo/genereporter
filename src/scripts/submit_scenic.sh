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
# ./my_program input.dat   # Lets simplfy your job script as following:
 
echo "Starting SCENIC Pipeline"

# Load modules
conda activate bachelor_env

# make f_db_names list

pyscenic ctx adj.tsv \ # TODO
    {f_db_names} \
    --annotations_fname {f_motif_path} \
    --expression_mtx_fname {f_loom_path_scenic} \
    --output SCENICfiles/reg.csv \
    --mask_dropouts \
    --num_workers 20

 
echo "Done with Steps 2-3"

pyscenic aucell \ #TODO
    {f_loom_path_scenic} \
    reg.csv \
    --output {f_pyscenic_output} \
    --num_workers 20


echo "Done with Step 4"

# TODO: Vizualize results in python script!! 