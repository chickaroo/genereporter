#!/bin/bash

#SBATCH --job-name=genie3_morememory
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=400GB
#SBATCH --time=240:00:00
#SBATCH -o genie3_upmem_%j.log
#SBATCH -e genie3_upmem_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_long
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate pyscenic2

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/

python src/scripts/run_grn_inference.py --input_loom SCENICfiles/data_filtered_scenic.loom --input_TFs SCENICfiles/genes.txt --output SCENICfiles/gene_gene_adj.csv

echo "Done!"

mamba deactivate