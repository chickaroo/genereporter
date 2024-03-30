#!/bin/bash

#SBATCH --job-name=grnboost
#SBATCH --nodes=1
#SBATCH --cpus-per-task=31
#SBATCH --mem=300GB
#SBATCH --time=72:00:00
#SBATCH -o grnboost_%j.log
#SBATCH -e grnboost_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate grnboost

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/src/scripts/

python run_grn_inference.py --celltype B Cell

echo "Done!"

mamba deactivate

