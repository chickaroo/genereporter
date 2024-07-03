#!/bin/bash

#SBATCH --job-name=GRNmstro
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=250GB
#SBATCH --time=12:00:00
#SBATCH -o grnbooststro_%j.log
#SBATCH -e grnbooststro_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate pyscenic_pipeline

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/src/scripts

python run_grn_inference.py --celltype Stroma
#python run_grn_tf.py

echo "Done!"

mamba deactivate

