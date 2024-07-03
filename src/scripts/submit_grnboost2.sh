#!/bin/bash

#SBATCH --job-name=GRNtcell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=250GB
#SBATCH --time=12:00:00
#SBATCH -o grnboosttcell_%j.log
#SBATCH -e grnboosttcell_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate pyscenic_pipeline

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/src/scripts

python run_grn_inference.py --celltype T Cell
#python run_grn_tf.py

echo "Done!"

mamba deactivate

