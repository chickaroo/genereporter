#!/bin/bash

#SBATCH --job-name=EP40grnb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300GB
#SBATCH --time=12:00:00
#SBATCH -o grnboostEP40_%j.log
#SBATCH -e grnboostEP40_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate grnboost

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/src/scripts/

python run_grn_inference.py --celltype B Cell
#python run_grn_tf.py

echo "Done!"

mamba deactivate

