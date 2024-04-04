#!/bin/bash

#SBATCH --job-name=aucell2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300GB
#SBATCH --time=12:00:00
#SBATCH -o aucell2_%j.log
#SBATCH -e aucell2_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba deactivate
mamba activate pyscenic_notebook

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/src/scripts/

python run_regulons.py

echo "AUCell done!"

mamba deactivate
