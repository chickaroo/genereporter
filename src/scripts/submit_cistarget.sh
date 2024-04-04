#!/bin/bash

#SBATCH --job-name=cistarget
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300GB
#SBATCH --time=12:00:00
#SBATCH -o cistarget_%j.log
#SBATCH -e cistarget_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate pyscenic_notebook

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/src/scripts/

python run_regulons.py

echo "Regulons done!"

mamba deactivate
