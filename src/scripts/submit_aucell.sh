#!/bin/bash

#SBATCH --job-name=aucell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --time=12:00:00
#SBATCH -o aucell_%j.log
#SBATCH -e aucell_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba deactivate
mamba activate decoupler_env

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/

python src/scripts/run_aucell.py

echo "AUCell done!"

mamba deactivate