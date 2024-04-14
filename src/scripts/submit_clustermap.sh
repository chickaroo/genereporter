#!/bin/bash

#SBATCH --job-name=clustermap
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=150GB
#SBATCH --time=01:00:00
#SBATCH -o clustermap_%j.log
#SBATCH -e clustermap_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate grnboost

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/

python src/scripts/clustermap.py

echo "Clustermap done!"

mamba deactivate