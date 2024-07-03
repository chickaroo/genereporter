#!/bin/bash

#SBATCH --job-name=regclean
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --time=12:00:00
#SBATCH -o regclean_%j.log
#SBATCH -e regclean_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba deactivate
mamba activate pyscenic_small

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/

python src/scripts/run_regulons.py

echo "Regulon cleaning done!"

mamba deactivate
