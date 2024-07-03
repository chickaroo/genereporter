#!/bin/bash

#SBATCH --job-name=makeloom
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=96GB
#SBATCH --time=4:00:00
#SBATCH -o loom_%j.log
#SBATCH -e loom_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate pyscenic_pipeline

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/

python make_loom.py

echo "Done!"

mamba deactivate