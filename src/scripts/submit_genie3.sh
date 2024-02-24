#!/bin/bash

#SBATCH --job-name=genie3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --time=240:00:00
#SBATCH -o genie3_%j.log
#SBATCH -e genie3_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_long
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate pyscenic

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/

python /lustre/groups/ml01/workspace/samantha.bening/Bachelor/src/scripts/run_genie3.py

echo "Done!"

mamba deactivate