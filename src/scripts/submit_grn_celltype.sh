#!/bin/bash

#SBATCH --job-name=optimize
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=5:00:00
#SBATCH -o optimize_%j.log
#SBATCH -e optimize_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

START_TIME=$(date +%s)

mamba activate pyscenic_pipeline

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/

echo "Environment activated, starting script..."

python src/scripts/run_grn_celltype.py \
    --celltype Stroma \
    --data veo_ibd_balanced.h5ad \
    --output src/SCENICfiles/tester \
    --cluster local \


# Print elapsed time
ELAPSED=$(($(date +%s) - START_TIME))
printf "Time elapsed: %s\n\n" "$(date -d@$ELAPSED -u +%H\ hours\ %M\ min\ %S\ sec)"


echo "DONE"

