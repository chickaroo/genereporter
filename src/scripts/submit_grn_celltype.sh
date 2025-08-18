#!/bin/bash

#SBATCH --job-name=disttest
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH -o disttest_%j.log
#SBATCH -e disttest_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

START_TIME=$(date +%s)

mamba activate pyscenic_pipeline

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/

python src/scripts/run_grn_celltype.py \
    --celltype Myeloid \
    --data veo_ibd_balanced.h5ad \
    --output src/SCENICfiles/tester \
    --cluster distributed \



# Print elapsed time
ELAPSED=$(($(date +%s) - START_TIME))
printf "Time elapsed: %s\n\n" "$(date -d@$ELAPSED -u +%H\ hours\ %M\ min\ %S\ sec)"


echo "DONE"

