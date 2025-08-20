#!/bin/bash

#SBATCH --job-name=SRunner
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH --time=3:00:00
#SBATCH -o SRunner_%j.log
#SBATCH -e SRunner_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

START_TIME=$(date +%s)

mamba activate pyscenic_pipeline

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/

echo "Environment activated, starting script..."

srun python grnboost2_slurmrunner.py \
  --data_file data2/veo_ibd_balanced.h5ad \
  --celltype Myeloid \
  --output src/SCENICFiles/tester/ \

# Print elapsed time
ELAPSED=$(($(date +%s) - START_TIME))
printf "Time elapsed: %s\n\n" "$(date -d@$ELAPSED -u +%H\ hours\ %M\ min\ %S\ sec)"


echo "DONE"