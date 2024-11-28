#!/bin/bash

#SBATCH --job-name=sce_pipe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=200GB
#SBATCH --time=12:00:00
#SBATCH -o sce_pipe_%j.log
#SBATCH -e sce_pipe_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
#SBATCH --array=0-4

source $HOME/.bashrc

START_TIME=$(date +%s)

mamba deactivate
mamba activate pyscenic_pipeline

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/

# 0. Load params

helpFunction()
{
   echo ""
   echo "Usage: $0 -a parameterA -o parameterB -s parameterC"
   echo "Please format all paths or file names without leading or trailing slashes."
   echo "The working directory is always /lustre/groups/ml01/workspace/christopher.lance/genereporter"
   echo "You may set your input adata file and output directory within the working directory here: "
   echo -e "\t-a Name of adata file in /lustre/groups/ml01/workspace/christopher.lance/data2 (e.g. veo_ibd_balanced.h5ad)"
   echo -e "\t-o Path to existing output directory in /lustre/groups/ml01/workspace/christopher.lance/genereporter (e.g. src/SCENICfiles/new)"
   echo -e "\t-s Number of cells in subset to calculate (e.g. 10000)"
   exit 1 # Exit script after printing help
}

while getopts "a:o:s:" opt
do
   case "$opt" in
      a ) parameterA="$OPTARG" ;;
      o ) parameterO="$OPTARG" ;;
      s ) parameterS="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ] || [ -z "$parameterO" ] || [ -z "$parameterS" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

X_PARAM=( 'T Cell' 'Stroma' 'Epithelium' 'B Cell' 'Myeloid' )
INDEX=${SLURM_ARRAY_TASK_ID}

# Begin script in case all parameters are correct

# 1. Run GRN inference

echo "STARTING GRNBoost2 CELL TYPE SPECIFIC PIPELINE"
echo "1. Running GRN inference"

python src/scripts/run_grn_celltype.py --data $parameterA --output $parameterO --subset $parameterS --celltype ${X_PARAM[$INDEX]}

echo -e "1. DONE: GRN Inference\n\n"