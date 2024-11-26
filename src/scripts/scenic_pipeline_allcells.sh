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

source $HOME/.bashrc

mamba deactivate
mamba activate pyscenic_pipeline

cd /lustre/groups/ml01/workspace/christopher.lance/genereporter/

# accept parameters: 
# Input adata file
# Directory for output files
# Number of cells in subset (or all cells)

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

# Begin script in case all parameters are correct



# 1. Run GRN inference

python src/scripts/run_grn_tf.py --data $parameterA --output $parameterO --subset $parameterS

echo "1. done: GRN inference"




# 2. Make loom file 

python src/scripts/make_loom.py --data $parameterA --output $parameterO 

echo "2. done: make loom file"




# 3. Run CisTarget

mamba deactivate
mamba activate pyscenic_small

pyscenic ctx "$parameterO"/TFtg_adj.csv \
    /lustre/groups/ml01/workspace/christopher.lance/genereporter/data/scenic_dbs/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /lustre/groups/ml01/workspace/christopher.lance/genereporter/data/scenic_dbs/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname /lustre/groups/ml01/workspace/christopher.lance/genereporter/data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname "$parameterO"/data_filtered_scenic.loom \
    --min_genes=15 \
    --output "$parameterO"/regulons_output.csv \
    --num_workers 29

echo "3. done: cisTarget regulon inference"




# 4. Run AUCell

mamba deactivate
mamba activate decoupler_env

python src/scripts/run_aucell.py --data $parameterA --output $parameterO 

echo "4. done: AUCell"
echo ""
echo "All SCENIC output files are in $parameterO"