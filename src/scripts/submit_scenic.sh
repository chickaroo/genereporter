#!/bin/bash

#SBATCH --job-name=pyscenic
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH -o pyscenic_%j.log
#SBATCH -e pyscenic_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate pyscenic

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/

pyscenic ctx SCENICfiles/adj.csv \
    /lustre/groups/ml01/workspace/samantha.bening/data/scenic_dbs/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /lustre/groups/ml01/workspace/samantha.bening/data/scenic_dbs/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname /lustre/groups/ml01/workspace/samantha.bening/data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname SCENICfiles/data_filtered_scenic.loom \
    --output SCENICfiles/reg.csv \
    --mask_dropouts \
    --num_workers 20

echo "Done!"

mamba deactivate

