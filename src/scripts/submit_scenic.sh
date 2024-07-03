#!/bin/bash

#SBATCH --job-name=ctx
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=150GB
#SBATCH --time=12:00:00
#SBATCH -o ctx_%j.log
#SBATCH -e ctx_%j.err
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

source $HOME/.bashrc

mamba activate pyscenic_notebook

cd /lustre/groups/ml01/workspace/samantha.bening/Bachelor/

#python src/scripts/run_grn_inference.py --input_loom SCENICfiles/data_filtered_scenic.loom --input_TFs SCENICfiles/allTFs_hg38.txt --output SCENICfiles/adj2.csv

#echo "GRN done!"

pyscenic ctx src/SCENICfiles/new/TFtg_adj.csv \
    /lustre/groups/ml01/workspace/samantha.bening/data/scenic_dbs/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /lustre/groups/ml01/workspace/samantha.bening/data/scenic_dbs/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname /lustre/groups/ml01/workspace/samantha.bening/data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname src/SCENICfiles/new/data_filtered_scenic.loom \
    --min_genes=15 \
    --output src/SCENICfiles/new/reg_full10k.csv \
    --num_workers 29

echo "cisTarget Done!"

mamba deactivate

