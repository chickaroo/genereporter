import subprocess
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Run GRN inference')
parser.add_argument('--input_loom', type=str, help='Input loom file')
parser.add_argument('--input_TFs', type=str, help='Input TF/gene file')
parser.add_argument('--output_adj', type=str, help='Output adjacencies .csv file')

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = -1

args = parser.parse_args()

subprocess.call(f"pyscenic grn {args.input_loom} {args.input_TFs} -o {args.output_adj} --num_workers 15", shell=True)