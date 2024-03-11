import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import anndata
import math
import scipy.stats as stats
from scipy import interpolate
from scipy.interpolate import BSpline, splrep
import warnings
warnings.filterwarnings('ignore')

plt.rcParams["font.family"] = "monospace"
plt.rcParams["font.size"] = 10
plt.rcParams['figure.figsize'] = (8,8)

# methods here

# set a working directory
wdir = "/Users/samibening/Projects/Bachelor/"
os.chdir(wdir)

adata = anndata.read_h5ad("data/output/adata.h5ad")

def get_adata() -> anndata.AnnData:
    return adata


def plot_umap(adata: anndata.AnnData = adata, color: str = "celltypist_cell_label_coarse") -> None:
    sc.pl.umap(
    adata, color=[color],
    frameon=False, legend_loc="on data", title=color, legend_fontsize=9, legend_fontweight=600, 
    legend_fontoutline=1
)
    

# Classify gene expression levels in discrete categories
    
def clean_data(df, col, threshold = 99.75):
    data = df[col].to_numpy().flatten()
    data = data[data <= np.percentile(data, threshold)]
    return data

def find_thresholds(filtered):
    range = (abs(max(filtered)) + abs(min(filtered)))
    step = range / float(4)
    very_low = min(filtered) + step
    low = min(filtered) + step*2
    middle = min(filtered) + step*3
    #print(round(very_low, 2), round(low, 2), round(middle, 2), round(max(filtered), 2))
    return very_low, low, middle, max(filtered)

def classify_exp_level(df, filtered, col):
    very_low, low, middle, high = find_thresholds(filtered)
    def func(x):
        if x <= very_low:
            return "very_low"
        elif very_low < x <= low:
            return "low"
        elif low < x <= middle:
            return "middle"
        elif middle < x <= high:
            return "high"
        else:
            return "very_high"
    df['expr_class'] = df[col].apply(func)

    summary = "Quantile thresholds: \n"
    summary += str("very low: " + str(round(stats.percentileofscore(df[col], very_low), 4)) + ", low: " + str(round(stats.percentileofscore(df[col], low), 4)) + ", middle: " + 
            str(round(stats.percentileofscore(df[col], middle), 4)) + ", high: " + str(round(stats.percentileofscore(df[col], high), 4)) + ", very high: 99.7500\n")
   
    summary += "\nNumber of genes per category: \n"
    summary += str("very_low: " + str(len(df[df['expr_class'] == 'very_low']))+ '\n')
    summary += str("low: " + str(len(df[df['expr_class'] == 'low']))+ '\n')
    summary += str("middle: " + str(len(df[df['expr_class'] == 'middle']))+ '\n')
    summary += str("high: " + str(len(df[df['expr_class'] == 'high'])))+ '\n'
    summary += str("very_high: " + str(len(df[df['expr_class'] == 'very_high']))+ '\n')
    return df, summary


def make_df(adata, col='log1p(means)', layer='log_norm'):
    sc.pp.highly_variable_genes(adata, layer=layer)
    df = adata.var.sort_values(['means'])
    df['gene_num'] = range(len(df))
    df['log(means)'] = np.log(df['means'])
    df['log1p(means)'] = np.log1p(df['means'])

    filtered = clean_data(df, col=col, threshold=99.75) # can choose upper outlier threshold here
    df_new, summary = classify_exp_level(df =df, filtered=filtered, col=col)
    return df_new, summary

def explain_expr_celltypes(GOI, adata=adata, col='log1p(means)', layer='log_norm'):
    out = pd.DataFrame()
    
    for cell_type in adata.obs['celltypist_cell_label_coarse'].unique():
        subset = adata[adata.obs['celltypist_cell_label_coarse'] == cell_type]
        df = make_df(subset, col)
        out = pd.concat([out, df.loc[df.index == GOI]])
        
    out['cell_type'] = adata.obs['celltypist_cell_label_coarse'].unique()
    out = out[['cell_type', 'expr_class', col]]
    out[col] = out[col].apply(lambda x: round(x, 3))
    out['expr_class'] = out['expr_class'].apply(lambda x: x.replace('_', ' '))
    out = out.rename(columns={'cell_type': 'Cell type', 'expr_class': 'Expression class', col: 'Avg. expression over cell type'})
    return out


# TODO threshold for expecting dropout zeros 

def dropout_threshold(adata):
    df = adata.to_df()
    # finalized threshold to come here
    return df


# find overall expression of all genes, highlighting GOI

def plot_expr_class(GOI, ax, adata=adata, cell_type=None, col='log1p(means)'):
    title = "Overall expression across all cell types"
    if cell_type != None:
        title = str("Overall expression in " + cell_type + " cells")
        if cell_type in adata.obs['celltypist_cell_label_coarse'].unique():
            adata = adata[adata.obs['celltypist_cell_label_coarse'] == cell_type]
        elif cell_type in adata.obs['celltypist_cell_label'].unique():
            adata = adata[adata.obs['celltypist_cell_label'] == cell_type]
        else:
            print("Cell type not found. Please check spelling.")
            return
    
    df_new, sum = make_df(adata)

    g = sns.scatterplot(data=df_new, ax=ax, x='gene_num',  y=col, hue='expr_class', linewidth=0)
    annotation = str(GOI + " (" + df_new.loc[df_new.index == GOI]['expr_class'].values[0] + ")")
    highlight_y = df_new.loc[df_new.index == GOI][col]
    highlight_x = df_new.loc[df_new.index == GOI]['gene_num']
    g.scatter(highlight_x, highlight_y, color = 'yellow', linewidth=1)
    props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
    g.annotate(annotation, (highlight_x, highlight_y), (highlight_x - 1000, highlight_y+0.7), arrowprops=props, ha='right') 
    #TODO ^ make arrow go down if gene is very high on plot
    g.legend([], [], frameon=False)
    g.set_title(title)
    g.set_xlabel("ranked genes")
    g.set_ylabel("log1p(mean expression)")
    return df_new, sum
    

def plot_expressions(GOI, cell_type='T cell', show_summary=False, adata=adata):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    df_allcells, sum_all = plot_expr_class(GOI, ax=axes[0])
    df_celltype, sum_cell = plot_expr_class(GOI, ax=axes[1], cell_type=cell_type)
    fig.tight_layout(pad=2.0)
    #fig.legend(['very_low', 'low', 'middle', 'high', 'very_high'], loc='center left', bbox_to_anchor=(0.5, 1.1), ncol=5) (doesn't work correctly atm)
    plt.show()
    if show_summary:
        print("Summary for all cells: \n" + sum_all)
        print("\nSummary for " + cell_type + " cells: \n" + sum_cell)
    