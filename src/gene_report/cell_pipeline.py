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
from itertools import chain
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
        df, sum = make_df(subset, col)
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
    

# standard scanpy GOI expression plotting 
        
def dotplot(GOI, cell_type=None, adata=adata):
    if cell_type == None:
        fig = sc.pl.dotplot(adata, var_names=GOI, groupby='celltypist_cell_label_coarse', 
                        return_fig=True, standard_scale='var', title="Coarse Cell Types", figsize=(3.5,8))
        fig.add_totals().show()
    else:
        fig = sc.pl.dotplot(adata[adata.obs['celltypist_cell_label_coarse'] == cell_type], 
                            var_names=GOI, groupby='celltypist_cell_label', return_fig=True,
                            standard_scale='var', title=str('Fine Cell Types: ' + cell_type), figsize=(4,3))
        fig.add_totals().show()

    plt.show()

def matrixplot(GOI, cell_type=None, adata=adata):
    if cell_type == None:
        fig = sc.pl.matrixplot(adata, var_names=GOI, groupby='celltypist_cell_label_coarse', layer='log_norm',
                        return_fig=True, standard_scale='var', title="Coarse Cell Types", figsize=(3.5,8))
        fig.add_totals().style(edge_color='black').show()
    else:
        fig = sc.pl.matrixplot(adata[adata.obs['celltypist_cell_label_coarse'] == cell_type], 
                            var_names=GOI, groupby='celltypist_cell_label', return_fig=True, layer='log_norm',
                            standard_scale='var', title=str('Fine Cell Types: ' + cell_type), figsize=(4,3))
        fig.add_totals().style(edge_color='black').show()
    plt.show()

def heatmap(GOI, layer='log_norm', adata=adata):
    sc.pl.heatmap(adata, GOI, groupby='celltypist_cell_label_coarse', swap_axes=True, figsize=(18,1.5), layer=layer, standard_scale='var')


# Expression vs. detection calculations and plots

# for cell type: mean expression of a gene (x) vs. percentace of cells where this gene is detected (y) (wihtin a cell type)

def expression_vs_detection(GOI, adata=adata, cell_type=None, layer='log_norm', col='log1p(means)', ax=None, return_df = False):
    
    title=str(GOI)
    if cell_type != None:
        title = str("Overall expression in " + cell_type + " cells")
        if cell_type in adata.obs['celltypist_cell_label_coarse'].unique():
            adata = adata[adata.obs['celltypist_cell_label_coarse'] == cell_type]
        elif cell_type in adata.obs['celltypist_cell_label'].unique():
            adata = adata[adata.obs['celltypist_cell_label'] == cell_type]
        else:
            print("Cell type not found. Please check spelling.")
            return
    else:
        title = "All cell types"
        
    df, sum = make_df(adata, col, layer=layer)

    # calculate percentage of cells (of the given cell type) where each gene is detected
    subset_df = adata.to_df(layer=layer)
    #nonzero_detected = pd.DataFrame(subset_df.astype(bool).sum(axis=0) / len(df) * 100, columns=['percent_detected'])
    nonzero_detected = pd.DataFrame(np.count_nonzero(subset_df, axis=0) / len(subset_df)*100, columns=['percent_detected'], index=subset_df.columns)
    df = df.join(nonzero_detected, how='left').sort_values(['gene_num'], ascending=True)
    if return_df:
        return df
        
    # plot mean expression of a gene (x) vs. percentage of cells where this gene is detected (y)
    ax = sns.scatterplot(data=df, x=col, y='percent_detected', hue='expr_class', linewidth=0)
    ax.set_title(str(title+ ": mean expression vs. percentage detected"))
    ax.legend(title='Expression Class', loc='lower left', bbox_to_anchor=(1, 0))

    annotation = str(GOI)
    highlight_y = df.loc[df.index == GOI]['percent_detected']
    highlight_x = df.loc[df.index == GOI][col]
    ax.scatter(highlight_x, highlight_y, color = 'yellow', linewidth=1)
    props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
    ax.annotate(annotation, (highlight_x, highlight_y), (highlight_x-0.1, highlight_y+20), arrowprops=props)
    plt.show()
    #return ax

def fit_spline(cell_type=None, adata=adata, plot=False):
    # make bins and take maximum of each bin to get (x,y) points to fit sigmoid curve to
    df = expression_vs_detection(adata, cell_type=cell_type, return_df=True)
    df['percent_detected'] = df['percent_detected']/100
    df_bins = df.groupby(pd.cut(df['log1p(means)'], np.arange(0, max(df['log1p(means)']), 0.05))).max() 
    df_bins = df_bins.dropna()
    # fit function here
    xdata = df_bins['log1p(means)'].values
    ydata = df_bins['percent_detected']
    for bin, value in ydata.items():
        if (bin.right > 1.25) and (value < 0.7):
            ydata.loc[bin] = 0.999
    ydata = ydata.values

    # fit and plot splines
    tck = splrep(xdata, ydata, s=0.0009, k=3) #TODO find optimal s value??? 
    #tck_s = splrep(xdata, ydata, s=len(xdata))

    # derivatives of spline function
    yders = interpolate.spalde(xdata, tck)
    yders_df = pd.DataFrame(yders, columns=['y', 'dy', 'd2y', 'd3y'], index= xdata)
    yders_df['d3y_diff'] = yders_df['d3y'].diff().fillna(0)
    # identify any changes in direction of the third derivative
    infls = yders_df[yders_df['d3y_diff'] != 0].index

    if plot:
        # plot data vs fitted spline with inflection points
        fig, ax2 = plt.subplots(1,1)
        ax2.scatter(data=df, x='log1p(means)', y='percent_detected', alpha=0.5, color='orange', s=2.5)
        ax2.plot(xdata, BSpline(*tck)(xdata), label='s=0.001')
        #ax2.plot(xdata, BSpline(*tck_s)(xdata), label=f's={len(xdata)}')
        ax2.set_title(str(cell_type + " Spline fit to expression vs. fraction detected"))
        ax2.set_xlabel("log1p(means)")
        ax2.set_ylabel("fraction detected")
        ax2.vlines(infls, 0, 1.05, color='red', linestyles='dashed', label='turning points')
        ax2.legend()
        plt.show()
    else:
        return tck, infls

# calculate the orthogonal distance of each point (gene) from the spline
# if greater than a certain value (threshold), then it is an outlier

# calculate the distance of each point from the linear interpolation of its respective section
def calc_distance_point(point, p1, p2):
    x0, y0 = point
    x1, y1 = p1
    x2, y2 = p2
    return abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) / math.sqrt((y2-y1)**2 + (x2-x1)**2) # orthogonal distance formula

def calc_distance(x, y):
    # make  linear interpolations of the sigmoid function, eg. one from 0 to inflection point, one from inflection point to max

    tck, inflection = fit_spline()

    lowest = [0,0]
    highest = [x[-1], BSpline(*tck)(x[-1])]
    
    distance = []
    complete_data = pd.DataFrame(data={'x': x.values, 'y': y.values})
    segment = complete_data[complete_data['x'] < inflection[0]]


    for i, j in zip(segment['x'], segment['y']):
        distance.append(calc_distance_point([i,j], lowest, [inflection[0], BSpline(*tck)(inflection[0])]))

    for infls in range(len(inflection)-1):
        segment = complete_data[(complete_data['x'] >= inflection[infls]) & (complete_data['x'] < inflection[infls+1])]
        for i, j in zip(segment['x'], segment['y']):
            distance.append(calc_distance_point([i,j], [inflection[infls], BSpline(*tck)(inflection[infls])], [inflection[infls+1], BSpline(*tck)(inflection[infls+1])]))
    
    segment = complete_data[(complete_data['x'] >= inflection[-1])]
    for i, j in zip(segment['x'], segment['y']):
        distance.append(calc_distance_point([i,j], [inflection[-1], BSpline(*tck)(inflection[-1])], highest))

    return distance    

def detect_outliers(cell_type=None):
    df = expression_vs_detection(adata, cell_type=cell_type, return_df=True)
    df['percent_detected'] = df['percent_detected']/100

    detecter = df[['log1p(means)', 'percent_detected']]
    detecter['distance'] = calc_distance(detecter['log1p(means)'], detecter['percent_detected'])
    detecter['is_outlier'] = detecter['distance'] > 0.15
    return detecter


def plot_outliers(GOI, cell_type=None):
    # make plot with outliers highlighted
    detecter = detect_outliers(cell_type=cell_type)
    tck, infls = fit_spline(cell_type=cell_type)
    mask = detecter['log1p(means)'].isin(infls)
    inflections = detecter[mask]
    inflections['spline'] = BSpline(*tck)(inflections['log1p(means)'])
    inflections['log1p(means)'][0]
    if cell_type == None:
        title = "expression vs. detected; outliers highlighted"
    else:
        title = str(cell_type + " expression vs. detected; outliers highlighted")

    fig, ax = plt.subplots(1,1)
    sns.scatterplot(data=detecter, x='log1p(means)', y='percent_detected', hue='is_outlier', alpha=1, ax=ax)
    x_list = [[0], [x for x in inflections['log1p(means)']], [detecter['log1p(means)'][-1]]]
    y_list = [[0], [y for y in inflections['spline']], [BSpline(*tck)(detecter['log1p(means)'][-1])]]
    x_list = list(chain(*x_list))
    y_list = list(chain(*y_list))
    ax.plot(x_list, y_list, color='red', linestyle='dashed')

    # highlight GOI 
    annotation = str(GOI)
    highlight_y = detecter.loc[detecter.index == GOI]['percent_detected']
    highlight_x = detecter.loc[detecter.index == GOI]['log1p(means)']
    ax.scatter(highlight_x, highlight_y, color = 'yellow', linewidth=1)
    props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
    ax.annotate(annotation, (highlight_x, highlight_y), (highlight_x-0.25, highlight_y+0.05), arrowprops=props)

    ax.set_title(title)
    ax.legend(title='Outlier; thresh 0.15', loc='lower left', bbox_to_anchor=(1, 0))
    plt.show()

def list_outliers(cell_type=None, head=5):
    detecter = detect_outliers(cell_type=cell_type)
    detecter = detecter.sort_values(by='distance', ascending=False)
    return detecter[detecter['is_outlier'] == True].head(head)