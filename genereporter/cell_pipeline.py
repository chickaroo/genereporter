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
from typing import Tuple, Union, List
import warnings


class CellPipeline:
    def __init__(self, wdir: str, data_file: str):
        """
        Initialize the CellPipeline class. Load the AnnData file and set the working directory.
        The AnnData given here will be used for all subsequent analyses.

        :param wdir: The working directory.
        :type wdir: str
        :param data_file: The directory to the AnnData file to load.
        :type data_file: str
        """
        warnings.filterwarnings('ignore')
        plt.rcParams["font.family"] = "monospace"
        plt.rcParams["font.size"] = 10
        plt.rcParams['figure.figsize'] = (8,6)
        self.wdir = wdir
        os.chdir(self.wdir)
        self.data_file = data_file
        self.adata = anndata.read_h5ad(self.data_file)


    def get_adata(self) -> anndata.AnnData:
        """
        This function returns the AnnData object loaded in this module.

        :return: The AnnData object loaded in this module.
        :rtype: anndata.AnnData
        """
        return self.adata


    def plot_umap(self, color: str = "celltypist_cell_label_coarse") -> None:
        """
        This function plots a UMAP of the given AnnData object, grouping by coarse cell type in the adata object.

        :param color: The color to use for the plot. Default is "celltypist_cell_label_coarse".
        :type color: str
        """
        sc.pl.umap(
            self.adata, color=[color],
            frameon=False, legend_loc="on data", title=color, legend_fontsize=9, legend_fontweight=600, 
            legend_fontoutline=1
        )


    # Classify gene expression levels in discrete categories
        
    def clean_data(self, df: pd.DataFrame, col: str, threshold = 99.75) -> np.ndarray:
        """
        This function cleans the data by removing outliers.

        :param df: The DataFrame containing the data to be cleaned.
        :type df: pandas.DataFrame
        :param col: The column in the DataFrame to be cleaned.
        :type col: str
        :param threshold: The percentile above which data points are considered outliers. Default is 99.75.
        :type threshold: float, optional
        :return: The cleaned data as a 1D numpy array.
        :rtype: numpy.ndarray
        """
        data = df[col].to_numpy().flatten()  # Convert the column to a 1D numpy array
        #print(type(data))
        data = data[data <= np.percentile(data, threshold)]  # Remove data points above the threshold percentile
        return data  # Return the cleaned data


    def find_thresholds(self, filtered):
        """
        This function calculates and returns four thresholds (very low, low, middle, high) based on the range of the input data.
        The very_high threshold is set to the 99.75th percentile of the input data.

        :param filtered: The input data for which to calculate the thresholds.
        :type filtered: list or numpy.ndarray
        :return: The calculated thresholds (very low, low, middle, high).
        :rtype: tuple
        """
        range = (abs(max(filtered)) + abs(min(filtered)))
        step = range / float(4)
        very_low = min(filtered) + step
        low = min(filtered) + step*2
        middle = min(filtered) + step*3
        return very_low, low, middle, max(filtered)


    def classify_exp_level(self, df, filtered, col) -> Tuple[pd.DataFrame, str]:
        """
        This function classifies the expression level of genes into five categories: very low, low, middle, high, very high.
        It also generates a summary of the quantile thresholds and the number of genes in each category.

        :param df: The input DataFrame containing the gene expression data.
        :type df: pandas.DataFrame
        :param filtered: The filtered gene expression data for which to calculate the thresholds.
        :type filtered: list or numpy.ndarray
        :param col: The column in df containing the gene expression data.
        :type col: str
        :return: The DataFrame with an additional column for the expression level category, and the summary string.
        :rtype: tuple(pandas.DataFrame, str)
        """
        very_low, low, middle, high = self.find_thresholds(filtered)

        def classify(x):
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
            
        df['expr_class'] = df[col].apply(classify)

        summary = str(f"Quantile thresholds: \nvery low: {round(stats.percentileofscore(df[col], very_low), 4)}, low: {round(stats.percentileofscore(df[col], low), 4)}, middle: {round(stats.percentileofscore(df[col], middle), 4)}, high: {round(stats.percentileofscore(df[col], high), 4)}, very high: 99.7500\n")
        summary += "\nNumber of genes per category: \n"
        categories = ['very_low', 'low', 'middle', 'high', 'very_high']
        for category in categories:
            summary += f"{category}: {len (df[df['expr_class'] == category])}\n"
        return df, summary


    def make_df(self, adata, threshold: float=99.75, col: str='log1p(means)', layer: str='log_norm') -> Tuple[pd.DataFrame, str]:
        """
        This function creates a DataFrame from the given AnnData object, where the genes' expression level
        is calculated and classified into five categories: very low, low, middle, high, very high.
        It also generates a summary of the quantile thresholds and the number of genes in each category.

        :param adata: The input AnnData object containing the gene expression data.
        :type adata: anndata.AnnData
        :param threshold: The percentile above which data points are considered outliers (i.e. very_high class). Default is 99.75.
        :type threshold: float, optional
        :param col: The column in the DataFrame to be used for the expression level classification. Default is 'log1p(means)'.
        :type col: str, optional
        :param layer: The layer of the AnnData object to be used for calculating highly variable genes. Default is 'log_norm'.
        :type layer: str, optional
        :return: The DataFrame with additional columns for the gene number and expression level category, and the summary string.
        :rtype: tuple(pandas.DataFrame, str)
        """
        sc.pp.highly_variable_genes(adata, layer=layer)
        df = adata.var.sort_values(['means'])
        df['gene_num'] = range(len(df))
        df['log(means)'] = np.log(df['means'])
        df['log1p(means)'] = np.log1p(df['means'])

        filtered = self.clean_data(df, col=col, threshold=threshold) # can choose upper outlier threshold here
        df_new, summary = self.classify_exp_level(df =df, filtered=filtered, col=col)
        return df_new, summary


    def explain_expr_celltypes(self, GOI: str, col: str='log1p(means)', layer='log_norm')-> pd.DataFrame:
        """
        This function explains the expression of cell types for a given gene of interest (GOI).
        It creates a DataFrame for each cell type, calculates the expression level of the GOI, and concatenates the results.
        The resulting DataFrame is then formatted for output.

        :param GOI: The gene of interest.
        :type GOI: str
        :param adata: The input AnnData object containing the gene expression data.
        :type adata: anndata.AnnData
        :param col: The column in the DataFrame to be used for the expression level calculation. Default is 'log1p(means)'.
        :type col: str, optional
        :param layer: The layer of the AnnData object to be used for calculating highly variable genes. Default is 'log_norm'.
        :type layer: str, optional
        :return: The DataFrame with the expression level of the GOI for each cell type.
        :rtype: pandas.DataFrame
        """
        out = pd.DataFrame()
        
        # Loop over each unique cell type
        for cell_type in self.adata.obs['celltypist_cell_label_coarse'].unique():
            subset = self.adata[self.adata.obs['celltypist_cell_label_coarse'] == cell_type]
            df, sum = self.make_df(adata=subset)
            out = pd.concat([out, df.loc[df.index == GOI]])
            
        # format the output
        out['cell_type'] = self.adata.obs['celltypist_cell_label_coarse'].unique()
        out = out[['cell_type', 'expr_class', col]]
        out[col] = out[col].apply(lambda x: round(x, 3))
        out['expr_class'] = out['expr_class'].apply(lambda x: x.replace('_', ' '))
        out = out.rename(columns={'cell_type': 'Cell type', 'expr_class': 'Expression class', col: 'Avg. expression over cell type'})
        return out


    # TODO threshold for expecting dropout zeros 
    #def dropout_threshold(adata):
        #df = adata.to_df()
        # finalized threshold to come here
        #return df


    # Plot overall expression of all genes, highlighting GOI


    def plot_expr_class(self, GOI: str, ax, adata, cell_type: str = None, col: str = 'log1p(means)') -> Tuple[pd.DataFrame, str]:
        """
        This function plots the expression class of a given gene of interest (GOI) across all cell types or a specific cell type.
        It creates a DataFrame, highlights the GOI on the plot, and annotates it with its expression class.

        :param GOI: The gene of interest.
        :type GOI: str
        :param ax: The axes object to draw the plot onto.
        :type ax: matplotlib.axes.Axes
        :param adata: The input AnnData object containing the gene expression data. Default is the global adata object.
        :type adata: anndata.AnnData, optional
        :param cell_type: The specific cell type to plot. If None, all cell types are plotted. Default is None.
        :type cell_type: str, optional
        :param col: The column in the DataFrame to be used for the y-axis. Default is 'log1p(means)'.
        :type col: str, optional
        :return: The DataFrame with the expression level of the GOI for each cell type, and the summary string.
        :rtype: tuple(pandas.DataFrame, str)
        """
        title = "Overall expression across all cell types"
        # If a specific cell type is provided, subset the data for that cell type
        if cell_type != None:
            title = str("Overall expression in " + cell_type + " cells")
            if cell_type in adata.obs['celltypist_cell_label_coarse'].unique():
                adata = adata[adata.obs['celltypist_cell_label_coarse'] == cell_type]
            elif cell_type in adata.obs['celltypist_cell_label'].unique():
                adata = adata[adata.obs['celltypist_cell_label'] == cell_type]
            else:
                print("Cell type not found. Please check spelling.")
                return
        
        # Create a DataFrame and calculate the expression level of the GOI
        df_new, sum = self.make_df(adata=adata)
        # Plot the data with seaborn
        g = sns.scatterplot(data=df_new, ax=ax, x='gene_num',  y=col, hue='expr_class', linewidth=0)
        # Highlight the GOI on the plot and annotate it with its expression class
        annotation = str(GOI + " (" + df_new.loc[df_new.index == GOI]['expr_class'].values[0] + ")")
        highlight_y = df_new.loc[df_new.index == GOI][col]
        highlight_x = df_new.loc[df_new.index == GOI]['gene_num']
        g.scatter(highlight_x, highlight_y, color = 'yellow', linewidth=1)
        props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
        g.annotate(annotation, (highlight_x, highlight_y), (highlight_x - 1000, highlight_y+0.7), arrowprops=props, ha='right')
        # Set the title and labels of the plot
        g.set_title(title)
        g.set_xlabel("ranked genes")
        g.set_ylabel("log1p(mean expression)")
        return df_new, sum


    def plot_expressions(self, GOI: str, cell_type: str = 'T cell', show_summary: bool = False) -> None:
        """
        This function plots the expression class of a given gene of interest (GOI) across all cell types and a specific cell type.
        It creates two subplots, one for all cell types and one for the specific cell type, and optionally prints a summary of the expression class for each plot.

        :param GOI: The gene of interest.
        :type GOI: str
        :param cell_type: The specific cell type to plot. Default is 'T cell'.
        :type cell_type: str, optional
        :param show_summary: Whether to print a summary of the expression class for each plot. Default is False.
        :type show_summary: bool, optional
        :param adata: The input AnnData object containing the gene expression data. Default is the global adata object.
        :type adata: anndata.AnnData, optional
        :return: None
        """
        # Create a figure with two subplots
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
        # Plot the expression class 
        df_allcells, sum_all = self.plot_expr_class(GOI, ax=axes[0], adata=self.adata)
        df_celltype, sum_cell = self.plot_expr_class(GOI, ax=axes[1], adata=self.adata, cell_type=cell_type)
        # Adjust and show
        fig.tight_layout(pad=2.0)
        plt.show()
        # If show_summary is True, print a summary of the expression class for each plot
        if show_summary:
            print("Summary for all cells: \n" + sum_all)
            print("\nSummary for " + cell_type + " cells: \n" + sum_cell)
        

    # standard scanpy GOI expression plotting 
            
    def dotplot(self, GOI: str, cell_type: str = None, **kwargs) -> None:
        """
        This function creates a dot plot of the expression of a given gene of interest (GOI) across all cell types or a specific cell type.

        :param GOI: The gene of interest.
        :type GOI: str
        :param cell_type: The specific cell type to plot. If None, all cell types are plotted. Default is None.
        :type cell_type: str, optional
        :param adata: The input AnnData object containing the gene expression data. Default is the global adata object.
        :type adata: anndata.AnnData, optional
        :param kwargs: Additional keyword arguments to be passed to the sc.pl.dotplot function.
        :return: None
        """
        if cell_type == None:
            fig = sc.pl.dotplot(self.adata, var_names=GOI, groupby='celltypist_cell_label_coarse', 
                            return_fig=True, standard_scale='var', title="Coarse Cell Types", figsize=(3.5,8), **kwargs)
            fig.add_totals().show()
        else:
            fig = sc.pl.dotplot(self.adata[self.adata.obs['celltypist_cell_label_coarse'] == cell_type], 
                                var_names=GOI, groupby='celltypist_cell_label', return_fig=True,
                                standard_scale='var', title=str('Fine Cell Types: ' + cell_type), figsize=(4,3), **kwargs)
            fig.add_totals().show()

        plt.show()

    def matrixplot(self, GOI: str, cell_type: str = None) -> None:
        """
        This function creates a matrix plot of the expression of a given gene of interest (GOI) across all cell types or a specific cell type.

        :param GOI: The gene of interest.
        :type GOI: str
        :param cell_type: The specific cell type to plot. If None, all cell types are plotted. Default is None.
        :type cell_type: str, optional
        :param adata: The input AnnData object containing the gene expression data. Default is the global adata object.
        :type adata: anndata.AnnData, optional
        :return: None
        """
        if cell_type == None:
            fig = sc.pl.matrixplot(self.adata, var_names=GOI, groupby='celltypist_cell_label_coarse', layer='log_norm',
                            return_fig=True, standard_scale='var', title="Coarse Cell Types", figsize=(3.5,8))
            fig.add_totals().style(edge_color='black').show()
        else:
            fig = sc.pl.matrixplot(self.adata[self.adata.obs['celltypist_cell_label_coarse'] == cell_type], 
                                var_names=GOI, groupby='celltypist_cell_label', return_fig=True, layer='log_norm',
                                standard_scale='var', title=str('Fine Cell Types: ' + cell_type), figsize=(4,3))
            fig.add_totals().style(edge_color='black').show()
        plt.show()

    def heatmap(self, GOI: str, layer: str = 'log_norm') -> None:
        """
        This function creates a heatmap of the expression of a given gene of interest (GOI) across all cell types.

        :param GOI: The gene of interest.
        :type GOI: str
        :param layer: The layer of the AnnData object to be used for calculating highly variable genes. Default is 'log_norm'.
        :type layer: str, optional
        :param adata: The input AnnData object containing the gene expression data. Default is the global adata object.
        :type adata: anndata.AnnData, optional
        :return: None
        """
        sc.pl.heatmap(self.adata, GOI, groupby='celltypist_cell_label_coarse', swap_axes=True, figsize=(18,1.5), layer=layer, standard_scale='var')


    ### Expression vs. detection calculations and plots ###
        

    # for cell type: mean expression of a gene (x) vs. percentace of cells where this gene is detected (y) (wihtin a cell type)
    def expression_vs_detection(self, GOI:str, adata = None, cell_type: str=None, layer:str='log_norm', col:str='log1p(means)', return_df:bool = False)-> Union[None, pd.DataFrame]:
        """
        This function creates a DataFrame, calculates the percentage of cells where each gene is detected, and optionally returns the DataFrame.
        If return_df=False, the function plots the mean expression of a given gene of interest (GOI) versus the percentage of cells where this gene is detected, either across all cell types or a specific cell type.

        :param GOI: The gene of interest.
        :type GOI: str
        :param adata: The input AnnData object containing the gene expression data. Default is the global adata object.
        :type adata: anndata.AnnData, optional
        :param cell_type: The specific cell type to plot. If None, all cell types are plotted. Default is None.
        :type cell_type: str, optional
        :param layer: The layer of the AnnData object to be used for calculating highly variable genes. Default is 'log_norm'.
        :type layer: str, optional
        :param col: The column in the DataFrame to be used for the x-axis. Default is 'log1p(means)'.
        :type col: str, optional
        :param return_df: Whether to return the DataFrame. If False, the function will plot the data and return None. Default is False.
        :type return_df: bool, optional
        :return: If return_df is True, the DataFrame with the mean expression and percentage detection of the GOI for each cell type. Otherwise, None.
        :rtype: Union[None, pandas.DataFrame]
        """
        if adata == None:
            adata = self.adata
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
            
        df, sum = self.make_df(adata, col=col, layer=layer)

        # calculate percentage of cells (of the given cell type) where each gene is detected
        subset_df = adata.to_df(layer=layer)
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

    def fit_spline(self, cell_type: str = None, plot: bool = False) -> Union[None, Tuple[np.ndarray, np.ndarray]]:
        """
        This function fits a spline to the expression versus detection data of a given cell type, and optionally plots the data and the fitted spline.

        :param cell_type: The specific cell type to plot. If None, all cell types are plotted. Default is None.
        :type cell_type: str, optional
        :param adata: The input AnnData object containing the gene expression data. Default is the global adata object.
        :type adata: anndata.AnnData, optional
        :param plot: Whether to plot the data and the fitted spline. If False, the function will return the spline parameters and inflection points. Default is False.
        :type plot: bool, optional
        :return: If plot is False, a tuple containing the spline parameters and the inflection points. Otherwise, None.
        :rtype: Union[None, Tuple[numpy.ndarray, numpy.ndarray]]
        """
        # make bins and take maximum of each bin to get (x,y) points to fit sigmoid curve to
        df = self.expression_vs_detection(self.adata, cell_type=cell_type, return_df=True)
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
    def calc_distance_point(self, point: Tuple[float, float], p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
        """
        This function calculates the orthogonal distance of a point from the line defined by two other points.

        :param point: The point for which the distance is to be calculated.
        :type point: Tuple[float, float]
        :param p1: The first point defining the line.
        :type p1: Tuple[float, float]
        :param p2: The second point defining the line.
        :type p2: Tuple[float, float]
        :return: The orthogonal distance of the point from the line.
        :rtype: float
        """
        x0, y0 = point
        x1, y1 = p1
        x2, y2 = p2
        return abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) / math.sqrt((y2-y1)**2 + (x2-x1)**2) # orthogonal distance formula

    def calc_distance(self, x: pd.Series, y: pd.Series, cell_type: str = None) -> List[float]:
        """
        This function calculates the orthogonal distance of each point in the given x and y series from the spline fitted to the data.

        :param x: The x-coordinates of the points.
        :type x: pandas.Series
        :param y: The y-coordinates of the points.
        :type y: pandas.Series
        :return: A list of the orthogonal distances of each point from the spline.
        :rtype: List[float]
        """
        # make  linear interpolations of the fit function, eg. from 0 to infl point and from infl point to max
        tck, inflection = self.fit_spline(cell_type=cell_type)

        lowest = [0,0]
        highest = [x[-1], BSpline(*tck)(x[-1])]
        
        distance = []
        complete_data = pd.DataFrame(data={'x': x.values, 'y': y.values})
        segment = complete_data[complete_data['x'] < inflection[0]]


        for i, j in zip(segment['x'], segment['y']):
            distance.append(self.calc_distance_point([i,j], lowest, [inflection[0], BSpline(*tck)(inflection[0])]))

        for infls in range(len(inflection)-1):
            segment = complete_data[(complete_data['x'] >= inflection[infls]) & (complete_data['x'] < inflection[infls+1])]
            for i, j in zip(segment['x'], segment['y']):
                distance.append(self.calc_distance_point([i,j], [inflection[infls], BSpline(*tck)(inflection[infls])], [inflection[infls+1], BSpline(*tck)(inflection[infls+1])]))
        
        segment = complete_data[(complete_data['x'] >= inflection[-1])]
        for i, j in zip(segment['x'], segment['y']):
            distance.append(self.calc_distance_point([i,j], [inflection[-1], BSpline(*tck)(inflection[-1])], highest))

        return distance    

    def detect_outliers(self, cell_type: str = None, outlier_threshold: float = 0.15) -> pd.DataFrame:
        """
        This function detects outliers in the expression versus detection data of a given cell type. 
        It calculates the orthogonal distance of each point from the spline fitted to the data, 
        and marks points with a distance greater than the specified outlier threshold as outliers.

        :param cell_type: The specific cell type to plot. If None, all cell types are plotted. Default is None.
        :type cell_type: str, optional
        :param outlier_threshold: The threshold for marking a point as an outlier. If the orthogonal distance of a point from the spline is greater than this threshold, it is marked as an outlier. Default is 0.15.
        :type outlier_threshold: float, optional
        :param adata: The input AnnData object containing the gene expression data. Default is the global adata object.
        :type adata: anndata.AnnData, optional
        :return: A DataFrame with the mean expression, percentage detection, orthogonal distance from the spline, and outlier status of each gene for the given cell type.
        :rtype: pandas.DataFrame
        """    
        df = self.expression_vs_detection(self.adata, cell_type=cell_type, return_df=True)
        df['percent_detected'] = df['percent_detected']/100

        detecter = df[['log1p(means)', 'percent_detected']]
        detecter['distance'] = self.calc_distance(detecter['log1p(means)'], detecter['percent_detected'], cell_type=cell_type)
        detecter['is_outlier'] = detecter['distance'] > outlier_threshold
        return detecter


    def plot_outliers(self, GOI: str, cell_type: str = None, outlier_threshold: float=0.15) -> None:
        """
        This function plots the expression versus detection data for a given cell type, with outliers highlighted. 
        It also highlights a gene of interest (GOI) in the plot.

        :param GOI: The gene of interest to be highlighted in the plot.
        :type GOI: str
        :param cell_type: The specific cell type to plot. If None, all cell types are plotted. Default is None.
        :type cell_type: str, optional
        :param outlier_threshold: The threshold for marking a point as an outlier. If the orthogonal distance of a point from the spline is greater than this threshold, it is marked as an outlier. Default is 0.15.
        :type outlier_threshold: float, optional
        :return: None
        """
        # make plot with outliers highlighted
        detecter = self.detect_outliers(cell_type=cell_type, outlier_threshold=outlier_threshold)
        tck, infls = self.fit_spline(cell_type=cell_type)
        mask = detecter['log1p(means)'].isin(infls)
        inflections = detecter[mask]
        inflections['spline'] = BSpline(*tck)(inflections['log1p(means)'])
        inflections['log1p(means)'][0]
        if cell_type == None:
            title = "expression vs. detected; outliers highlighted"
        else:
            title = str(cell_type + " expression vs. detected; outliers highlighted")

        fig, ax = plt.subplots(1,1)
        sns.scatterplot(data=detecter, x='log1p(means)', y='percent_detected', hue='is_outlier', linewidth=0, alpha=1, ax=ax)
        x_list = [[0], [x for x in inflections['log1p(means)']], [detecter['log1p(means)'][-1]]]
        y_list = [[0], [y for y in inflections['spline']], [BSpline(*tck)(detecter['log1p(means)'][-1])]]
        x_list = list(chain(*x_list))
        y_list = list(chain(*y_list))
        ax.plot(x_list, y_list, color='red', linestyle='dashed')

        # highlight GOI 
        annotation = str(GOI)
        highlight_y = detecter.loc[detecter.index == GOI]['percent_detected']
        highlight_x = detecter.loc[detecter.index == GOI]['log1p(means)']
        ax.scatter(highlight_x, highlight_y, color = 'yellow', linewidth=1, alpha=0.8)
        props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
        ax.annotate(annotation, (highlight_x, highlight_y), (highlight_x-0.25, highlight_y+0.05), arrowprops=props)

        ax.set_title(title)
        ax.legend(title='Outlier', loc='lower left', bbox_to_anchor=(1, 0))
        plt.show()

    def list_outliers(self, cell_type: str = None, head: int = 5) -> pd.DataFrame:
        """
        This function lists the top outliers in the expression versus detection data for a given cell type. 
        It sorts the data by the orthogonal distance from the spline, in descending order, and returns the top 'head' number of outliers.

        :param cell_type: The specific cell type to plot. If None, all cell types are plotted. Default is None.
        :type cell_type: str, optional
        :param head: The number of top outliers to return. Default is 5.
        :type head: int, optional
        :return: A DataFrame with the top 'head' number of outliers, sorted by the orthogonal distance from the spline in descending order.
        :rtype: pandas.DataFrame
        """
        detecter = self.detect_outliers(cell_type=cell_type)
        detecter = detecter.sort_values(by='distance', ascending=False)
        return detecter[detecter['is_outlier'] == True].head(head)