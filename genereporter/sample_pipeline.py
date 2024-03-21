import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import warnings

class SamplePipeline():
    def __init__(self, wdir, adata):
        """
        Initialize the SamplePipeline class. This class is used to generate plots for the GOI across different samples (patients).

        :param wdir: The working directory.
        :type wdir: str
        :param adata: The AnnData object.
        :type adata: AnnData
        """
        warnings.filterwarnings('ignore')
        pd.options.mode.chained_assignment = None  # default='warn'
        # set a working directory
        os.chdir( wdir )
        self.adata = sc.read_h5ad(adata)

    def get_adata(self):
        return self.adata

    def sample_v_celltype(self, GOI) -> pd.DataFrame:
        """
        Get expression data for a gene of interest (GOI), calculate mean expression per celltype and sample,
        and count the number of cells for each cell type.

        :param GOI: The gene of interest.
        :type GOI: str
        :return: DataFrame with mean expression per celltype and sample, and cell counts.
        :rtype: DataFrame
        """
        # get data for GOI
        goi_adata = self.adata[:, GOI]
        goi_df = goi_adata.to_df(layer='log_norm')
        goi_df['celltype'] = goi_adata.obs['celltypist_cell_label_coarse']
        goi_df['sample'] = goi_adata.obs['sampleID']
        # calculate mean expression per celltype and sample
        goi_df_mean = goi_df.groupby(['celltype', 'sample'], observed=True).mean()
        goi_df_mean = goi_df_mean.pivot_table(index='celltype', columns='sample', values=GOI)
        goi_df_mean = goi_df_mean.fillna(0)
        # count the number of cells for each cell type and add as a new column
        cell_counts = goi_df['celltype'].value_counts()
        goi_df_mean['cell_count'] = goi_df_mean.index.map(cell_counts.to_dict())

        return goi_df_mean
    
    def pl_sample_celltype(self, GOI):
        """
        Generate a clustermap of GOI expression data across patients and cell types. 
        The clustermap is accompanied by a barplot showing the number of cells per cell type.

        :param GOI: The gene of interest.
        :type GOI: str
        """

        # Assuming goi_expr is your DataFrame with the gene expression data
        # and it has a column 'cell_count' with the cell counts
        goi_expr = self.sample_v_celltype(GOI)
        counts_df = goi_expr['cell_count']
        goi_expr = goi_expr.drop(columns='cell_count')

        # Generate the clustermap
        g = sns.clustermap(goi_expr, cmap='YlGnBu', cbar_kws={'label': 'mean expression in group'}, 
                        dendrogram_ratio=(0, 0.08), standard_scale=1, figsize=(9,6))

        # Rotate x-axis labels
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right') 

        # Add title
        g.ax_heatmap.set_title(f'{GOI} expression per sample per celltype', y=1.1)
        # move y axis to the left
        g.ax_heatmap.yaxis.set_label_position("left")
        g.ax_heatmap.yaxis.tick_left()

        g.figure.subplots_adjust(right=0.8)
        total_barplot_ax = g.figure.add_axes([0.8, 0.1, 0.25, 0.805])

        # match the order of the barplot with the order of the y-axis of the heatmap
        y_axis = list(reversed([x.get_text() for x in g.ax_heatmap.get_yticklabels()]))
        counts_df = counts_df.reindex(y_axis)

        counts_df.plot(
            kind="barh",
            color="salmon",
            position=-0.2,
            ax=total_barplot_ax,
            edgecolor="black",
            width=0.5,
        )

        # edit the colorbar
        g.cax.set_visible(False)
        cbar_ax = g.figure.add_axes([0.85, 0, 0.1, 0.05])
        g.figure.colorbar(g.ax_heatmap.get_children()[0], cax=cbar_ax, orientation='horizontal', label='mean expression in group')

        # add numbers to the right of the bars
        max_x = max([p.get_width() for p in total_barplot_ax.patches])
        for p in total_barplot_ax.patches:
            if p.get_width() >= 1000:
                display_number = f"{np.round(p.get_width() / 1000, decimals=1)}k"
            else:
                display_number = np.round(p.get_width(), decimals=1)
            total_barplot_ax.annotate(
                display_number,
                ((p.get_width()), p.get_y() + p.get_height()),
                ha="center",
                va="top",
                xytext=(10, 0),
                fontsize="x-small",
                textcoords="offset points",
            )
        total_barplot_ax.set_xlim(0, max_x * 1.4)
        total_barplot_ax.grid(False)
        total_barplot_ax.axis("off")

        plt.show()


    def pl_violin(self, GOI, celltype=None):
        """
        Generate a violin plot of gene of interest expression across all or a specific cell type of interest.

        :param GOI: The gene of interest.
        :type GOI: str
        :param celltype: The cell type of interest, defaults to None.
        :type celltype: str, optional
        """
        # violin plot of gene of interest expression across all or a specific cell type of interest
        temp = self.adata
        title = f"{GOI} expression per sample across all cell types"
        if celltype != None:
            temp = self.adata[self.adata.obs['celltypist_cell_label_coarse'] == celltype]
            title = f'{GOI} expression per sample in {celltype} cells'
        fig, ax1 = plt.subplots(figsize=(10, 5))
        sc.pl.violin(temp, keys=GOI, groupby='sampleID', rotation=90, stripplot=True, jitter=True, show=False, ax=ax1)
        ax1.set_title(title)
        ax1.set_ylabel(f'{GOI} Expression')
        ax1.set_xlabel('Sample ID')
        plt.show()