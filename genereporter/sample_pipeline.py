import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import warnings

class SamplePipeline:
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

        plt.rcParams["font.size"] = 9
        plt.rcParams['figure.figsize'] = (10,6)
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
        goi_df = goi_adata.to_df()
        goi_df['celltype'] = goi_adata.obs['celltype_l2']
        goi_df['sample'] = goi_adata.obs['PatientID_genotype']
        # calculate mean expression per celltype and sample
        goi_df_mean = goi_df.groupby(['celltype', 'sample'], observed=True).mean()
        goi_df_mean = goi_df_mean.pivot_table(index='celltype', columns='sample', values=GOI)
        #goi_df_mean = goi_df_mean.fillna(0)
        goi_df_mean = goi_df_mean.apply(lambda row: row.fillna(row.mean()), axis=1) # fill missing boxes with mean of cell type
        # count the number of cells for each cell type and add as a new column
        cell_counts = goi_df['celltype'].value_counts()
        goi_df_mean['cell_count'] = goi_df_mean.index.map(cell_counts.to_dict())

        return goi_df_mean
    
    def pl_sample_celltype(self, GOI, standard_scale=1, z_score=False):
        """
        Generate a clustermap of GOI expression data across patients and cell types. 
        The clustermap is accompanied by a barplot showing the number of cells per cell type.

        :param GOI: The gene of interest.
        :type GOI: str
        """

        # Assuming goi_expr is your DataFrame with the gene expression data
        # and it has a column 'cell_count' with the cell counts
        goi_expr = self.sample_v_celltype(GOI)
        counts_df = pd.DataFrame(goi_expr['cell_count'])
        goi_expr = goi_expr.drop(columns='cell_count')

        # Generate the clustermap
        if z_score:
            label = "Z-Score in group"
            g = sns.clustermap(goi_expr, mask=False, cmap='BrBG', z_score=1, center=0, 
                        xticklabels=True, dendrogram_ratio=(0, 0.09), figsize=(12,8))
        else:
            label = "Mean expression in group"
            g = sns.clustermap(goi_expr, mask=False, cmap='YlGnBu', standard_scale=1, 
                        xticklabels=True, dendrogram_ratio=(0, 0.09), figsize=(12,8))

        # Rotate x-axis labels
        #plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right') 

        # Add title
        g.ax_heatmap.set_title(f'{GOI} Expression per Sample per Celltype', y=1.1)
        # move y axis to the left
        g.ax_heatmap.yaxis.set_label_position("left")
        g.ax_heatmap.yaxis.tick_left()
        g.ax_heatmap.set_xticklabels((str(x.get_text()).replace("PATIENT", "") for x in g.ax_heatmap.get_xticklabels()), ha='center')
        g.ax_heatmap.set_xlabel('Patient ID')

        g.figure.subplots_adjust(right=0.85)
        total_barplot_ax = g.figure.add_axes([0.85, 0.26, 0.2, 0.66])

        # match the order of the barplot with the order of the y-axis of the heatmap
        y_axis = list(reversed([x.get_text() for x in g.ax_heatmap.get_yticklabels()]))
        counts_df = counts_df.reindex(y_axis)
        counts_df['cell_count'] = counts_df['cell_count'].astype(int)
        #print(counts_df.head())
        #print(counts_df.columns)
        counts_df.plot(
            kind="barh",
            color="salmon",
            position=-0.2,
            ax=total_barplot_ax,
            edgecolor="black",
            width=0.5,
            legend=False
        )

        # edit the colorbar
        g.cax.set_visible(False)
        cbar_ax = g.figure.add_axes([0.9, 0.05, 0.1, 0.05])
        g.figure.colorbar(g.ax_heatmap.get_children()[0], cax=cbar_ax, orientation='horizontal', label=label)

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

        goi_adata = self.adata[:, GOI]
        
        title = f"{GOI} Expression per Sample Across all Cell Types"
        if celltype != None:
            goi_adata = self.adata[self.adata.obs['celltype_l2'] == celltype]
            title = f'{GOI} Expression per Sample in {celltype} Cells'
        
        # sort samples by mean expression in group (all cells or specific cell type)
        goi_df = goi_adata.to_df()
        goi_df['celltype'] = goi_adata.obs['celltype_l2']
        goi_df['sample'] = goi_adata.obs['PatientID_genotype']
        # calculate mean expression per sample
        goi_df_mean = goi_df.groupby(['sample'], observed=True).mean(numeric_only=True)
        goi_df_mean.sort_values(by=GOI, axis=0, ascending=False, inplace=True)
        sorted_patients = list(goi_df_mean.index)

        plt.rcParams["font.size"] = 9
        plt.rcParams['figure.figsize'] = (14,6)
        fig, ax1 = plt.subplots()
        sc.pl.violin(goi_adata, keys=GOI, groupby='PatientID_genotype', rotation=90, stripplot=True, 
                     jitter=False, size=1.5, log=False, linewidth=0, order=sorted_patients, 
                     show=False, ax=ax1)
        ax1.set_title(title)
        ax1.set_ylabel(f'{GOI} Expression')
        ax1.set_xlabel('Patient ID')
        ax1.set_xticklabels((str(x.get_text()).replace("PATIENT", "") for x in ax1.get_xticklabels()), ha='center', fontsize=7)
        plt.show()