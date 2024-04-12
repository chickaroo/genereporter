import warnings
import os
import pandas as pd
from pyvis.network import Network
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import time
import xmltodict
from collections import defaultdict
from Bio import Entrez
from pathlib import Path
import requests
import plotly.express as px
from itertools import chain, repeat
from IPython.display import display, HTML
from typing import List
import warnings

class GRNPipeline():
    def __init__(self, wdir, adata, f_adj, f_reg, dir_gg_adj):
        warnings.filterwarnings('ignore')
        pd.options.mode.chained_assignment = None  # default='warn'

        # set a working directory
        os.chdir( wdir )

        # load the adata object
        adata = sc.read_h5ad(adata)

        # read the adjacency and regulon data
        adjacencies = pd.read_csv(f_adj)
        regulon = pd.read_csv(f_reg)
        def clean_target_genes(row: pd.Series) -> list:
            return eval(row['TargetGenes'])
        regulon.apply(clean_target_genes, axis=1)
        

        # assign to self
        self.adata = adata
        self.adj_df = adjacencies
        self.reg_df = regulon
        self.reactome = self.get_reactome()
        self.regulon_geneset = self.get_regulon_genesets()
        self.geneset_df = self.get_genesets()

        self.bcell_adj = pd.read_csv(os.path.join(dir_gg_adj, 'gg_adj_B_Cell.csv'), names=['gene','target','importance'])
        self.bcell_adj['cell_lineage'] = 'b_cell'
        self.epithelium_adj = pd.read_csv(os.path.join(dir_gg_adj, 'gg_adj_Epithelium.csv'), names=['gene','target','importance'])
        self.epithelium_adj['cell_lineage'] = 'epithelium'
        self.myeloid_adj = pd.read_csv(os.path.join(dir_gg_adj, 'gg_adj_Myeloid.csv'), names=['gene','target','importance'])
        self.myeloid_adj['cell_lineage'] = 'myeloid'
        self.stroma_adj = pd.read_csv(os.path.join(dir_gg_adj, 'gg_adj_Stroma.csv'), names=['gene','target','importance'])
        self.stroma_adj['cell_lineage'] = 'stroma'
        self.tcell_adj = pd.read_csv(os.path.join(dir_gg_adj, 'gg_adj_T_Cell.csv'), names=['gene','target','importance'])
        self.tcell_adj['cell_lineage'] = 't_cell'
        self.gene_gene_adj = [self.bcell_adj, self.epithelium_adj, self.myeloid_adj, self.stroma_adj, self.tcell_adj]


    def get_reactome(self) -> pd.DataFrame:
        """
        This function retrieves the reactome data, filters it based on geneset size, and returns it as a DataFrame.

        :return: A DataFrame with the filtered reactome data.
        :rtype: pandas.DataFrame
        """
        reactome = self.gmt_to_decoupler("data2/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
        # Filtering genesets to match behaviour of fgsea
        geneset_size = reactome.groupby("geneset").size()
        gsea_genesets = geneset_size.index[(geneset_size > 15) & (geneset_size < 300)]
        reactome = reactome[reactome.geneset.isin(gsea_genesets)]
        return reactome

    def get_regulon_genesets(self) -> pd.DataFrame:
        """
        Returns a DataFrame with all genesets for each unique transcription factor (TF) in the given regulon data.

        :param reg_df: The DataFrame containing the regulon data.
        :type reg_df: pandas.DataFrame
        :return: A DataFrame with all genesets for each unique TF.
        :rtype: pandas.DataFrame
        """
        df = pd.DataFrame()
        for TF in self.reg_df['TF'].unique():
            df = pd.concat([df, self.get_regulon_genes(self.reg_df, TF)], axis=0)
        df = df.drop(columns='importance')
        df = df.rename(columns={'target': 'genesymbol', 'TF': 'geneset'})
        df = df.reset_index(drop=True)
        return df
    
    def get_genesets(self) -> pd.DataFrame:
        """
        This function concatenates the reactome and regulon dataframes, resets the index, and returns the result.

        :param reactome: The DataFrame containing the reactome data. Default is reactome.
        :type reactome: pandas.DataFrame, optional
        :param reg_df: The DataFrame containing the regulon data. Default is regulon_geneset.
        :type reg_df: pandas.DataFrame, optional
        :return: A DataFrame with the concatenated reactome and regulon data.
        :rtype: pandas.DataFrame
        """
        geneset_df = pd.concat([self.reactome, self.regulon_geneset], axis=0)
        geneset_df = geneset_df.reset_index(drop=True)
        return geneset_df



    # find all regulons that have GOI (CASP8) in their target genes
    def find_TFs(self, df: pd.DataFrame, GOI: str) -> np.array:
        """
        This function finds the transcription factors (TFs) that regulate a given gene of interest (GOI).

        :param df: The DataFrame containing the regulon data.
        :type df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        :return: An array of TFs that regulate the GOI.
        :rtype: numpy.array
        """
        goi_regulons = df[df['TargetGenes'].str.contains(str(GOI+'\''))]
        return goi_regulons['TF'].values


    def make_regulon_dataframe(self, df: pd.DataFrame, TF: str) -> pd.DataFrame:
        """
        This function creates a DataFrame for a given transcription factor (TF) with its target genes and their importance.

        :param df: The DataFrame containing the regulon data.
        :type df: pandas.DataFrame
        :param TF: The transcription factor.
        :type TF: str
        :return: A DataFrame with the target genes of the TF, their importance, the TF, and the group (TF_regulon).
        :rtype: pandas.DataFrame
        """
        reg_df = pd.DataFrame()
        for i in df[df['TF'] == TF]['TargetGenes']:
            temp = eval(i)
            reg_df = pd.concat([reg_df, pd.DataFrame(temp)], axis=0)
        reg_df = reg_df.reset_index(drop=True)
        reg_df = reg_df.drop_duplicates()
        reg_df['TF'] = TF
        reg_df = reg_df.rename(columns={0: 'target', 1: 'importance', 2: 'TF'})
        reg_df = reg_df.sort_values(by='importance', ascending=False)
        #reg_df = reg_df.head(500) # keep top x target genes per regulon here? supported by motif analysis
        reg_df['group'] = str(TF + "_regulon")
        return reg_df


    def make_adj_df(self, adj_df: pd.DataFrame, GOI: str) -> pd.DataFrame:
        """
        This function creates a DataFrame for a given gene of interest (GOI) with its adjacencies sorted by importance.

        :param adj_df: The DataFrame containing the adjacency data.
        :type adj_df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        :return: A DataFrame with the adjacencies of the GOI, sorted by importance, and a group column set to 'adjacencies'.
        :rtype: pandas.DataFrame
        """
        adj_interest = adj_df[adj_df['target'] == GOI]
        adj_interest = adj_interest.sort_values(by='importance', ascending=False)
        #adj_interest = adj_interest.head(15) # select top n 'important' TFs, threshold can be adjusted
        adj_interest['group'] = 'adjacencies'
        return adj_interest


    def make_goi_grn(self, GOI: str) -> pd.DataFrame:
        """
        This function creates a gene regulatory network (GRN) for a given gene of interest (GOI).

        :param GOI: The gene of interest.
        :type GOI: str
        :param df: The DataFrame containing the regulon data. Default is reg_df.
        :type df: pandas.DataFrame, optional
        :return: A DataFrame representing the GRN of the GOI, sorted by importance and only including target genes that appear more than once.
        :rtype: pandas.DataFrame
        """
        goi_regulons = self.find_TFs(self.reg_df, GOI)
        goi_grn = pd.DataFrame()
        for i in goi_regulons:
            goi_grn = pd.concat([goi_grn, self.make_regulon_dataframe(df=self.reg_df, TF=i)], axis=0)
        goi_grn.drop_duplicates(subset=['importance', 'TF', 'group'], keep="first", inplace=True)
        goi_grn = goi_grn.sort_values(by='importance', ascending=False)
        goi_grn = goi_grn[goi_grn.duplicated(subset=['target'], keep=False)]
        return goi_grn



    def get_entrez_gene_summary(self, gene_name: str, email: str, organism: str = "human", max_gene_ids: int = 10) -> dict:
        """
        Returns the 'Summary' contents for provided input gene from the Entrez Gene database.

        :param gene_name: Official (HGNC) gene name (e.g., 'KAT2A')
        :type gene_name: str
        :param email: Required email for making requests
        :type email: str
        :param organism: Filters results only to match organism. Set to None to return all organism unfiltered. Default is 'human'.
        :type organism: str, optional
        :param max_gene_ids: Sets the number of Gene ID results to return (absolute max allowed is 10K). Default is 100.
        :type max_gene_ids: int, optional
        :return: Summaries for all gene IDs associated with gene_name (where: keys → [orgn][gene name], values → gene summary)
        :rtype: dict
        """
        Entrez.email = email

        query = (
            f"{gene_name}[Gene Name]"
            if not organism
            else f"({gene_name}[Gene Name]) AND {organism}[Organism]"
        )
        handle = Entrez.esearch(db="gene", term=query, retmax=max_gene_ids)
        record = Entrez.read(handle)
        handle.close()

        gene_summaries = defaultdict(dict)
        gene_ids = record["IdList"]

        #print(f"{len(gene_ids)} gene IDs returned associated with gene {gene_name}.")

        for gene_id in gene_ids:
            #print(f"\tRetrieving summary for {gene_id}...")
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="docsum")
            gene_dict = xmltodict.parse(
                "".join([x.decode(encoding="utf-8") for x in handle.readlines()]),
                dict_constructor=dict,
            )
            gene_docsum = gene_dict["eSummaryResult"]["DocumentSummarySet"][
                "DocumentSummary"
            ]
            name = gene_docsum.get("Name")
            summary = gene_docsum.get("Summary")
            gene_organism = gene_docsum.get("Organism")["CommonName"]
            gene_summaries[gene_organism][name] = summary
            handle.close()
            time.sleep(0.34)  # Requests to NCBI are rate limited to 3 per second

        return gene_summaries


    # method to format gene summary from NCBI to print readably
    def format_gene_summary(self, df: pd.DataFrame, GOI: str) -> None:
        """
        This function formats and prints the gene summary from NCBI in a readable way.

        :param df: The DataFrame containing the gene data.
        :type df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        """
        gene_summary = self.get_gene_summary(df, GOI)
        i = 0
        for gene_name, summary in gene_summary.items():
            if i == 0:
                # skip the first gene summary
                i += 1
                continue
            print(f"\n{gene_name}: \n")
            text = summary
            # split gene summary text into separate strings after each 10th white space
            n = 12
            try:
                groups = text.split(' ')
                text = [' '.join(groups[i:i+n]) for i in range(0, len(groups), n)]
                for t in text:
                    print(f"\t{t}")
            except:
                print(f"\tNo summary available for {gene_name}.")


    def get_gene_summary(self, df: pd.DataFrame, GOI: str, email: str = 'samantha.bening@helmholtz-munich.de') -> dict:
        """
        This function gets the gene summary from NCBI for a given gene of interest (GOI) and its transcription factors (TFs).

        :param df: The DataFrame containing the gene data.
        :type df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        :param email: The email to use for making requests to NCBI. Default is 'samantha.bening@helmholtz-munich.de'.
        :type email: str, optional
        :return: A dictionary with the gene summaries, where the keys are the gene names and the values are the summaries.
        :rtype: dict
        """
        gene_summaries = {'Gene': 'Summary'}

        for gene in ([GOI] + df[df['target'] == GOI]['TF'].unique().tolist()):
            gene_summaries[gene] = self.get_entrez_gene_summary(gene, email)['human'][gene]
        return gene_summaries


    # summary of GOI and its regulons
    def GOI_network_stats(self, df: pd.DataFrame, GOI: str) -> None:
        """
        This function prints a summary of a given gene of interest (GOI) and its regulons.

        :param df: The DataFrame containing the gene data.
        :type df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        """
        df = df[df['target'] == GOI]
        print(f"Summary of {GOI}:\n")
        print(f"There are {len(df)} regulons that have {GOI} in their target genes.\n")
        print(f"Regulons that have {GOI} in their target genes:\n")
        print(f"\t(TF: GRNBoost2 Importance Score)")
        for row, col in df.iterrows():
            print(f"\t{col['TF']}: {round(col['importance'], 3)}")
        print("\n")
        df_adj = self.adj_df
        df_adj = self.make_adj_df(df_adj, GOI)
        # filter df_adj for any TFs that are NOT in df
        df_adj = df_adj[~df_adj['TF'].isin(df['TF'])]
        print(f"There are {len(df_adj)} TFs for {GOI} that were NOT supported by a regulon (motif analysis),")
        print(f"here are the top 10:\n")
        print(f"\t(TF: GRNBoost2 Importance Score)")
        counter = 0
        for row, col in df_adj.iterrows():
            if counter > 10:
                break
            print(f"\t{col['TF']}: {round(col['importance'], 3)}")
            counter += 1


    # cell-type specific analysis of TFs and regulons
    def plot_regulon_expression(self, df: pd.DataFrame, GOI: str) -> None:
        """
        This function plots the expression of a given gene of interest (GOI) and its regulons.

        :param df: The DataFrame containing the gene data.
        :type df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        :param adata: The AnnData object containing the single-cell data. Default is adata.
        :type adata: anndata.AnnData, optional
        """
        genesets = self.adata.obsm['aucell_estimate'].columns
        # find any genesets that are regulons for GOI
        regulons_in_genesets = []
        TFs = df[df['target'] == GOI]
        TFs
        for TF in TFs['TF']:
            for geneset in genesets:
                if (str(TF) + "_REGULON") in geneset:
                    regulons_in_genesets.append(geneset)

        self.adata.obs[regulons_in_genesets] = self.adata.obsm["aucell_estimate"][regulons_in_genesets]

        print(f"The following regulons are present in relevant AUCell cellular-level analysis: {regulons_in_genesets}")
        print(f"Plotting cell-type specific expression or AUCell score for GOI and regulons")
        sc.pl.umap(
            self.adata,
            color=["celltype_l2", GOI] + regulons_in_genesets,
            frameon=False,
            ncols=2,
            wspace=0.4,
        )

    def genegene_importance_histograms(self, log_scale=False, xlim=10):
        fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows = 1, ncols = 5, figsize=(14,7), sharey=True)
        fig.suptitle("Gene-Gene Adj. Importance Score Distributions")
        fig.supxlabel("Importance score")
        if log_scale:
            fig.supylabel("log(Number of genes)")
        else:
            fig.supylabel("Number of genes")
        
        ax1.hist(self.bcell_adj['importance'], bins=15, log=log_scale, rwidth = 0.9, range=(0, xlim))
        ax1.set_title("B Cells")
        ax2.hist(self.epithelium_adj['importance'], bins=15, log=log_scale, rwidth = 0.9, range=(0, xlim))
        ax2.set_title("Epithelium Cells")
        ax3.hist(self.myeloid_adj['importance'], bins=15, log=log_scale, rwidth = 0.9, range=(0, xlim))
        ax3.set_title("Myeloid Cells")
        ax4.hist(self.stroma_adj['importance'], bins=15, log=log_scale, rwidth = 0.9, range=(0, xlim))
        ax4.set_title("Stroma Cells")
        ax5.hist(self.tcell_adj['importance'], bins=15, log=log_scale, rwidth = 0.9, range=(0, xlim))
        ax5.set_title("T Cells")

        fig.tight_layout()
        plt.show()
            



    def make_network(self, df: pd.DataFrame, GOI: str, direct_TF: bool = True, top_n: int = None, out_file: str = 'src/SCENICfiles/network') -> str:
        """
        This function creates a network visualization of a given gene of interest (GOI) and its regulons.

        :param df: The DataFrame containing the gene data.
        :type df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        :param direct_TF: If True, only direct neighbors of the GOI are included. If False, all neighbors are included. Default is True.
        :type direct_TF: bool, optional
        :param top_n: The number of top neighbors of each regulon (other than GOI) to include. If None, no regulon neighbors (other than GOI) are included. Default is None.
        :type top_n: int, optional
        :param out_file: The path to the output file where the network visualization will be saved. Default is 'src/gene_report/goi_network.html'.
        :type out_file: str, optional
        """
        net = Network(notebook=True, height='600px', width='700px', bgcolor="#FFFFFF", 
                    font_color="black", cdn_resources='in_line')

        # filter for only direct neighbors of GOI
        if direct_TF and top_n == None:
            df = df[df['target'] == GOI]

        if direct_TF and top_n != None:
            temp = df.groupby('TF').head(top_n)
            df = pd.concat([temp, df[df['target'] == GOI]], axis=0)
        
        if not direct_TF and top_n != None:
            df = df.groupby('TF').head(top_n)

        # assign colors to regulons
        groups = {}
        colors = sns.color_palette("Set2", len(df['group'].unique())).as_hex()
        i = 0
        for group in df['group'].unique():
            groups[group] = colors[i]
            i = i + 1

        # build network
        sources = df['TF']
        targets = df['target']
        weights = df['importance']
        group = df['group']

        edge_data = zip(sources, targets, weights, group)

        for e in edge_data:
                        src = e[0]
                        dst = e[1]
                        w = e[2]

                        net.add_node(src, src, title=str(e[3]), color=groups[e[3]])
                        net.add_node(dst, dst, title=str(e[3]), color=groups[e[3]])
                        net.add_edge(src, dst, value=w, label=w, title=w, color=groups[e[3]])

        neighbor_map = net.get_adj_list()

        # add neighbor data to node hover data
        for node in net.nodes:
            node["value"] = len(neighbor_map[node["id"]])
            if node['id'] == GOI:
                    node['color'] = '#001861'


        net.toggle_physics(False)
        network_file = f"{out_file}_regulons_{GOI}.html"
        net.save_graph(network_file)
        return network_file


        
    def make_gene_gene_network(self, GOI: str, top_n: int = None, out_file: str = 'src/SCENICfiles/network') -> str:
        """
        This function creates a network visualization of a given gene of interest (GOI) and its co-expressed genes.

        :param df: The DataFrame containing the gene data.
        :type df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        :param top_n: The number of top neighbors of each regulon (other than GOI) to include. If None, no regulon neighbors (other than GOI) are included. Default is None.
        :type top_n: int, optional
        :param out_file: The path to the output file where the network visualization will be saved. 
        :type out_file: str, optional
        """
        net = Network(notebook=True, height='600px', width='700px', bgcolor="#FFFFFF", 
                    font_color="black", cdn_resources='in_line')

        df = pd.DataFrame()
        for df_i in self.gene_gene_adj:
             df_i = df_i[(df_i['target'] == GOI) | (df_i['gene'] == GOI)]
             df = pd.concat([df, df_i.head(top_n)], axis=0)

        # assign colors to regulons
        groups = {}
        colors = sns.color_palette("hls", len(df['cell_lineage'].unique())).as_hex()
        i = 0
        for group in df['cell_lineage'].unique():
            groups[group] = colors[i]
            i = i + 1

        # build network
        sources = df['gene']
        targets = df['target']
        weights = df['importance']
        group = df['cell_lineage']

        edge_data = zip(sources, targets, weights, group)

        for e in edge_data:
                        src = e[0]
                        dst = e[1]
                        w = e[2]

                        net.add_node(src, title=str(e[3]), color=groups[e[3]])
                        net.add_node(dst, title=str(e[3]), color=groups[e[3]])
                        net.add_edge(src, dst, title=w, color=groups[e[3]], label=w, value=w)

        neighbor_map = net.get_adj_list()

        # add neighbor data to node hover data
        for node in net.nodes:
            node["value"] = len(neighbor_map[node["id"]])
            if node['id'] == GOI:
                    node['color'] = '#001861'

        net.toggle_physics(False)
        network_file1 = f"{out_file}_genegene_{GOI}.html"
        net.save_graph(network_file1)
        return network_file1



    def show_network(self, GOI: str, type: str ='gene_gene', top_n: int = 5) -> None:
        """
        This function displays the HTML content of a given file.
        """
        if type == 'regulon':
             df = self.make_goi_grn(GOI)
             network_file = self.make_network(GOI=GOI, df=df, direct_TF=True, top_n=top_n)
        elif type =='gene_gene':
             network_file = self.make_gene_gene_network(GOI=GOI, top_n=top_n)
                # Read the contents of the HTML file
        
        with open(network_file, 'r') as file:
            html_content = file.read()
        # Display the HTML content
        display(HTML(html_content))



    ### GSEA stuff here ###
        
    @staticmethod
    def gmt_to_decoupler(pth: Path) -> pd.DataFrame:
        """
        Parse a gmt file to a decoupler pathway dataframe.

        :param pth: The path to the gmt file.
        :type pth: Path
        :return: A DataFrame with the geneset and genesymbol from the gmt file.
        :rtype: pandas.DataFrame
        """
        pathways = {}
        with Path(pth).open("r") as f:
            for line in f:
                name, _, *genes = line.strip().split("\t")
                pathways[name] = genes

        return pd.DataFrame.from_records(
            chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
            columns=["geneset", "genesymbol"],
        )



    def get_regulon_genes(self, reg_df: pd.DataFrame, TF: str) -> pd.DataFrame:
        """
        Returns a DataFrame with all target genes and importance scores for a given transcription factor (TF).

        :param reg_df: The DataFrame containing the regulon data.
        :type reg_df: pandas.DataFrame
        :param TF: The name of the transcription factor.
        :type TF: str
        :return: A DataFrame with all target genes and importance scores for the given TF.
        :rtype: pandas.DataFrame
        """
        all_targets = reg_df[reg_df['TF'] == TF]['TargetGenes']
        df = pd.DataFrame()
        for i in all_targets:
            i = eval(i)
            df = pd.concat([df, pd.DataFrame(i)], axis=0)
            
        df = df.drop_duplicates()
        df = df.rename(columns={0: 'target', 1: 'importance'})
        df = df.sort_values(by='importance', ascending=False)
        df['TF'] = TF + "_REGULON"
        return df # returns a dataframe with all target genes and importance scores for a given TF



    def get_goi_pathways(self, GOI, method='spearman'):
        """
        This function calculates and ranks the correlation between each pathway's AUCell score and the expression of the gene of interest (GOI).

        :param GOI: The gene of interest.
        :type GOI: str
        :param geneset_df: The DataFrame containing the geneset data. Default is geneset_df.
        :type geneset_df: pandas.DataFrame, optional
        :param adata: The AnnData object containing the single-cell data. Default is adata.
        :type adata: anndata.AnnData, optional
        :param method: The method to use for calculating correlation. Default is 'spearman'.
        :type method: str, optional
        :return: A DataFrame with the pathways for the GOI, sorted by absolute correlation value.
        :rtype: pandas.DataFrame
        """
        # correlation between each pathways AUCell score and GOI expression (ONLY GENESETS AFTER FILTERING)
        pathways_goi = self.geneset_df[(self.geneset_df['genesymbol'] == GOI)]

        # get the AUCell scores for the genes in the pathway
        goi_expression = self.adata.to_df()[GOI]

        def calc_correlation(column): # Spearman
            return column.corr(goi_expression, method='spearman')

        pathways_goi.loc[:, 'correlation'] = self.adata.obsm['aucell_estimate'][pathways_goi['geneset']].apply(calc_correlation).values
        pathways_goi = pathways_goi.sort_values('correlation', key=pd.Series.abs, ascending=False)
        pathways_goi = pathways_goi.reset_index(drop=True)
        return pathways_goi

    # find gene set of specific regulon (for now automatically set to the first regulon in the list)
    def get_regulon_geneset(self, regulon=None) -> List[str]:
        """
        This function returns the gene set of a specific regulon. If no regulon is specified, it defaults to the first regulon in the list.

        :param df: The DataFrame containing the geneset data. Default is geneset_df.
        :type df: pandas.DataFrame, optional
        :param regulon: The name of the regulon. If None, the first regulon in the list is used. Default is None.
        :type regulon: str, optional
        :return: A list of genes in the specified regulon.
        :rtype: List[str]
        """
        if regulon == None:
            regulon = self.geneset_df[self.geneset_df['geneset'].str.match(r'(.*)_REGULON')].iloc[0]['geneset']
        
        return list(self.geneset_df[self.geneset_df['geneset'] == regulon]['genesymbol'].values)


    def plot_pathways(self, df: pd.DataFrame, GOI: str) -> None:
        """
        This function plots the UMAP of the top 5 pathways along with the cell type and the gene of interest (GOI).

        :param df: The DataFrame containing the pathway data.
        :type df: pandas.DataFrame
        :param GOI: The gene of interest.
        :type GOI: str
        :param adata: The AnnData object containing the single-cell data. Default is adata.
        :type adata: anndata.AnnData, optional
        """
        # TODO: sort types of pathways (e.g. Reactome vs regulon) in columns under cell type and GOI expression reference
        top_pathways = list(df['geneset'].values[:4])
        self.adata.obs[top_pathways] = self.adata.obsm["aucell_estimate"][top_pathways]

        sc.pl.umap(
            self.adata,
            color=["celltype_l2", GOI] + top_pathways,
            frameon=False,
            ncols=2,
            wspace=0.4,
        )


    def gGOSt(self, regulon: str) -> pd.DataFrame:
        """
        This function performs g:Profiler g:GOSt analysis for pathways in a given regulon, plots the results, and returns the result DataFrame.
        Currently only searching in the REAC, KEGG, and GO:BP databases.
        Some regulons may not have significant pathways, in which case the function will return a message.

        :param regulon: The name of the regulon.
        :type regulon: str
        :return: A DataFrame with the g:Profiler g:GOSt analysis results for the given regulon.
        :rtype: pandas.DataFrame
        """
        reg_geneset = self.get_regulon_geneset(regulon=regulon)
        r = requests.post(
            url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
            json={
                'organism':'hsapiens',
                'query':reg_geneset,
                'sources': ['REAC', 'KEGG', 'GO:BP']
            }
            )
        result = r.json()['result']
        result = pd.DataFrame(result)
        try:
            # manipulate result df for plotting
            result['nlog(p)'] = -np.log10(result['p_value'])
            result = result.sort_values('source', ascending=False)
            result = result.reset_index(drop=True)
            result = result[result['term_size'] < 500] # filter out large pathways in databases
            # visualize p-values in dot plot nicely
            fig = px.scatter(result, x=result.index, y="nlog(p)", color="source", hover_name='name', hover_data={'source': False, 'native': True})
            fig.update_layout(title=f"g:Profiler g:GOSt analysis for pathways in {regulon}", xaxis_title="Pathway", yaxis_title="-log10(p-value)")
            fig.update_xaxes(ticktext=result[result['source'].duplicated() == False]['source'], tickvals=result[result['source'].duplicated() == False].index)
            fig.update_traces(marker=dict(size=11, opacity=0.8), showlegend=False )
            fig.show()
            return result
        except:
            print('No significant pathways found for this regulon. Try another regulon or gene set.')


    def gGOSt_listed(self, GOI: str) -> str:
        """
        Print list of top three significantly differentially expressed in each 
        """
        print(f"Top 3 significantly differentially expressed pathways in genes co-expressed with {GOI}, per cell lineage.\nThis is in REACTOME and KEGG databases.\n")
        for lineage in self.gene_gene_adj:
            # get list of gene names in each gene_gene adjacencies file
            lister = lineage[(lineage['gene'] == GOI) | (lineage['target'] == GOI)]
            gene_list = lister['gene'].unique()
            target_list = lister['target'].unique()
            new_list = np.append(gene_list, target_list)
            new_list = np.unique(new_list).tolist()
            r = requests.post(
                url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
                json={
                    'organism':'hsapiens',
                    'query':new_list,
                    'sources': ['REAC', 'KEGG']
                }
                )
            result = r.json()['result']
            result = pd.DataFrame(result)
            # sort by p-value
            result = result.sort_values('p_value', ascending=True)
            result = result[result['term_size'] < 500]
            result = result.reset_index(drop=True)
            result = result.head(3)
            print(f"In {lineage['cell_lineage'].unique().tolist()[0]} cells:")
            for index, row in result.iterrows():
                print(f"{row['name']} (p-value: {row['p_value']})")
            if len(result) == 0:
                print("No significant pathways found.")
            print("\n")



