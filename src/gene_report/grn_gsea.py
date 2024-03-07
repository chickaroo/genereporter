import warnings
warnings.filterwarnings('ignore')
import os
import pandas as pd
from pyvis.network import Network
import seaborn as sns
import scanpy as sc
import numpy as np
import time
import xmltodict
from collections import defaultdict
from Bio import Entrez
from pathlib import Path
import requests
import plotly.express as px
from IPython.display import display, HTML

# set a working directory
wdir = "/Users/samibening/Projects/Bachelor/"
os.chdir( wdir )

adata = sc.read_h5ad('data/output/adata_aucell.h5ad')

def clean_target_genes(row):
    return eval(row['TargetGenes'])

def read_grn(f_adj='SCENICfiles/adj.csv', f_reg='SCENICfiles/reg.csv'):
    adjacencies = pd.read_csv(f_adj)
    regulon = pd.read_csv(f_reg)
    regulon.apply(clean_target_genes, axis=1)
    return adjacencies, regulon

# find all regulons that have GOI (CASP8) in their target genes
def find_TFs(df, GOI):
    goi_regulons = df[df['TargetGenes'].str.contains(str(GOI+'\''))]
    return goi_regulons['TF'].values

def make_regulon_dataframe(df, TF):
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

def make_adj_df(adj_df, GOI): # question is if we even want the adjacencies in the netowork? 
    adj_interest = adj_df[adj_df['target'] == GOI]
    adj_interest = adj_interest.sort_values(by='importance', ascending=False)
    #adj_interest = adj_interest.head(15) # select top n 'important' TFs, threshold can be adjusted
    adj_interest['group'] = 'adjacencies'
    return adj_interest

def make_goi_grn(df, GOI):
    goi_regulons = find_TFs(df, GOI)
    goi_grn = pd.DataFrame()
    for i in goi_regulons:
        goi_grn = pd.concat([goi_grn, make_regulon_dataframe(df=df, TF=i)], axis=0)
    goi_grn.drop_duplicates(subset=['importance', 'TF', 'group'], keep="first", inplace=True)
    goi_grn = goi_grn.sort_values(by='importance', ascending=False)
    goi_grn = goi_grn[goi_grn.duplicated(subset=['target'], keep=False)]
    return goi_grn


def get_entrez_gene_summary(
    gene_name, email, organism="human", max_gene_ids=100
):
    """Returns the 'Summary' contents for provided input
    gene from the Entrez Gene database. 
    
    Args:
        gene_name (string): Official (HGNC) gene name 
           (e.g., 'KAT2A')
        email (string): Required email for making requests
        organism (string, optional): defaults to human. 
           Filters results only to match organism. Set to None
           to return all organism unfiltered.
        max_gene_ids (int, optional): Sets the number of Gene
           ID results to return (absolute max allowed is 10K).
        
    Returns:
        dict: Summaries for all gene IDs associated with 
           gene_name (where: keys → [orgn][gene name],
                      values → gene summary)
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
def format_gene_summary(df, GOI):
    gene_summary = get_gene_summary(df, GOI)
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

def get_gene_summary(df, GOI, email='samantha.bening@helmholtz-munich.de'):
    gene_summaries = {'Gene': 'Summary'}

    for gene in ([GOI] + df[df['target'] == GOI]['TF'].unique().tolist()):
        gene_summaries[gene] = get_entrez_gene_summary(gene, email)['human'][gene]
    return gene_summaries

# summary of GOI and its regulons
def GOI_network_stats(df, GOI):
    df = df[df['target'] == GOI]
    print(f"Summary of {GOI}:\n")
    print(f"There are {len(df)} regulons that have {GOI} in their target genes.\n")
    print(f"Regulons that have {GOI} in their target genes:\n")
    print(f"\t(TF: GENIE3 Importance Score)")
    for row, col in df.iterrows():
        print(f"\t{col['TF']}: {round(col['importance'], 3)}")
    print("\n")
    df_adj, temp = read_grn()
    df_adj = make_adj_df(df_adj, GOI)
    # filter df_adj for any TFs that are NOT in df
    df_adj = df_adj[~df_adj['TF'].isin(df['TF'])]
    print(f"There are {len(df_adj)} TFs for {GOI} that were NOT supported by a regulon (motif analysis),")
    print(f"here are the top 10:\n")
    print(f"\t(TF: GENIE3 Importance Score)")
    counter = 0
    for row, col in df_adj.iterrows():
        if counter > 10:
            break
        print(f"\t{col['TF']}: {round(col['importance'], 3)}")
        counter += 1


# cell-type specific analysis of TFs and regulons
def plot_regulon_expression(df, GOI, adata=adata):
    genesets = adata.obsm['aucell_estimate'].columns
    # find any genesets that are regulons for GOI
    regulons_in_genesets = []
    TFs = df[df['target'] == GOI]
    TFs
    for TF in TFs['TF']:
        for geneset in genesets:
            if (str(TF) + "_REGULON") in geneset:
                regulons_in_genesets.append(geneset)

    adata.obs[regulons_in_genesets] = adata.obsm["aucell_estimate"][regulons_in_genesets]

    print(f"The following regulons are present in relevant AUCell cellular-level analysis: {regulons_in_genesets}")
    print(f"Plotting cell-type specific expression or AUCell score for GOI and regulons")
    sc.pl.umap(
        adata,
        color=["celltypist_cell_label_coarse", GOI] + regulons_in_genesets,
        frameon=False,
        ncols=2,
        wspace=0.4,
    )

def make_network(df, GOI, out_file='src/gene_report/goi_network.html'):
    net = Network(notebook=True, height='500px', width='400px', bgcolor="#222222", font_color="white", cdn_resources='remote')

    # filter for only direct neighbors of GOI
    df = df[df['target'] == GOI]

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

                    net.add_node(src, src, title=src, color=groups[e[3]])
                    net.add_node(dst, dst, title=dst, color=groups[e[3]])
                    net.add_edge(src, dst, value=w, color=groups[e[3]])

    neighbor_map = net.get_adj_list()

    # add neighbor data to node hover data
    for node in net.nodes:
                    node["title"] = node['id']
                    node["value"] = len(neighbor_map[node["id"]])
                    node['label'] = node['id']

    net.toggle_physics(False)
    net.save_graph(out_file)


def show_network(out_file='src/gene_report/goi_network.html'):
    # Read the contents of the HTML file
    with open(out_file, 'r') as file:
        html_content = file.read()

    # Display the HTML content
    display(HTML(html_content))




### GSEA stuff here ###
    


def gmt_to_decoupler(pth: Path) -> pd.DataFrame:
    """
    Parse a gmt file to a decoupler pathway dataframe.
    """
    from itertools import chain, repeat

    pathways = {}

    with Path(pth).open("r") as f:
        for line in f:
            name, _, *genes = line.strip().split("\t")
            pathways[name] = genes

    return pd.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=["geneset", "genesymbol"],
    )

def get_reactome():
    reactome = gmt_to_decoupler("data/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
    # Filtering genesets to match behaviour of fgsea
    geneset_size = reactome.groupby("geneset").size()
    gsea_genesets = geneset_size.index[(geneset_size > 15) & (geneset_size < 300)]
    reactome = reactome[reactome.geneset.isin(gsea_genesets)]
    return reactome


def get_regulon_genes(reg_df, TF): # give TF name as string, e.g. 'KLF5'

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

def get_regulon_genesets(reg_df):
    df = pd.DataFrame()
    for TF in reg_df['TF'].unique():
        df = pd.concat([df, get_regulon_genes(reg_df, TF)], axis=0)
    df = df.drop(columns='importance')
    df = df.rename(columns={'target': 'genesymbol', 'TF': 'geneset'})
    df = df.reset_index(drop=True)
    return df


def get_genesets(reactome, reg_df):
    geneset_df = pd.concat([reactome, reg_df], axis=0)
    geneset_df = geneset_df.reset_index(drop=True)
    return geneset_df

def get_goi_pathways(geneset_df, GOI, adata=adata, method='spearman'):
    # correlation between each pathways AUCell score and CASP8 expression (ONLY GENESETS AFTER FILTERING!!)
    pathways_goi = geneset_df[(geneset_df['genesymbol'] == GOI)]

    # get the AUCell scores for the genes in the pathway
    goi_expression = adata.to_df()[GOI]

    def calc_correlation(column): # Spearman
        return column.corr(goi_expression, method='spearman')

    pathways_goi.loc[:, 'correlation'] = adata.obsm['aucell_estimate'][pathways_goi['geneset']].apply(calc_correlation).values
    pathways_goi = pathways_goi.sort_values('correlation', key=pd.Series.abs, ascending=False)
    pathways_goi = pathways_goi.reset_index(drop=True)
    return pathways_goi

# find gene set of specific regulon (for now automatically set to the first regulon in the list)
def get_regulon_geneset(df, regulon=None):
    if regulon == None:
        regulon = df[df['geneset'].str.match(r'(.*)_REGULON')].iloc[0]['geneset']
    
    return list(df[df['geneset'] == regulon]['genesymbol'].values)


def plot_pathways(df, GOI, adata=adata):
    top_pathways = list(df['geneset'].values[:5])
    adata.obs[top_pathways] = adata.obsm["aucell_estimate"][top_pathways]

    sc.pl.umap(
        adata,
        color=["celltypist_cell_label_coarse", GOI] + top_pathways,
        frameon=False,
        ncols=2,
        wspace=0.4,
    )


def gGOSt(reg_geneset, GOI):
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
        fig.update_layout(title=f"g:Profiler g:GOSt analysis for {GOI} pathways", xaxis_title="Pathway", yaxis_title="-log10(p-value)")
        fig.update_xaxes(ticktext=result[result['source'].duplicated() == False]['source'], tickvals=result[result['source'].duplicated() == False].index)
        fig.update_traces(marker=dict(size=11, opacity=0.8), showlegend=False )
        fig.show()
        return result
    except:
        print('No significant pathways found for this regulon. Try another regulon or gene set.')