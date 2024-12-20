from pathlib import Path
import os
import decoupler
import pandas as pd 
import anndata as ad
pd.options.mode.chained_assignment = None  # default='warn'
from pathlib import Path
from itertools import chain, repeat
from argparse import ArgumentParser

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


def get_regulon_genes(reg_df: pd.DataFrame, TF: str) -> pd.DataFrame:
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
    return df


def get_reactome(gmt_file: str) -> pd.DataFrame:
    """
    Retrieves the reactome data, filters it based on geneset size, and returns it as a DataFrame.

    :param gmt_file: The path to the gmt file.
    :return: A DataFrame with the filtered reactome data.
    """
    reactome = gmt_to_decoupler(Path(gmt_file))
    geneset_size = reactome.groupby("geneset").size()
    gsea_genesets = geneset_size.index[(geneset_size > 15) & (geneset_size < 300)]
    reactome = reactome[reactome.geneset.isin(gsea_genesets)]
    return reactome


def get_regulon_genesets(reg_df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a DataFrame with all genesets for each unique transcription factor (TF) in the given regulon data.

    :param reg_df: The DataFrame containing the regulon data.
    :return: A DataFrame with all genesets for each unique TF.
    """
    df = pd.DataFrame()
    for TF in reg_df['TF'].unique():
        df = pd.concat([df, get_regulon_genes(reg_df, TF)], axis=0)
    df = df.drop(columns='importance')
    df = df.rename(columns={'target': 'genesymbol', 'TF': 'geneset'})
    df = df.reset_index(drop=True)
    return df


def get_genesets(reactome: pd.DataFrame, reg_geneset: pd.DataFrame) -> pd.DataFrame:
    """
    Concatenates the reactome and regulon dataframes, resets the index, and returns the result.

    :param reactome: The DataFrame containing the reactome data.
    :param reg_df: The DataFrame containing the regulon data.
    :return: A DataFrame with the concatenated reactome and regulon data.
    """
    geneset_df = pd.concat([reactome, reg_geneset], axis=0)
    geneset_df = geneset_df.reset_index(drop=True)
    return geneset_df

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("--data", type=str, default = 'veo_ibd_balanced.h5ad')
    parser.add_argument("--output", type=str, default = 'src/SCENICfiles/new')
    args = parser.parse_args()
    data_file = f'data2/{args.data}'
    output_dir = args.output
    print("\tArgs read in")

    wdir = "/lustre/groups/ml01/workspace/christopher.lance/genereporter/"
    os.chdir( wdir )

    # read regulon data
    regulon = pd.read_csv(f'{output_dir}/regulons_output.csv', skiprows=3, header=0,
                           names=['TF','MotifID','AUC','NES','MotifSimilarityQvalue','OrthologousIdentity','Annotation',
                           'Context','TargetGenes','RankAtMax'])
    def clean_target_genes(row: pd.Series) -> list:
        return eval(row['TargetGenes'])
    regulon.apply(clean_target_genes, axis=1)
    print('\tCleaned regulon')

    # usage:
    reactome_df = get_reactome("data2/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
    regulon_genesets_df = get_regulon_genesets(regulon)  
    geneset_df = get_genesets(reactome_df, regulon_genesets_df)

    print('\tRead in genesets (regulons and reactome)')

    adata = ad.read_h5ad(data_file)


    decoupler.run_aucell(
        adata,
        geneset_df,
        source="geneset",
        target="genesymbol",
        use_raw=False,
    )

    print("\tAUCell decoupler done")


    print(adata.obs.columns)

    adata.write_h5ad(f"{output_dir}/veo_ibd_balanced_aucell_final.h5ad")

    print("\tAUCell written to adata file done")