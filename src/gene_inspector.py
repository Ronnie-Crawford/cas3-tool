import numpy as np
import pandas as pd
from pyfaidx import Fasta

from config_loader import config

def find_size(genes_df : pd.DataFrame):

    """
    Finds the size of each gene, irrespective of gene direction.

    Arguments:
        - genes_df (pd.DataFrame): DataFrame containing metrics on genes.

    Returns:
        - genes_df (pd.DataFrame): DataFrame containing metrics on genes.
    """

    genes_df["size"] = genes_df.apply(lambda gene : gene["end"] - gene["start"] if gene["end"] > gene["start"] else gene["start"] - gene["end"], axis = 1)

    return genes_df

def find_isolation(genes_df):

    """
    Find the name and position of genes nearest upstream and downstream to each gene,
    then finds the distance to the closest neighbour of each gene.

    Arguments:
        - genes_df (pd.DataFrame): DataFrame containing metrics on genes.

    Returns:
        - genes_df (pd.DataFrame): DataFrame containing metrics on genes.
    """

    chromosome_dfs = []

    for chromosome in config["CHROMOSOMES"]:

        print(f"Finding attributes of genes on chromosome {chromosome}.")

        chromosome_df = genes_df[genes_df["seqname"] == chromosome]

        if chromosome_df.empty:

            print(f"No essential genes found on chromosome {chromosome}.")
            continue

        chromosome_df = chromosome_df.sort_values(by = 'start').reset_index(drop = True)
        chromosome_df[[
            'nearest_upstream_name', 'nearest_upstream_start', 'nearest_upstream_end',
            'nearest_downstream_name', 'nearest_downstream_start', 'nearest_downstream_end'
            ]] = chromosome_df.apply(lambda gene : find_nearest_genes(gene, chromosome_df), axis = 1)

        chromosome_df["distance_to_upstream"] = chromosome_df["start"] - chromosome_df["nearest_upstream_end"]
        chromosome_df["distance_to_downstream"] = chromosome_df["nearest_downstream_start"] - chromosome_df["end"]
        chromosome_df["distance_to_nearest_gene"] = chromosome_df[['distance_to_upstream', 'distance_to_downstream']].min(axis=1)

        chromosome_df = chromosome_df[chromosome_df["distance_to_upstream"].notna()]
        chromosome_df = chromosome_df[chromosome_df["distance_to_downstream"].notna()]

        chromosome_dfs.append(chromosome_df)

    genes_df = pd.concat(chromosome_dfs, ignore_index=True)

    return genes_df

def find_nearest_genes(gene, chromosome_df):

    """
    The function used by "find_isolation()" to find neighbouring genes for a given gene

    Arguments:
        - gene (pd.Series): A row from the chromosome_df DataFrame.
        - chromosome_df (pd.DataFrame): A DataFrame of all available genes for a given chromosome.

    Returns:
        - Neighbouring gene values (pd.Series): The name, start and end of upstream and downstream neighbouring genes.
    """

    nearest_upstream_name = None
    nearest_upstream_start = None
    nearest_upstream_end = None

    nearest_downstream_name = None
    nearest_downstream_start = None
    nearest_downstream_end = None

    upstream_candidates_df = chromosome_df[chromosome_df['end'] < gene['start']]
    upstream_candidates_df = upstream_candidates_df[upstream_candidates_df["strand"] == gene["strand"]]
    upstream_candidates_df =upstream_candidates_df.sort_values(by='end').reset_index(drop=True)

    downstream_candidates_df = chromosome_df[chromosome_df['start'] > gene['end']]
    downstream_candidates_df = downstream_candidates_df[downstream_candidates_df["strand"] == gene["strand"]]
    downstream_candidates_df = downstream_candidates_df.sort_values(by='start').reset_index(drop=True)

    if not upstream_candidates_df.empty:

        nearest_upstream_name = upstream_candidates_df.iloc[-1]['gene']
        nearest_upstream_start = upstream_candidates_df.iloc[-1]['start']
        nearest_upstream_end = upstream_candidates_df.iloc[-1]['end']

    if not downstream_candidates_df.empty:

        nearest_downstream_name = downstream_candidates_df.iloc[0]['gene']
        nearest_downstream_start = downstream_candidates_df.iloc[0]['start']
        nearest_downstream_end = downstream_candidates_df.iloc[0]['end']

    return pd.Series([nearest_upstream_name, nearest_upstream_start, nearest_upstream_end, nearest_downstream_name, nearest_downstream_start, nearest_downstream_end])

def find_search_windows(genes_df):

    genome = Fasta(config["GENOME_SEQUENCE_PATH"])
    chromosome_dfs = []

    for chromosome in config["CHROMOSOMES"]:

        print(f"Finding FASTA sequence around genes in chromosome {chromosome}.")

        chromosome_df = genes_df[genes_df["seqname"] == chromosome]

        if chromosome_df.empty:

            print(f"No essential genes found in chromosome {chromosome}.")
            continue

        chromosome_df = chromosome_df.sort_values(by = 'start').reset_index(drop = True)

        chromosome_df["upstream_window_start"] = np.where(
            (chromosome_df["start"] - chromosome_df["nearest_upstream_end"]) < config["SEARCH_WINDOW_UPSTREAM_SIZE_MAX"],
            chromosome_df["nearest_upstream_end"],
            chromosome_df["start"] - config["SEARCH_WINDOW_UPSTREAM_SIZE_MAX"]
            )
        chromosome_df["upstream_window_end"] = chromosome_df["start"] - 1
        chromosome_df["downstream_window_start"] = chromosome_df["end"] + 1
        chromosome_df["downstream_window_end"] = np.where(
            (chromosome_df["nearest_downstream_start"] - chromosome_df["end"]) < config["SEARCH_WINDOW_DOWNSTREAM_SIZE_MAX"],
            chromosome_df["nearest_downstream_start"],
            chromosome_df["end"] + config["SEARCH_WINDOW_DOWNSTREAM_SIZE_MAX"]
            )
        chromosome_df["number_motifs_upstream"] = chromosome_df.apply(
            lambda gene : fetch_fasta_sequence(genome, chromosome, gene["upstream_window_start"], gene["upstream_window_end"]), axis = 1
        )
        chromosome_df["number_motifs_downstream"] = chromosome_df.apply(
            lambda gene : fetch_fasta_sequence(genome, chromosome, gene["downstream_window_start"], gene["downstream_window_end"]), axis = 1
        )

        chromosome_dfs.append(chromosome_df)

    genes_df = pd.concat(chromosome_dfs, ignore_index=True)

    return genes_df

def fetch_fasta_sequence(genome, chromosome, start, end):

    chromosome = f"{chromosome}"
    start = int(start)
    end = int(end)

    searched_sequence = genome[chromosome][start:end].seq
    number_of_motifs = 0

    for motif in config["MOTIFS"]:

        number_of_motifs += searched_sequence.count(motif)

    return int(number_of_motifs)
