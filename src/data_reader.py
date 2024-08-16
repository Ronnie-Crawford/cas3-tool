# Third party modules
import numpy as np
import pandas as pd
from pyfaidx import Fasta

# Local modules
from config_loader import config

def read_genome_sequence():

    try:

        print("Reading genome sequence file.")

        genome_sequence = Fasta(config["GENOME_SEQUENCE_PATH"])

    except:

        raise

    return genome_sequence

def read_essentiality():

    try:

        print("Reading essentiality file.")

        essentiality_df = pd.read_csv(config["ESSENTIALITY_PATH"], sep = "\t")
        essentiality_df = essentiality_df[essentiality_df["lfc"] < config["GENE_ESSENTIALITY_THRESHOLD"]]
        essentiality_df["gene"] = essentiality_df["gene"].str.upper()

    except:

        raise

    return essentiality_df

def read_genome_annotations():

    try:

        print("Reading annotations file.")

        annotation_column_names = [
            'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'
        ]
        annotation_column_dtypes = dtype = {
            'seqname': str, 'source': str, 'feature': str, 'start': int, 'end': int, 'score': str, 'strand': str, 'frame': str, 'attribute': str
        }
        annotations_df = pd.read_csv(
            config["GENOME_ANNOTATIONS_PATH"], sep = '\t', comment = '#', header = None, names = annotation_column_names, dtype = annotation_column_dtypes
        )
        annotations_df = annotations_df[annotations_df["feature"] == "gene"]
        annotations_df = annotations_df[annotations_df['seqname'].isin(config["CHROMOSOMES"])]
        annotations_df[[
            "gene_id", "gene_version", "gene", "gene_source", "gene_biotype"
        ]] = annotations_df["attribute"].apply(lambda attributes : separate_attributes(attributes))
        annotations_df = annotations_df[annotations_df["gene"] != None]
        annotations_df["gene"] = annotations_df["gene"].str.upper()
        annotations_df["gene"] = annotations_df["gene"].str.strip('"')
        annotations_df = annotations_df[annotations_df["gene"].notna()]

    except:

        raise

    return annotations_df

def separate_attributes(attributes):

    attributes = attributes[:-1].split("; ")
    gene_id = None
    gene_version = None
    gene_name = None
    gene_source = None
    gene_biotype = None

    for attribute in attributes:

        key, value = attribute.split(" ")

        if key == "gene_id" : gene_id = value
        if key == "gene_version" : gene_version = value
        if key == "gene_name" : gene_name = value
        if key == "gene_source" : gene_source = value
        if key == "gene_biotype" : gene_biotype = value

    return pd.Series([gene_id, gene_version, gene_name, gene_source, gene_biotype])

def read_expression():

    try:

        print("Reading genome sequence file.")

        expression_df = pd.read_csv(config["EXPRESSION_PATH"], sep = "\t")
        expression_df["gene"] = expression_df["gene"].str.upper()

    except:

        raise

    return expression_df

def combine_data(annotations_df, essentiality_df, expression_df):

    genes_df = pd.merge(annotations_df, essentiality_df, on = "gene")
    genes_df = pd.merge(genes_df, expression_df, on = "gene")

    return genes_df
