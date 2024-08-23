# Cas3 Gene Ranker and Motif Finder
# For the Cas3 project in Parts Group, Wellcome Sanger Institute
# Made by Ronnie Crawford
# Started: 14th July 2024
# Last Modified: 12th August 2024

# To-Do:
# Test with genuine data
# Review isolation metric
# Identify strand of genes
# Fetch gene name of upstream and downstream genes

# Third-party modules
import numpy as np
import pandas as pd
from pyfaidx import Fasta

# Configuration
GENOME_SEQUENCE_PATH = '/Users/rc30/Documents/projects/cas3/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
GENOME_ANNOTATIONS_PATH = '/Users/rc30/Documents/projects/cas3/Homo_sapiens.GRCh38.112.gtf'
ESSENTIALITY_PATH = "/Users/rc30/Documents/projects/cas3/HAP1_essentiality_Toronto (1).tsv"

CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
PAM_MOTIFS = ["AAG"]

GENE_ESSENTIALITY_THRESHOLD = -2
FILTER_GENES_CUTOFF = 100

SEARCH_WINDOW_UPSTREAM_SIZE_MAX = 100000
SEARCH_WINDOW_DOWNSTREAM_SIZE_MAX = 100000
GUIDE_LENGTH = 30

def find_nearest_genes(gene, chromosome_df):
    
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
        
        nearest_upstream_start = upstream_candidates_df.iloc[-1]['start']
        nearest_upstream_end = upstream_candidates_df.iloc[-1]['end']
      
    if not downstream_candidates_df.empty:

        nearest_downstream_start = downstream_candidates_df.iloc[0]['start']
        nearest_downstream_end = downstream_candidates_df.iloc[0]['end']
        
    return pd.Series([nearest_upstream_name, nearest_upstream_start, nearest_upstream_end, nearest_downstream_name, nearest_downstream_start, nearest_downstream_end])
    
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

def fetch_fasta_sequence(genome, chromosome, start, end):
    
    chromosome = f"{chromosome}"
    start = int(start)
    end = int(end)
    
    searched_sequence = genome[chromosome][start:end].seq
    
    number_of_motifs = 0
    
    for motif in PAM_MOTIFS:
        
        number_of_motifs += searched_sequence.count(motif)
    
    return int(number_of_motifs)

def find_search_windows(genes_df):
    
    genome = Fasta(GENOME_SEQUENCE_PATH)
    
    for chromosome in CHROMOSOMES:
        
        print(f"Finding FASTA sequence around genes in chromosome: {chromosome}.")
        
        chromosome_df = genes_df[genes_df["seqname"] == chromosome]
        
        if chromosome_df.empty:
        
            print(f"No essential genes found on chromosome {chromosome}")
            continue
        
        chromosome_df = chromosome_df.sort_values(by = 'start').reset_index(drop = True)
        
        chromosome_df["upstream_window_start"] = np.where(
            (chromosome_df["start"] - chromosome_df["nearest_upstream_end"]) < SEARCH_WINDOW_UPSTREAM_SIZE_MAX,
            chromosome_df["nearest_upstream_end"],
            chromosome_df["start"] - SEARCH_WINDOW_UPSTREAM_SIZE_MAX
            )
        chromosome_df["upstream_window_end"] = chromosome_df["start"] - 1
        chromosome_df["downstream_window_start"] = chromosome_df["end"] + 1
        chromosome_df["downstream_window_end"] = np.where(
            (chromosome_df["nearest_downstream_start"] - chromosome_df["end"]) < SEARCH_WINDOW_DOWNSTREAM_SIZE_MAX,
            chromosome_df["nearest_downstream_start"],
            chromosome_df["end"] + SEARCH_WINDOW_DOWNSTREAM_SIZE_MAX
            )
        chromosome_df["number_motifs_upstream"] = chromosome_df.apply(lambda gene : fetch_fasta_sequence(genome, chromosome, gene["upstream_window_start"], gene["upstream_window_end"]), axis = 1)
        chromosome_df["number_motifs_downstream"] = chromosome_df.apply(lambda gene : fetch_fasta_sequence(genome, chromosome, gene["downstream_window_start"], gene["downstream_window_end"]), axis = 1)

        chromosome_dfs.append(chromosome_df)
    
    genes_df = pd.concat(chromosome_dfs, ignore_index=True)
    
    return genes_df

# Get list of essential genes, with threshold
essentiality_df = pd.read_csv(ESSENTIALITY_PATH, sep = '\t')
essentiality_df = essentiality_df[essentiality_df["lfc"] < GENE_ESSENTIALITY_THRESHOLD]
essentiality_df["gene"] = essentiality_df["gene"].str.upper()

# Get list of all genes, chromosome, start and end
annotation_column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annotation_column_dtypes = dtype = {'seqname': str, 'source': str, 'feature': str, 'start': int, 'end': int, 'score': str, 'strand': str, 'frame': str, 'attribute': str}
annotations_df = pd.read_csv(GENOME_ANNOTATIONS_PATH, sep = '\t', comment = '#', header = None, names = annotation_column_names, dtype = annotation_column_dtypes)
annotations_df = annotations_df[annotations_df["feature"] == "gene"]
annotations_df = annotations_df[annotations_df['seqname'].isin(CHROMOSOMES)]
annotations_df[["gene_id", "gene_version", "gene", "gene_source", "gene_biotype"]] = annotations_df["attribute"].apply(lambda attributes : separate_attributes(attributes))
annotations_df = annotations_df[annotations_df["gene"] != None]
annotations_df["gene"] = annotations_df["gene"].str.upper()
annotations_df["gene"] = annotations_df["gene"].str.strip('"')
annotations_df = annotations_df[annotations_df["gene"].notna()]

# Combine data
genes_df = pd.merge(annotations_df, essentiality_df, on = "gene")
print("Number of essential genes: ", len(essentiality_df))
print("Number of annotated genes: ", len(annotations_df[annotations_df["gene"].notna()]))
print("How many genes are shared?: ", len(essentiality_df["gene"][essentiality_df["gene"].isin(annotations_df["gene"])]), "\n")

if genes_df.empty:
    
    print("No essential genes found in annotations file, ignoring essentiality.")
    genes_df = annotations_df

# Find isolation of genes
chromosome_dfs = []

for chromosome in CHROMOSOMES:
    
    print(f"Finding attributes of genes on chromosome : {chromosome}.")
    
    essential_chromosome_df = genes_df[genes_df["seqname"] == chromosome]
    annotations_chromosome_df = annotations_df[annotations_df["seqname"] == chromosome]
    
    if essential_chromosome_df.empty:
        
        print(f"No essential genes found on chromosome {chromosome}")
        continue
    
    annotations_chromosome_df = annotations_chromosome_df.sort_values(by = 'start').reset_index(drop = True)
    essential_chromosome_df[[
        'nearest_upstream_name', 'nearest_upstream_start', 'nearest_upstream_end',
        'nearest_downstream_name', 'nearest_downstream_start', 'nearest_downstream_end'
        ]] = essential_chromosome_df.apply(lambda gene : find_nearest_genes(gene, annotations_chromosome_df), axis = 1)

    essential_chromosome_df["distance_to_upstream"] = essential_chromosome_df["start"] - essential_chromosome_df["nearest_upstream_end"]
    essential_chromosome_df["distance_to_downstream"] = essential_chromosome_df["nearest_downstream_start"] - essential_chromosome_df["end"]
    essential_chromosome_df.loc[:, "min_distance_to_gene"] = essential_chromosome_df[['distance_to_upstream', 'distance_to_downstream']].min(axis=1)
    
    chromosome_dfs.append(essential_chromosome_df)
    
genes_df = pd.concat(chromosome_dfs, ignore_index=True) 

# Get sequence near each gene
genes_df = find_search_windows(genes_df)

# Sort genes based on desired metrics
genes_df = genes_df[["gene", "seqname", "start", "end", "min_distance_to_gene", "number_motifs_upstream", "number_motifs_downstream"]].dropna()
genes_df = genes_df.sort_values(by = 'min_distance_to_gene', ascending = False).reset_index(drop=True)

print("Number of essential genes: ", len(essentiality_df))
print("Number of annotated genes: ", len(annotations_df[annotations_df["gene"].notna()]))
print("How many genes are shared?: ", len(essentiality_df["gene"][essentiality_df["gene"].isin(annotations_df["gene"])]), "\n")


print(genes_df)
genes_df.to_csv("ranked_essential_genes.tsv", sep="\t")