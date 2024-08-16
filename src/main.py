# Cas3 Gene Ranker and Motif Finder
# For the Cas3 project in Parts Group, Wellcome Sanger Institute
# Made by Ronnie Crawford
# Started: 14th July 2024
# Last Modified: 12th August 2024

# To-Do:
# Finish README
# Separate gene ranking and motif finding
# Find guides from motifs
# Find exons for motifs
# Add data from other files

# Standard modules
import json

# Third-party modules
import numpy as np
import pandas as pd

# Local modules
from config_loader import config
from data_reader import read_genome_sequence, read_genome_annotations, read_essentiality, read_expression, combine_data
from gene_inspector import find_isolation, find_nearest_genes, find_search_windows, fetch_fasta_sequence, find_size

# Read in and preprocess data
essentiality_df = read_essentiality()
annotations_df = read_genome_annotations()
expression_df = read_expression()

# Generate gene metrics
genes_df = combine_data(essentiality_df, annotations_df, expression_df)
genes_df = find_size(genes_df)
genes_df = find_isolation(genes_df)

#
#genes_df = find_search_windows(genes_df)

print(genes_df)
genes_df.to_csv("./results/gene_metrics.tsv", sep = "\t")
# Sort genes based on desired metrics
#genes_df = genes_df[["gene", "seqname", "start", "end", "isolation", "number_motifs_upstream", "number_motifs_downstream"]].dropna()
#genes_df = genes_df.sort_values(by = 'isolation', ascending = False)
