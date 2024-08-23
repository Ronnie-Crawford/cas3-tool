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
from gene_inspector import find_isolation, find_nearest_genes, find_size, find_search_windows, find_motifs_for_all_genes

# Read in and preprocess data
essentiality_df = read_essentiality()
annotations_df = read_genome_annotations()
expression_df = read_expression()

# Generate gene metrics
genes_df = combine_data(essentiality_df, annotations_df, expression_df)
genes_df = find_size(genes_df)
genes_df = find_isolation(genes_df)

# Find motifs
genes_df = find_search_windows(genes_df)
motifs_df = find_motifs_for_all_genes(genes_df)

print(genes_df[config["GENE_METRIC_COLUMNS"]])
genes_df.to_csv(f"{config["PATHS"]["RESULTS"]}/gene_metrics.tsv", sep = "\t", columns = config["GENE_METRIC_COLUMNS"])

print(motifs_df)
motifs_df.to_csv(f"{config["PATHS"]["RESULTS"]}/motif_metrics.tsv", sep = "\t")
# Sort genes based on desired metrics
#genes_df = genes_df[["gene", "seqname", "start", "end", "isolation", "number_motifs_upstream", "number_motifs_downstream"]].dropna()
#genes_df = genes_df.sort_values(by = 'isolation', ascending = False)
