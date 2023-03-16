import os
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import networkx as nx
import matplotlib.pyplot as plt

# Define the input file paths and output directories
genotype_dir = "/path/to/genotype_files/"
output_dir = "/path/to/output_dir/"

# Create the output directories if they don't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Define the threshold for sequence similarity
similarity_threshold = 0.8

# Define a function to read in the protein sequences from a FASTA file
def read_sequences(fasta_file):
    seq_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_records.append(record)
    return seq_records

# Define a function to perform a pairwise alignment of two protein sequences and return the similarity score
def get_similarity_score(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
    if len(alignments) > 0:
        best_alignment = alignments[0]
        alignment_length = max(len(best_alignment[0]), len(best_alignment[1]))
        similarity_score = float(best_alignment[2]) / float(alignment_length)
        return similarity_score
    else:
        return 0

# Read in the NLR protein sequences for all genotypes
sequences = {}
for genotype_file in os.listdir(genotype_dir):
    genotype_name = genotype_file.split(".")[0]
    genotype_path = os.path.join(genotype_dir, genotype_file)
    sequences[genotype_name] = read_sequences(genotype_path)

# Perform pairwise alignments of all NLR protein sequences
similarity_matrix = pd.DataFrame(index=sequences.keys(), columns=sequences.keys())
for genotype1 in sequences.keys():
    for genotype2 in sequences.keys():
        if genotype1 == genotype2:
            similarity_matrix.loc[genotype1, genotype2] = 1
        else:
            seq1 = sequences[genotype1]
            seq2 = sequences[genotype2]
            max_similarity_score = 0
            for nlr1 in seq1:
                for nlr2 in seq2:
                    similarity_score = get_similarity_score(nlr1, nlr2)
                    if similarity_score > max_similarity_score:
                        max_similarity_score = similarity_score
            similarity_matrix.loc[genotype1, genotype2] = max_similarity_score

# Generate a network graph of orthogroup families
G = nx.Graph()
for i, row in similarity_matrix.iterrows():
    for j, value in row.iteritems():
        if value >= similarity_threshold:
            G.add_edge(i, j)

# Write out the orthogroup families to a file
orthogroup_families_file = os.path.join(output_dir, "orthogroup_families.txt")
with open(orthogroup_families_file, "w") as outfile:
    for component in nx.connected_components(G):
        outfile.write("\t".join(component) + "\n")

# Generate a network graph of the orthogroups
pos = nx.spring_layout(G)
plt.figure(figsize=(12, 12))
nx.draw_networkx_nodes(G, pos, node_size=1000, node_color="lightblue")
nx.draw_networkx_edges(G, pos, width=1, alpha=0.7)
nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif")
plt.axis("off")
plt.savefig(os.path.join(output
