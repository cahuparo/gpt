import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Define input and output files
input_dir = "path/to/input/directory"
output_dir = "path/to/output/directory"
genotype_list = "path/to/genotype_list.txt"
nlr_seq_file = "path/to/NLR_sequences.fasta"
evalue_cutoff = 1e-10

# Read genotype list
with open(genotype_list, "r") as f:
    genotypes = [line.strip() for line in f]

# Read NLR protein sequences
nlr_seqs = SeqIO.parse(nlr_seq_file, "fasta")

# Initialize dictionary to store NLR sequences by genotype
genotype_nlr_seqs = {}
for genotype in genotypes:
    genotype_nlr_seqs[genotype] = []

# Store NLR sequences in dictionary by genotype
for nlr_seq in nlr_seqs:
    for genotype in genotypes:
        if genotype in nlr_seq.id:
            genotype_nlr_seqs[genotype].append(nlr_seq)

# Perform all-vs-all BLASTP
for i, genotype1 in enumerate(genotypes):
    for j, genotype2 in enumerate(genotypes[i+1:]):
        blastp_output = f"{output_dir}/{genotype1}_vs_{genotype2}.blastp"
        if not os.path.isfile(blastp_output):
            cmd = f"blastp -query {input_dir}/{genotype1}.fasta -subject {input_dir}/{genotype2}.fasta -evalue {evalue_cutoff} -out {blastp_output}"
            os.system(cmd)

# Parse BLASTP output and store results in DataFrame
df = pd.DataFrame(columns=["Genotype1", "Genotype2", "E-value", "Orthogroup"])
for i, genotype1 in enumerate(genotypes):
    for j, genotype2 in enumerate(genotypes[i+1:]):
        blastp_output = f"{output_dir}/{genotype1}_vs_{genotype2}.blastp"
        with open(blastp_output, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                evalue = float(fields[-2])
                if evalue <= evalue_cutoff:
                    seq1 = SeqRecord(Seq(fields[0]), id=genotype1, description="")
                    seq2 = SeqRecord(Seq(fields[1]), id=genotype2, description="")
                    orthogroup = f"{genotype1}_{genotype2}"
                    df = df.append({"Genotype1": genotype1, "Genotype2": genotype2, "E-value": evalue, "Orthogroup": orthogroup}, ignore_index=True)

# Define function to find orthogroups
def find_orthogroups(df, genotype_nlr_seqs, orthogroups=None):
    if orthogroups is None:
        orthogroups = {}
    for genotype in genotypes:
        if genotype not in orthogroups:
            orthogroups[genotype] = []
    for i, row in df.iterrows():
        seq1 = SeqRecord(Seq(row["Orthogroup"].split("_")[0]), id=row["Genotype1"], description="")
        seq2 = SeqRecord(Seq(row["Orthogroup"].split("_")[1]), id=row["Genotype2"], description="")
        if seq1 in genotype_nlr_seqs[row["Genotype1"]] and seq2 in genotype_nlr_seqs[row["
