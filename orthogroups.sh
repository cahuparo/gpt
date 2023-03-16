#!/usr/bin/env bash

# Load required modules
module load orthofinder/2.5.4

# Define input and output directories
input_dir="path/to/NLR_protein_sequences"
output_dir="path/to/orthogroups"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Run OrthoFinder
orthofinder \
    -f $input_dir \
    -t 4 \
    -S diamond \
    -A msa \
    -M mafft \
    -T fasttree \
    -o $output_dir

# Extract orthogroup families
awk '/^OG/{print ">"$0;next}{print}' $output_dir/OrthoFinder/Results*/Orthogroups.txt \
    | sed -e 's/ //g' -e 's/:/\n/g' \
    | awk '/^>/ {gsub("OG", ">Orthogroup"); print; next} {print}' \
    > $output_dir/orthogroups.fasta

