#!/bin/bash

# Set genome
GENOME="mm10"

# Directories (adjust if needed)
INPUT_DIR="./dars"
BACKGROUND_DIR="./backgrounds"
OUTPUT_DIR="./motif_results"

# Create output directory if not exists
mkdir -p $OUTPUT_DIR

# List of directions and file prefix pattern
for direction in up down; do
  for bed_file in ${INPUT_DIR}/${direction}_*.bed; do
    # Extract the base name (e.g., "Oligo_NN_2vs18" from up_Oligo_NN_2vs18.bed)
    base=$(basename "$bed_file" .bed)
    celltype=$(echo $base | sed "s/${direction}_//")

    # Define background file
    background_file="${BACKGROUND_DIR}/background_peaks_${celltype}.bed"

    # Define output directory
    out_dir="${OUTPUT_DIR}/${celltype}/${direction}"
    mkdir -p "$out_dir"

    # Run HOMER motif enrichment
    findMotifsGenome.pl "$bed_file" $GENOME "$out_dir" -bg "$background_file" -size given -p 8

    echo "Finished $direction for $celltype"
  done
done
