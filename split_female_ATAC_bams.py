#!/usr/bin/env python

"""
This script splits female BAM files by cell type using barcode metadata.
Requires: pysam, pandas, numpy, multiprocessing.
"""

import pandas as pd
import os
import subprocess

# Load the tracksheet with relevant sample info
track = pd.read_csv("tracksheet_with_paths_oct1723.csv")

# Iterate over each sample to split BAM files by barcode
for _, row in track.iterrows():
    path = row["Path-luisa"]
    sample_name = row["Sample_name"]
    bam_file = os.path.join(path, "outs", "possorted_bam.bam")
    barcodes_file = os.path.join(path, "outs", "filtered", "feature_bc_matrix", "barcodes.tsv.gz")

    print("Processing:", sample_name)

    # Output BAM file path
    output_bam = os.path.join("female_split_bams", f"{sample_name}_filtered.bam")

    # Run samtools to filter the BAM file using barcodes
    cmd = [
        "samtools", "view", "-h", "-b",
        f"-N", barcodes_file,
        "-o", output_bam,
        bam_file
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"Successfully created {output_bam}\n")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {sample_name}: {e}\n")
