#!/usr/bin/env python

"""
This script splits male BAM files by cell type using barcode metadata.
Requires: pysam, pandas, numpy, multiprocessing.
"""

import os
import pandas as pd
import numpy as np
import pysam
from multiprocessing import Pool

# === Load metadata and BAM paths ===
meta = pd.read_csv("final_meta.csv", sep=",")
male_bams = pd.read_csv("male_files.csv", sep=",")
threads = 14

# Ensure sample naming is consistent
male_bams["Sample"] = "Male:" + male_bams["Sample_id"]
male_bams["path"] = male_bams["path"].str.replace(
    "/projects/ps-renlab/lamaral/",
    "/tscc/projects/ps-renlab/lamaral/",
    regex=False
)

# Extract barcodes from cell IDs in metadata
meta['barcode'] = meta['cell_id'].apply(lambda x: x.split(':')[2])

# === Get all unique cell types ===
celltypes = np.unique(meta['celltype_final'].astype(str))

# === Function to write BAMs by cell type ===
def write_bam_by_celltype(ct):
    print(f"Processing cell type: {ct}")
    filtered = meta[meta['celltype_final'] == ct]
    if filtered.shape[0] == 0:
        print(f"No entries for {ct}, skipping.")
        return

    sample_list = pd.unique(filtered['sample'])
    for sample in sample_list:
        subset = filtered[filtered['sample'] == sample]
        if subset.shape[0] == 0:
            continue

        barcode_list = list(subset['barcode'])
        output_dir = os.path.join("split_bams", ct.replace(" ", "_"))
        os.makedirs(output_dir, exist_ok=True)

        bam_path = male_bams[male_bams['Sample'] == sample]['path'].values[0]
        outbam_path = os.path.join(output_dir, sample.replace(":", "_") + ".bam")

        if os.path.exists(outbam_path):
            print(f"{outbam_path} exists, skipping.")
            continue

        print(f"Writing {outbam_path}")
        
        with pysam.AlignmentFile(bam_path, "rb") as infile, \
             pysam.AlignmentFile(outbam_path, "wb", header=infile.header) as out:
            for read in infile:
                bc = read.get_tag("BX") if read.has_tag("BX") else None
                if bc and bc in barcode_list:
                    out.write(read)

        # Index the BAM file
        pysam.index(outbam_path)

# === Run in parallel ===
if __name__ == "__main__":
    with Pool(threads) as p:
        p.map(write_bam_by_celltype, celltypes)