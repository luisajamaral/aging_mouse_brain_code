#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import subprocess
import snapatac2 as snap
import scanpy as scs
import sys


# In[ ]:


data=snap.read_dataset("../h5ads_final/combined_save2-Copy1.h5ads")


# In[ ]:


data


# In[ ]:


peak_mat = snap.read("../h5ads_final/pmat_celltype_final_rep_01-Copy1.h5ad")


# In[ ]:


peaks = snap.tl.merge_peaks(data.uns['celltype_final_rep_01'], chrom_sizes=snap.genome.mm10)


# In[ ]:


# Define the min_pct threshold
min_pct = 0.01

# Add necessary columns to data.obs
data.obs['celltype_age'] = data.obs['celltype_final'] + ":" + data.obs['age']

# Get unique cell types
uniq_cts = data.obs['celltype_final'].unique()


# In[ ]:


# Define a function to process each cell type
def process_cell_type(ct):
    file = "diff_peaks_" + ct.replace(" ", "_") + "_2vs18.csv"
    background_file_csv = "background_peaks_" + ct.replace(" ", "_") + "_2vs18.csv"
    background_file_bed = "background_peaks_" + ct.replace(" ", "_") + "_2vs18.bed"
    bed_up_file = "up_peaks_" + ct.replace(" ", "_") + ".bed"
    bed_down_file = "down_peaks_" + ct.replace(" ", "_") + ".bed"
    homer_up_dir = bed_up_file.replace(".bed", "_motifs")
    homer_down_dir = bed_down_file.replace(".bed", "_motifs")
    anno_up_file = bed_up_file.replace(".bed", ".ann.txt")
    anno_down_file = bed_down_file.replace(".bed", ".ann.txt")

    # Check if all output files exist
    if (os.path.exists(file) and os.path.exists(background_file_csv) and os.path.exists(background_file_bed) and
        os.path.exists(bed_up_file) and os.path.exists(bed_down_file) and
        os.path.exists(homer_up_dir) and os.path.exists(homer_down_dir) and
        os.path.exists(anno_up_file) and os.path.exists(anno_down_file)):
        print(f"Skipping {ct}, all output files already exist.\n")
        return

    print(ct)
    print(file)
    group1 = ct + ":2mo"
    group2 = ct + ":18mo"
    g1 = data.obs['celltype_age'] == group1
    g2 = data.obs['celltype_age'] == group2
    
    if g1.sum() < 250 or g2.sum() < 250: 
        print("Not enough cells\n")
        return
    
    celltype = ct.split(":")[0]
    
    # Calculate presence of peaks in g1 and g2
    presence_g1 = (peak_mat.X[g1, :] > 0).sum(axis=0) / g1.sum()
    presence_g2 = (peak_mat.X[g2, :] > 0).sum(axis=0) / g2.sum()
    
    # Ensure the logical OR operation result is a flat boolean array
    background_mask = np.ravel((presence_g1 >= min_pct) | (presence_g2 >= min_pct))
    
    # Convert peak_mat.var_names to a NumPy array for proper indexing
    var_names_array = np.array(peak_mat.var_names)
    background_indices = np.where(background_mask)[0]
    background_peaks = var_names_array[background_indices]
    
    # Save background peaks to a DataFrame and then to a CSV and BED file
    background_peaks_df = pd.DataFrame(background_peaks, columns=['background_peaks'])
    background_peaks_df['pct_cells_2mo'] = np.ravel(presence_g1)[background_indices] * 100
    background_peaks_df['pct_cells_18mo'] = np.ravel(presence_g2)[background_indices] * 100
    background_peaks_df.to_csv(background_file_csv, index=False)
    
    # Convert "chr:start-end" to BED format
    bed_background_peaks = background_peaks_df['background_peaks'].str.split(r'[:-]', expand=True)
    bed_background_peaks.columns = ['chr', 'start', 'end']
    bed_background_peaks.to_csv(background_file_bed, sep="\t", header=False, index=False)
    print("DAR testing")
    # Perform differential accessibility testing
    diff_peaks = snap.tl.diff_test(
        peak_mat,
        cell_group1=g2,
        cell_group2=g1,
        min_log_fc=0.25, min_pct=min_pct,
        features=background_mask
    )
    
    if diff_peaks.shape[0] == 0:
        return
    column_names = diff_peaks.columns
    diff_peaks_df = pd.DataFrame(diff_peaks, columns=column_names)

    # Save differential peaks to CSV with percentage columns
    diff_peaks_df.to_csv(file, index=False)
    print("DONE")
    
    # Filter for significant peaks with adjusted p-value < 0.05
    significant_peaks = diff_peaks_df[diff_peaks_df['adjusted p-value'] < 0.05]
    
    # Separate up-regulated and down-regulated peaks
    up_peaks = significant_peaks[significant_peaks['log2(fold_change)'] > 0]
    down_peaks = significant_peaks[significant_peaks['log2(fold_change)'] < 0]
    
    # Convert "chr:start-end" to BED format for up and down peaks
    if up_peaks.empty:
        print(f"No significant up-regulated peaks for {ct}. Skipping up-regulated analysis.")
    else:
        # Convert "chr:start-end" to BED format for up peaks
        up_peaks_bed = up_peaks['feature name'].str.split(r'[:-]', expand=True)
        up_peaks_bed.columns = ['chr', 'start', 'end']
        up_peaks_bed.to_csv(bed_up_file, sep="\t", header=False, index=False)
        print(f"Up-regulated peaks saved to {bed_up_file}")

        # Run HOMER analysis for up-regulated peaks
        print("Homer up")
        run_homer_analysis(bed_up_file, background_file_bed)
        print("Homer up anno")
        run_homer_annotatepeaks(bed_up_file, background_file_bed)
    
    if down_peaks.empty:
        print(f"No significant down-regulated peaks for {ct}. Skipping down-regulated analysis.")
    else:
        # Convert "chr:start-end" to BED format for down peaks
        down_peaks_bed = down_peaks['feature name'].str.split(r'[:-]', expand=True)
        down_peaks_bed.columns = ['chr', 'start', 'end']
        down_peaks_bed.to_csv(bed_down_file, sep="\t", header=False, index=False)
        print(f"Down-regulated peaks saved to {bed_down_file}")

        # Run HOMER analysis for down-regulated peaks
        print("Homer down")
        run_homer_analysis(bed_down_file, background_file_bed)
        print("Homer down anno")
        run_homer_annotatepeaks(bed_down_file, background_file_bed)

def run_homer_analysis(bed_file, background_file):
    output_dir = bed_file.replace(".bed", "_motifs")
    command = f"findMotifsGenome.pl {bed_file} mm10 {output_dir} -size 200 -p 4 -bg {background_file} -nomotif"
    subprocess.run(command, shell=True, check=True)
    print(f"HOMER analysis completed for {bed_file}")

def run_homer_annotatepeaks(bed_file, background_file):
    output_dir = bed_file.replace(".bed", "")
    command = f"annotatePeaks.pl {bed_file} mm10 -go {output_dir}_go -genomeOntology {output_dir}_genOn -p 4 > {output_dir}.ann.txt"
    subprocess.run(command, shell=True, check=True)
    print(f"HOMER analysis completed for {bed_file}")
# Use ThreadPoolExecutor to parallelize the processing

#with ThreadPoolExecutor(max_workers=8) as executor:
#    executor.map(process_cell_type, uniq_cts)

#print("All tasks completed.")


# In[ ]:


for ct in uniq_cts: 
    print(ct)
    process_cell_type(ct)


# In[ ]:


data.close()
peak_mat.close()



# In[ ]:


uniq_cts


# In[ ]:




