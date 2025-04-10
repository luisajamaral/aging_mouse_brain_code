#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import snapatac2 as snap
import scanpy as scs
import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import numpy as np


os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/')


# In[ ]:


def plot_and_save_doublet_scores(simulated_scores, observed_scores, filename):
    # Define score bins
    score_bins = np.arange(0, 1, 0.01)  # Adjust the range and bin width as needed

    # Create a histogram of doublet scores
    plt.figure(figsize=(10, 6))  # Set the figure size

    plt.hist(simulated_scores, bins=score_bins, edgecolor='black', color='skyblue', label='Simulated Doublet Scores', alpha=0.7)
    plt.hist(observed_scores, bins=score_bins, edgecolor='black', color='orange', label='Observed Doublet Scores', alpha=0.7)

    # Add labels and title
    plt.xlabel('Doublet Score')
    plt.ylabel('Count')
    plt.title('Doublet Score Distribution')
    plt.legend()

    # Customize x-axis tick marks
    plt.xticks(score_bins)
    plt.xticks(rotation=45)
    plt.gca().get_xaxis().set_major_locator(plt.MaxNLocator(nbins=25))  # Adjust the number of tick marks

    # Save the plot to the specified filename
    plt.tight_layout()  # Adjust spacing for labels
    plt.savefig(filename)  # Specify the filename and format
    plt.show()


# In[ ]:


track = pd.read_csv("tracksheet_with_paths_oct1723.csv", sep = ",") 
def preprocess_fragment(fragment_file, sample_name,nfrag,tsse,nfrag_max):
    
    frag_plot_file = ''.join(["female_h5ads/",sample_name, '.frag.pdf'])
    tsse_file = ''.join(["female_h5ads/",sample_name, '.tsse.pdf'])
    h5ad_file = ''.join(["female_h5ads/",sample_name, '.h5ad'])
    doublet_file = ''.join(["female_h5ads/",sample_name, '.doublet.png'])

    if not os.path.exists(h5ad_file):
        print(sample_name, "\n")
        file_path = ''.join(["female_h5ads/",sample_name, "_out.txt"])
        file = open(file_path, "w")

        data = snap.pp.import_data(
        fragment_file,
        chrom_sizes=snap.genome.mm10,
        file=''.join(["female_h5ads/",sample_name, ".h5ad"]),  # Optional
        sorted_by_barcode=False,chunk_size = 250000)
        snap.pl.frag_size_distr(data, interactive=False,out_file = frag_plot_file)
        snap.metrics.tsse(data, snap.genome.mm10)
        snap.pl.tsse(data, interactive=False, out_file = tsse_file)
        snap.pp.filter_cells(data, min_counts=nfrag,max_counts=nfrag_max, min_tsse = tsse)
        snap.pp.add_tile_matrix(data,bin_size=500)
        snap.pp.select_features(data)
        snap.pp.scrublet(data)
        file.write(str(data))
        plot_and_save_doublet_scores( data.uns["scrublet_sim_doublet_score"], data.obs["doublet_score"], doublet_file)

        #snap.pp.filter_doublets(data)
        #file.write(str(data))

        data.close()
    

for row in range(40,len(track)):
    fragment_file = "/" + track["Path-luisa"][row] + "/outs/atac_fragments.tsv.gz"
    sample_name = track["Sample_name"][row]
    nfrag = track["nfrag_cutoff"][row]
    tsse = track["TSS_cutoff"][row]
    nfrag_max = track["nfrag_max"][row]

    #print(fragment_file,sample_name, "\n")
    preprocess_fragment(fragment_file, sample_name,nfrag,tsse,nfrag_max)
    


# In[ ]:


data = snap.read("female_h5ads/AMY_2mo_1.h5ad")


# In[ ]:


data


# In[ ]:


doublet_file = ''.join(["female_h5ads/","AMY_2mo_1", '.doublet.png'])
plot_and_save_doublet_scores( data.uns["scrublet_sim_doublet_score"], data.obs["doublet_score"], doublet_file)


# In[ ]:


snap.pp.filter_doublets(data, score_threshold=0.08,probability_threshold=None)


# In[ ]:




