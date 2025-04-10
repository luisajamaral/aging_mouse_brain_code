#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
import sys

adata = sc.read_h5ad("AIT21_10Xv3.h5ad")

# Display the number of cells before filtering
num_cells_before_filtering = adata.shape[0]
print(f"Number of cells before filtering: {num_cells_before_filtering}")

# Step 1: Filter cells based on 'sex' and keep only "F"
filtered_sex = adata[adata.obs['sex'] == "F"]
num_cells_after_sex_filter = filtered_sex.shape[0]
print(f"Number of cells after filtering based on 'sex': {num_cells_after_sex_filter}")

# Step 2: Filter cells based on 'cl' and randomly sample up to 1000 cells from each cluster
max_cells_per_cluster = 1000
filtered_indices = []

cluster_indices = {}

for cluster in filtered_sex.obs['cl'].unique():
    cluster_indices[cluster] = np.where(filtered_sex.obs['cl'] == cluster)[0]
    if len(cluster_indices[cluster]) > max_cells_per_cluster:
        sampled_indices = np.random.choice(cluster_indices[cluster], max_cells_per_cluster, replace=False)
        filtered_indices.extend(sampled_indices)
    else:
        filtered_indices.extend(cluster_indices[cluster])

filtered_cl = filtered_sex[filtered_indices]
num_cells_after_cl_filter = filtered_cl.shape[0]
print(f"Number of cells after filtering based on 'cl' (up to 1000 per cluster): {num_cells_after_cl_filter}")

# Save the filtered AnnData to an H5AD file
output_file = "Female_10xV3_1000_per_cl.h5ad"
filtered_cl.write(output_file)

adata = sc.read_h5ad("Female_10xV3_1000_per_cl.h5ad")


# Read the key file
key_file = "AIT21_annotation_freeze_081523.tsv"
key_df = pd.read_csv(key_file, sep='\t')

# Merge the key information into the AnnData object based on the 'cl' column
adata.obs['subclass_label'] = adata.obs['cl'].map(key_df.set_index('cl')['subclass_label'])

max_cells_per_cluster = 1000
filtered_indices = []

cluster_indices = {}

for cluster in adata.obs['subclass_label'].unique():
    cluster_indices[cluster] = np.where(adata.obs['subclass_label'] == cluster)[0]
    if len(cluster_indices[cluster]) > max_cells_per_cluster:
        sampled_indices = np.random.choice(cluster_indices[cluster], max_cells_per_cluster, replace=False)
        filtered_indices.extend(sampled_indices)
    else:
        filtered_indices.extend(cluster_indices[cluster])

filtered_adata = adata[filtered_indices]
num_cells_after_cl_filter = filtered_adata.shape[0]
print(f"Number of cells after filtering based on 'cl' (up to 2000 per cluster): {num_cells_after_cl_filter}")






with open("AIT21_k8_markers.txt", "r") as file:
    gene_list = file.read().split("\n")

# Remove empty strings if they exist
gene_list = [gene for gene in gene_list if gene]

# Subset the AnnData to include only the genes in the list
filtered_cl_sub = filtered_cl[:, filtered_cl.var.index.isin(gene_list)]

output_file = "Female_10xV3_1000_per_cl_8k_gene.h5ad"
filtered_cl_sub.write(output_file)

max_cells_per_cluster = 100
filtered_indices = []

cluster_indices = {}

for cluster in filt.obs['cl'].unique():
    cluster_indices[cluster] = np.where(filt.obs['cl'] == cluster)[0]
    if len(cluster_indices[cluster]) > max_cells_per_cluster:
        sampled_indices = np.random.choice(cluster_indices[cluster], max_cells_per_cluster, replace=False)
        filtered_indices.extend(sampled_indices)
    else:
        filtered_indices.extend(cluster_indices[cluster])

filtered_cl = filtered_sex[filtered_indices]
num_cells_after_cl_filter = filtered_cl.shape[0]
print(f"Number of cells after filtering based on 'cl' (up to 1000 per cluster): {num_cells_after_cl_filter}")


# In[ ]:




