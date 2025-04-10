#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import snapatac2 as snap
import scanpy as scs
import pandas as pd
import numpy as np
import os
import sys

#os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/')


# In[ ]:


import scipy


# In[ ]:


#os.chdir('h5ads_final')

data=snap.read_dataset("combined.h5ads")


# In[ ]:


data.obs['celltype_batch_age_region']= data.obs['celltype_batch_age']+"_"+data.obs['region']


# In[ ]:


agc= snap.tl.aggregate_X(adata, groupby="celltype_batch_age_region", file="celltype_batch_age_region_count.h5ad")


# In[ ]:


ag = snap.tl.aggregate_X(adata, groupby="celltype_batch_age_region", file="celltype_batch_age_region_RPM.h5ad", normalize = "RPM")


# In[ ]:


data.close()


# In[ ]:


data.close()


# In[ ]:


data.obs['celltype_batch_age_region'].head()


# In[ ]:


get_ipython().run_cell_magic('time', '', 'adata = data.to_adata()\nadata\n')


# In[ ]:


1133694


# In[ ]:


data.close()


# In[ ]:


peaks = snap.tl.merge_peaks(data.uns['celltype_final_rep_01'], chrom_sizes=snap.genome.mm10)


# In[ ]:


get_ipython().run_cell_magic('time', '', "peak_mat = snap.pp.make_peak_matrix(data, use_rep=peaks['Peaks'])\npeak_mat\n")


# In[ ]:


for unique_value in data.obs['celltype_final'].unique():
    p = peak_mat[peak_mat.obs['celltype_final'] == unique_value].copy()
    scipy.io.mmwrite(f"{unique_value}.mtx", p.X[:],field='integer')


# In[ ]:





# In[ ]:


data.close()


# In[ ]:


peak_mat


# In[ ]:


adata


# In[ ]:


p = peak_mat[adata.obs['celltype_final'] == unique_value].copy()


# In[ ]:


peak_mat


# In[ ]:


#adata = data.to_adata()

for unique_value in adata.obs['celltype_final'].unique():
    adata_c0 = adata[adata.obs['celltype_final'] == unique_value].copy()
    adata_c0
    # Perform your desired operations on adata_c0
    snap.pp.select_features(adata_c0, n_features=300000)
    snap.tl.spectral(adata_c0)
    snap.pp.harmony(adata_c0, batch="batch",  key_added='X_spectral_harmony')
    snap.pp.knn(adata_c0, use_rep="X_spectral_harmony")
    snap.tl.leiden(adata_c0, resolution = 1.5)
    snap.tl.umap(adata_c0, use_rep="X_spectral_harmony")
    # Save UMAP figures and HTML files
    snap.pl.umap(adata_c0, color="leiden", interactive=True, out_file=f"{unique_value}_subcluster_leiden.html")
    snap.pl.umap(adata_c0, color="age", interactive=True, out_file=f"{unique_value}_subcluster_age.html")
    snap.pl.umap(adata_c0, color="batch", interactive=False, out_file=f"{unique_value}_subcluster_batch.png")
    snap.pl.umap(adata_c0, color="region", interactive=True, out_file=f"{unique_value}_subcluster_region.html")
    snap.pl.umap(adata_c0, color="sample", interactive=True, out_file=f"{unique_value}_subcluster_sample.html")
    snap.pl.umap(adata_c0, color="tsse_max", interactive=False, out_file=f"{unique_value}_subcluster_tsse.png")
    snap.pl.umap(adata_c0, color="n_fragment_max", interactive=False, out_file=f"{unique_value}_subcluster_nfrag.png")
    snap.pl.umap(adata_c0, color="doublet_probability", interactive=False, out_file=f"{unique_value}_subcluster_doublet_prob.png")
    snap.pl.umap(adata_c0, color="CellType_1127", interactive=True, out_file=f"{unique_value}_CellType_1127.html")
    snap.pl.umap(adata_c0, color="celltypes_ext", interactive=True, out_file=f"{unique_value}_celltypes_ext.html")
    snap.pl.umap(adata_c0, color="predicted_id_sep", interactive=True, out_file=f"{unique_value}predicted_id_sep.html")

    # Create a DataFrame with metadata
    df = pd.DataFrame()
    df["cell_id"] = adata_c0.obs_names
    for element in ['sample', 'doublet_probability', 'tsse','n_fragment', 'age','leiden',
                    'L2Annot.rough', 'subclass_label_v3','best_celltype','region','batch','leiden_1',
                   'leiden-2', 'leiden-3', 'subclass_label_v3']:
        df[element] = adata_c0.obs[element].tolist()

    df['umap_x'] = adata_c0.obsm['X_umap'][:,0]
    df['umap_y'] = adata_c0.obsm['X_umap'][:,1]

    # Save the metadata table to CSV
    df.to_csv(f"{unique_value}_subcluster_meta.csv")

