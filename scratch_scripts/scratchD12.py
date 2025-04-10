#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import snapatac2 as snap
import pandas as pd
import numpy as np
import os
import sys

os.chdir('/tscc/nfs/home/lamaral/ps-renlab2/projects/combined_all/h5ads_final')
data=snap.read_dataset("combined_save2-Copy1.h5ads")
adata = data.to_adata()
data.close()
adata_D12 = adata[adata.obs['celltype_final'] == "STR D12 Gaba"].copy()
snap.tl.macs3(adata_D12, n_jobs = 12, replicate = 'rep',replicate_qvalue=0.05, 
              qvalue=0.01,tempdir='/tscc/lustre/ddn/scratch/lamaral',
              groupby='celltype_batch_age',key_added = 'celltype_age_batch_01',
              blacklist="/tscc/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed")


# In[ ]:


import snapatac2 as snap
import pandas as pd
import numpy as np
import os
import sys

os.chdir('/tscc/nfs/home/lamaral/ps-renlab2/projects/combined_all/h5ads_final')


# In[ ]:


data=snap.read_dataset("combined_save2-Copy1.h5ads")
data
#data.close()


# In[ ]:


data.obsm


# In[ ]:


data.close()


# In[ ]:


adata


# In[ ]:


data.close()


# In[ ]:


adata_D12



# In[ ]:


import scanpy as sc


# In[ ]:


category_key = 'celltype_batch_age'
target_cells_per_category = 25000

# Get indices for each category
indices_to_keep = []
unique_categories = np.unique(adata_D12.obs[category_key])
for cat in unique_categories:
    cat_indices = np.where(adata_D12.obs[category_key] == cat)[0]
    np.random.shuffle(cat_indices)
    indices_to_keep.extend(cat_indices[:target_cells_per_category])


# In[ ]:


adata_D12_subsampled = adata_D12[indices_to_keep].copy()


# In[ ]:


adata_D12_subsampled.obs['celltype_batch_age'].value_counts()


# In[ ]:


adata_D12 = adata[adata.obs['celltype_final'] == "STR D12 Gaba"].copy()


# In[ ]:


adata_D12.obs['celltype_batch_age'].value_counts()


# In[ ]:


snap.tl.macs3(adata_D12, n_jobs = 12, replicate = 'rep',replicate_qvalue=0.05, 
              qvalue=0.01,tempdir='/tscc/lustre/ddn/scratch/lamaral',
              groupby='celltype_batch_age',key_added = 'celltype_age_batch_01',
              blacklist="/tscc/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed")


# In[ ]:


data.close()


# In[ ]:


adata.uns['celltype_final_rep_01']


# In[ ]:




