#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from gprofiler import gprofiler
import logging,matplotlib,os,sys

from rpy2.robjects import pandas2ri

pandas2ri.activate()
# anndata2ri.activate()
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')

plt.rcParams['figure.figsize']=(6,6) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=300, dpi_save=300)
sc.logging.print_versions()

matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10

os.getcwd()
os.chdir('/public/home/jphe/ATACseq/mm10/sc/E-MTAB-6714/bt2/scTE/scTE/scanpy/cellLines')
os.getcwd()


# In[ ]:


import glob, os, sys

adata = sc.read('../../c1_mESC.csv',cache=False)
sc.pp.filter_cells(adata, min_genes = 500)
adata.obs['time'] ='c1_mESC'

for f in sorted(glob.glob('../../*csv')):
    head = os.path.split(f)[1].replace('.csv','')
    if head == 'c1_mESC' or head == 'mSp':
        continue
    print(head)
    tmp = sc.read(f, cache = False)
    sc.pp.filter_cells(tmp, min_genes = 5)
    
    tmp.obs['time'] = head
    adata = adata.concatenate([tmp])

adata


# In[ ]:


adata.obs.index = [ k.split('-')[0] for k in adata.obs.index ]
adata.obs.head(5)


# In[ ]:


sc.pp.filter_genes(adata, min_cells=10)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
adata


# In[ ]:


adata.X = adata.X.astype('float64')
adata.obs['time'].value_counts()


# In[ ]:


adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes'] = (adata.X > 0).sum(1)

#Thresholding decision: counts
p3 = sb.distplot(adata.obs['n_counts'], kde=False)
plt.show()

#Thresholding decision: genes
p6 = sb.distplot(adata.obs['n_genes'], kde=False, bins=60)
plt.show()


# In[ ]:


p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes')
p2 = sc.pl.violin(adata, ['n_genes','n_counts'], groupby='time', size=2, log=True, cut=0)


# In[ ]:


sc.pp.filter_cells(adata, max_counts = 1e6)


# In[ ]:


adata.raw=adata


# In[ ]:


from __future__ import division
import glob,os

size_factors = {}
for f in glob.glob('../../../../*/zzresults.tsv'):
    head = os.path.split(f)[0]
    if 'mSp' in head:
            continue
    print(head)
    
    o = open(f,'rU')
    for l in o:
        if 'Sample' in l:
            continue
        l=l.strip().split('\t')
        sample=l[0].replace('mouse_skin_fibr','mSF')
        nf=int(l[1])/1e6
        
        size_factors[sample]=nf


# In[ ]:


size = [size_factors[k] for k in adata.obs.index ]


# In[ ]:


adata.X=adata.raw.X


# In[ ]:


adata.obs['size_factors'] = size

sc.pl.scatter(adata, 'size_factors', 'n_counts')
sc.pl.scatter(adata, 'size_factors', 'n_genes')

sb.distplot(size, bins=50, kde=False)
plt.show()

#Keep the count data in a counts layer
adata.layers["counts"] = adata.X.copy()

#Normalize adata 
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)

# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata


# In[ ]:


sc.pp.highly_variable_genes(adata, flavor='seurat')
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pl.highly_variable_genes(adata)

sc.pp.pca(adata, n_comps=20, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)

sc.tl.tsne(adata, n_jobs=20)
sc.tl.umap(adata)
sc.tl.diffmap(adata)

sc.pl.pca(adata,color=['time'])
sc.pl.tsne(adata,color=['time'])
sc.pl.umap(adata,color=['time'])


# In[ ]:


sc.pp.pca(adata, n_comps=20, use_highly_variable=False, svd_solver='arpack')
sc.pp.neighbors(adata)

sc.tl.tsne(adata, n_jobs=20)
sc.tl.umap(adata)

plt.rcParams['figure.figsize']=(6,6) 
sc.pl.pca(adata,color=['time'])
sc.pl.tsne(adata,color=['time'])
sc.pl.umap(adata,color=['time'])


# In[ ]:


adata = sc.read('write/scATAC.XC.cellLine.h5ad')


# In[ ]:


sc.set_figure_params(dpi=200, dpi_save=300)
plt.rcParams['figure.figsize']=(4,4) 
sc.pl.umap(adata,color=['time'], frameon=False, save='.sample.pdf',size=20, palette='tab10')


# In[ ]:


sc.set_figure_params(dpi=100, dpi_save=300)
sc.pl.umap(adata,color=['RLTR4_Mm'],frameon=True,vmin=3,vmax=6,save='.RLTR4_Mm.pdf')


# In[ ]:


sc.pl.umap(adata,color=['RMER19B'],frameon=False)


# In[ ]:


adata.obs.head()
cd4 = [k for k in adata.obs.index if adata.obs['leiden_r0.2'][k] in ['1']]
cd4 =[k for k in cells if 'cd4' in k]
len(cd4)

other =  [k for k in adata.obs.index if adata.obs['leiden_r0.2'][k] not in ['1']]
other[:5]

filt = cd4+other
filt = adata[filt,]


# In[ ]:


filt


# In[ ]:


sc.pp.pca(filt, n_comps=20, use_highly_variable=False, svd_solver='arpack')
sc.pp.neighbors(filt)

sc.tl.tsne(filt, n_jobs=20)
sc.tl.umap(filt)


# In[ ]:


sc.pl.umap(filt,color=['time'],frameon=False,size=20, palette='tab10',save='.sample.pdf')


# In[ ]:


adata = filt.copy()


# In[ ]:


sc.tl.leiden(adata, key_added='leiden_r1')
sc.tl.leiden(adata, resolution=0.2, key_added='leiden_r0.2')
sc.tl.leiden(adata, resolution=0.1, key_added='leiden_r0.1')


# In[ ]:


sc.pl.umap(adata,color=['time','leiden_r1','leiden_r0.2','leiden_r0.1'],save='.time.leiden.pdf', palette='Set1', frameon=False)


# In[ ]:


adata.obs.to_csv('write/obs.csv',sep=',')
adata.write('write/scATAC.XC.cellLine.h5ad')
adata


# In[ ]:





# In[ ]:


sc.tl.rank_genes_groups(adata,'leiden_r0.2', n_genes=100, method='t-test')

df=pd.DataFrame( {group + '_' + key[:1]: adata.uns['rank_genes_groups'][key][group]  for group in  set(adata.obs['leiden_r0.2']) for key in ['names','logfoldchanges']})
df.to_csv('write/rank.leiden_r0.2.tsv',sep='\t')

df.head(20)


# In[ ]:


plt.rcParams['figure.figsize']=(4,4) 
sc.pl.umap(adata,color=['RLTR13A','RLTR13A1','RLTR13A2','RLTR13A3','RLTR1F_Mm','RLTR1D'],save='CPC.TEs.pdf', frameon=True)


# In[ ]:


sc.pl.rank_genes_groups_heatmap(adata,n_genes=50,swap_axes=True, use_raw=False, standard_scale='var',dendrogram=False)


# In[ ]:


adataz=adata.copy()
sc.pp.scale(adataz)
sc.pl.rank_genes_groups_heatmap(adataz,n_genes=50,swap_axes=True, use_raw=False,dendrogram=False,vmax=1,vmin=-1, save='.zscore.pdf')


# In[ ]:


sc.pl.umap(adata, color =['RLTR13G','RLTR9E','RLTR9D',], save='.ESC.pdf', vmin=4, vmax=8)
sc.pl.umap(adata, color =['L1Md_T','L1Md_A'], save='.L1Md.pdf',vmin=5, vmax=10)


# In[ ]:


sc.pl.umap(adata, color =['RLTR13A','RLTR13A1','RLTR13A2','RLTR13A3'], save='.CPC.pdf',vmin=2, vmax=6)


# In[ ]:


sc.pl.umap(adata, color =['IAPLTR2_Mm','IAPLTR1_Mm','IAPLTR4','RMER6B','MER95'])


# In[ ]:


sc.pl.umap(adata, color =['RMER19B','RMER19C','MER121','MER125','RMER15'])


# In[ ]:




