#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import snapatac2 as snap
import scanpy as scs
import pandas as pd
import numpy as np
import os
import sys

os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/')


# In[ ]:


print("IN")


# ### Convert male bam files to fragment files + h5ad

# In[ ]:


track = pd.read_csv("male_files.csv", sep = ",")
def preprocess_fragment(bam_file, sample_name):
    frag_file =''.join(["male_h5ads/",sample_name, '.frag'])
    h5ad_file = ''.join(["male_h5ads/",sample_name, '.h5ad'])
    tsse_file = ''.join(["male_h5ads/",sample_name, '.tsse.pdf'])
    frag_plot_file = ''.join(["male_h5ads/",sample_name, '.frag.pdf'])
    if not os.path.exists(frag_file):
        data = snap.pp.make_fragment_file(
            bam_file, barcode_tag = 'BX',
            output_file=frag_file, chunk_size = 50000000)
        print(sample_name,"created\n")
        data
    if not os.path.exists(h5ad_file):
        file_path = ''.join(["male_h5ads/",sample_name, "_out.txt"])
        file = open(file_path, "w")

        print(sample_name,"h5ad\n")

        data = snap.pp.import_data(frag_file,
                   chrom_sizes=snap.genome.mm10,
        file=h5ad_file, chunk_size = 250000,
        sorted_by_barcode=False)
        snap.pl.frag_size_distr(data, interactive=False,out_file = frag_plot_file)
        snap.metrics.tsse(data, snap.genome.mm10)
        snap.pl.tsse(data, interactive=False, out_file = tsse_file)
        file.write(str(data))

        snap.pp.filter_cells(data, min_counts=1000,max_counts=100000, min_tsse = 7)
        snap.pp.add_tile_matrix(data,bin_size=500)
        snap.pp.select_features(data)
        snap.pp.scrublet(data)
        file.write(str(data))

        snap.pp.filter_doublets(data)
        file.write(str(data))

        data.close()
        
for row in range(len(track)):
    bam_file = track["path"][row] 
    sample_name = track["Sample_id"][row]
    print(sample_name, "\n")
    preprocess_fragment(bam_file, sample_name)


# In[ ]:





# ### Convert Female fragments to h5ad

# In[ ]:


track = pd.read_csv("tracksheet_with_paths_oct1723.csv", sep = ",") 
def preprocess_fragment(fragment_file, sample_name,nfrag,tsse,nfrag_max):
    
    frag_plot_file = ''.join(["female_h5ads_filt_oct1723/",sample_name, '.frag.pdf'])
    tsse_file = ''.join(["female_h5ads_filt_oct1723/",sample_name, '.tsse.pdf'])
    h5ad_file = ''.join(["female_h5ads_filt_oct1723/",sample_name, '.h5ad'])

    if not os.path.exists(h5ad_file):
        print(sample_name, "\n")
        file_path = ''.join(["female_h5ads_filt_oct1723/",sample_name, "_out.txt"])
        file = open(file_path, "w")

        data = snap.pp.import_data(
        fragment_file,
        chrom_sizes=snap.genome.mm10,
        file=''.join(["female_h5ads_filt_oct1723/",sample_name, ".h5ad"]),  # Optional
        sorted_by_barcode=False,chunk_size = 250000)
        snap.pl.frag_size_distr(data, interactive=False,out_file = frag_plot_file)
        snap.metrics.tsse(data, snap.genome.mm10)
        snap.pl.tsse(data, interactive=False, out_file = tsse_file)
        snap.pp.filter_cells(data, min_counts=nfrag,max_counts=nfrag_max, min_tsse = tsse)
        snap.pp.add_tile_matrix(data,bin_size=500)
        snap.pp.select_features(data)
        snap.pp.scrublet(data)
        file.write(str(data))
        snap.pp.filter_doublets(data)
        file.write(str(data))

        data.close()
    

for row in range(len(track)):
    fragment_file = "/" + track["Path-luisa"][row] + "/outs/atac_fragments.tsv.gz"
    sample_name = track["Sample_name"][row]
    nfrag = track["nfrag_cutoff"][row]
    tsse = track["TSS_cutoff"][row]
    nfrag_max = track["nfrag_max"][row]

    #print(fragment_file,sample_name, "\n")
    preprocess_fragment(fragment_file, sample_name,nfrag,tsse,nfrag_max)
    


# ### Creating anndataset
# #### Male TSSe 5 , frag 1000
# #### Female TSSe 5, frag 1000

# In[ ]:





# In[ ]:


from pathlib import Path
from anndata import AnnData

directory = "h5ads_filter/anndatas/" 
adatas = []

for filename in os.listdir(directory):
    if filename.endswith(".h5ad"):
        key = os.path.splitext(filename)[0]
        file_path = os.path.join(directory, filename)
        adatas.append((key, Path(file_path)))


# In[ ]:


from pathlib import Path
from anndata import AnnData

directory = "male_h5ads/" 

adatas = []

for filename in os.listdir(directory):
    if filename.endswith(".h5ad"):
        key = "Male:"+os.path.splitext(filename)[0]
        file_path = os.path.join(directory, filename)
        adatas.append((key, Path(file_path)))

#directory = "/home/lamaral/scratch/combined/female/" 
directory = "female_h5ads_filt_oct1723/" 

for filename in os.listdir(directory):
    if filename.endswith(".h5ad"):
        key = "Female:"+os.path.splitext(filename)[0]
        file_path = os.path.join(directory, filename)
        adatas.append((key, Path(file_path)))


# In[ ]:


get_ipython().run_cell_magic('time', '', 'data = snap.AnnDataSet(\n    adatas=adatas,\n    filename="combined.h5ads"\n)\ndata\ndata.close()\n')


# In[ ]:


data= snap.read_dataset("combined.h5ads")
sample_names = np.array(data.obs['sample']).astype(str)
obs_names = np.array(data.obs_names).astype(str)
result = np.char.add(sample_names, ":")
result= np.char.add(result, obs_names)
data.obs_names = result
data.close()


# In[ ]:


data = snap.read_dataset("combined.h5ads")


# ### Select features, run spectral, plot UMAP (slow)

# In[ ]:


get_ipython().run_cell_magic('time', '', 'data= snap.read_dataset("combined.h5ads")\nsnap.pp.add_tile_matrix(data,bin_size=500)\ndata.close()\n')


# In[ ]:


get_ipython().run_cell_magic('time', '', 'data= snap.read_dataset("combined.h5ads")\nsnap.pp.select_features(data, 700000)\ndata.close()\n')


# In[ ]:


get_ipython().run_cell_magic('time', '', 'data= snap.read_dataset("combined.h5ads")\nsnap.tl.spectral(data)\ndata.close()\n')


# In[ ]:


data= snap.read_dataset("combined.h5ads")
data


# In[ ]:


get_ipython().run_cell_magic('time', '', 'data= snap.read_dataset("combined.h5ads")\nsnap.tl.umap(data)\nsnap.pl.umap(data, color="sample", interactive=True, out_file = "combined_before_sample.html", width =700, height = 550)\ndata.close()\n')


# In[ ]:


data


# ### Run batch effect correction + cluster

# In[ ]:


get_ipython().run_cell_magic('time', '', 'data=snap.read_dataset("combined.h5ads")\ndata.obs[\'batch\'] = [item.split(\':\')[0] for item in data.obs_names]\nsnap.pl.umap(data,color="batch", interactive=True, out_file = "combined_before_batch.html", width =700, height = 550)\n\n# mnc\nsnap.pp.mnc_correct(data, batch="batch")\nsnap.tl.umap(data, use_rep="X_spectral_mnn")\nsnap.pl.umap(data, color="sample", interactive=True, out_file = "combined_mnc.html", width =700, height = 550)\nsnap.pl.umap(data, color="batch", interactive=True, out_file = "combined_mnc_batch.html", width =700, height = 550)\n\ndata.close()\n')


# In[ ]:


data= snap.read_dataset("combined.h5ads")
snap.pp.harmony(data, batch = "batch")
snap.tl.umap(data, use_rep="X_spectral_harmony")
snap.pl.umap(data, color="sample", interactive=True, out_file = "combined_harmony.html", width =700, height = 550)
snap.pl.umap(data, color="batch", interactive=True, out_file = "combined_harmony_batch.html", width =700, height = 550)
data.close()


# In[ ]:


data=snap.read_dataset("combined.h5ads")
snap.pp.knn(data, use_rep="X_spectral_harmony")
snap.tl.leiden(data)
snap.pl.umap(data, color="leiden", interactive=True, out_file = "combined_harmony_leiden.html", width =700, height = 550)
data.obs['leiden_1'] = data.obs['leiden']
data.close()


# In[ ]:


data=snap.read_dataset("combined.h5ads")
snap.tl.leiden(data, resolution = 2)
data.obs['leiden-2'] = data.obs['leiden']
snap.pl.umap(data, color="leiden", interactive=True, out_file = "combined_harmony_leiden_2.html", width =700, height = 550)
data.close()


# In[ ]:


data=snap.read_dataset("combined.h5ads")
snap.tl.leiden(data, resolution = 3)
data.obs['leiden-3'] = data.obs['leiden']
snap.pl.umap(data, color="leiden", interactive=True, out_file = "combined_harmony_leiden_3.html", width =700, height = 550)
data.close()


# In[ ]:


data


# ### Annotate with Male data

# In[ ]:


data=snap.read_dataset("combined_save.h5ads")


# In[ ]:


meta = pd.read_csv("meta_with_rachel_anno.txt", sep = "\t")


# In[ ]:


meta.index


# In[ ]:


meta['rachel_celltypes'].fillna('none', inplace=True)

data.obs['rachel_celltypes'] = meta['rachel_celltypes']


# In[ ]:


meta['L2Annot.rough'].fillna('none', inplace=True)

data.obs['L2Annot.rough'] = meta['L2Annot.rough']


# In[ ]:


snap.pl.umap(data, color='L2Annot.rough', interactive=True, width =800, height = 550)


# In[ ]:


snap.pl.umap(data, color='L2Annot.rough', interactive=True, width =800, height = 550, out_file = "L2Annot.rough.html")


# In[ ]:


snap.pl.umap(data, color='rachel_celltypes', interactive=True, width =800, height = 550, out_file = "rachel_celltypes.html")


# In[ ]:


data.close()


# In[ ]:


meta['subclass_label_v3'].fillna('none', inplace=True)

meta['subclass_label_v3']


# In[ ]:


data.obs['subclass_label_v3']= meta['subclass_label_v3']


# In[ ]:


#data=snap.read_dataset("combined.h5ads")
meta = pd.read_csv("/home/lamaral/projects/brain_aging_mouse/analysis/integrtaion_adata/meta_after.txt", sep = "\t")
track = pd.read_csv("male_files.csv", sep = ",")
track['sample'] = track['path'].str.split('/').str[-1]
track['sample'] = track['sample'].str.split('.').str[0]

# Merge 'meta' and 'track' DataFrames based on 'sample' column
merged_df = meta.merge(track[['sample', 'Sample_id']], on='sample', how='left')

# Rename the merged column to 'sample_id' in 'meta' DataFrame
meta['Sample_id'] = merged_df['Sample_id']
meta['id'] = "Male:"+ meta['Sample_id'] + ':' + meta['barcode']

id_df = pd.DataFrame({'id': data.obs_names})

# Merge the 'table' DataFrame with the 'id_df' DataFrame based on the 'id' column
merged_df = pd.merge(meta, id_df, on='id', how='right' )

# Get the resulting clusters
filtered_clusters = merged_df['Before']


# In[ ]:


snap.pl.umap(data, color=filtered_clusters, interactive=True, out_file = "combine_ct_annotation.html", width =700, height = 550)


# In[ ]:


data.obs['frac_mito'] = pd.to_numeric(data.adatas.obs['frac_mito'])
data.obs['n_fragment'] = pd.to_numeric(data.adatas.obs['n_fragment'])
data.obs['log_n_fragment'] = np.log(pd.to_numeric(data.adatas.obs['n_fragment']))
data.obs['doublet_probability'] = pd.to_numeric(data.adatas.obs['doublet_probability'])
data.obs['tsse'] = pd.to_numeric(data.adatas.obs['tsse'])
data.obs['n_fragment_max'] = data.obs['n_fragment']
data.obs['n_fragment_max'] = data.obs['n_fragment_max'].apply(lambda x: 15000 if x > 15000 else x)
data.obs['tsse_max'] = data.obs['tsse'].apply(lambda x: 20 if x > 20 else x)


# In[ ]:


snap.pl.umap(data, color='tsse_max', interactive=False, out_file = "combined_tsse_max.html", width =700, height = 550)


# In[ ]:


snap.pl.umap(data, color='n_fragment_max', interactive=False, out_file = "combined_nfrag_15kmax.html", width =700, height = 550)
snap.pl.umap(data, color='tsse', interactive=False, out_file = "combined_tsse.html", width =700, height = 550)
snap.pl.umap(data, color='doublet_probability', interactive=False, out_file = "combined_doublet_probability.html", width =700, height = 550)


# In[ ]:


data.close()


# In[ ]:


data=snap.read_dataset("combined.h5ads")

df = pd.DataFrame()
df["cell_id"] = data.obs_names
for element in ['sample', 'doublet_probability', 'n_fragment']:
    df[element] = data.obs[element]


for element in ['tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_probability', 'doublet_score']:
    df[element] = data.adatas.obs[element]

df['male_ct'] = filtered_clusters
df['batch'] = data.obs['batch']
df['leiden_1'] = data.obs['leiden_1']
df['leiden_2'] = data.obs['leiden-2']
df['leiden_3'] = data.obs['leiden-3']

df['umap_x'] = data.obsm['X_umap'][:,0]
df['umap_y'] = data.obsm['X_umap'][:,1]
df.to_csv("combined_meta.csv") 
data.close()


# In[ ]:


data=snap.read_dataset("combined.h5ads")


# In[ ]:


meta = pd.read_csv("meta.csv")


# In[ ]:


(meta["keep"]=="yes")


# In[ ]:


data.subset(obs_indices = (meta["keep"]=="yes"),out="h5ads_filter")


# In[ ]:


data


# In[ ]:


data.close()


# In[ ]:


print("out")


# In[ ]:


data=snap.read_dataset("combined.h5ads")
meta = pd.read_csv("meta.csv", sep = ",")


# In[ ]:


data=snap.read_dataset("combined.h5ads")


# In[ ]:


snap.pl.umap(data, color=meta['best_celltype'], interactive=True, out_file = "combine_best_annotation.html", width =700, height = 550)


# In[ ]:


snap.pl.umap(data, color=meta['best_celltype'], interactive=True, width =700, height = 550)


# In[ ]:


snap.pl.umap(data, color=meta['age'], interactive=True, width =700, height = 550)


# In[ ]:


data.obs['age'] = meta['age']
data.obs['best_celltype'] = meta['best_celltype']
data.obs['region'] = meta['region']


# In[ ]:


data.close()


# ### Subclustering

# In[ ]:


data=snap.read_dataset("combined_save-Copy3.h5ads")


# In[ ]:


#data.subset(obs_indices = (meta["doub"]=="no"))


# In[ ]:


os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/subclustering3')


# In[ ]:


snap.pl.umap(data, color='leiden-.25', interactive=True, 
             width =700, height = 550)


# In[ ]:


get_ipython().run_cell_magic('time', '', 'adata = data.to_adata()\nadata\n')


# In[ ]:


data.close()


# In[ ]:


for unique_value in adata.obs['leiden-.25'].unique():
    adata_c0 = adata[adata.obs['leiden-.25'] == unique_value].copy()
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
    snap.pl.umap(adata_c0, color="best_celltype", interactive=True, out_file=f"{unique_value}_subcluster_best_celltype.html")
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
    df.to_csv(f"{unique_value}_combined_meta.csv")


# In[ ]:


data = snap.read_dataset("combined_save-Copy3.h5ads")


# In[ ]:


meta = pd.read_csv("final_meta.csv", sep = ",")


# In[ ]:


meta.head()


# In[ ]:


data


# In[ ]:





# In[ ]:


data.close()


# In[ ]:


#data.obs['best_subclass_label_V3'] = meta['best_subclass_label_V3']
#data.obs['best_L2Annot.rough'] = meta['best_L2Annot.rough']
#data.obs['best_celltype_fixed2'] = meta['best_celltype_fixed']
#meta['celltypes_ext'].fillna('none', inplace=True)
#data.obs['celltypes_ext'] = meta['celltypes_ext']

#meta['celltype'].fillna('none', inplace=True)
#data.obs['celltype_RNA'] = meta['celltype']



#meta['CellType_1127'].fillna('none', inplace=True)
#data.obs['CellType_1127'] = meta['CellType_1127']




#meta['predicted_id_sep'].fillna('none', inplace=True)
#data.obs['predicted_id_sep'] = meta['predicted_id_sep']


meta['celltype_final'].fillna('none', inplace=True)
data.obs['celltype_final'] = meta['celltype_final']


# In[ ]:


data.close()


# ### Call Peaks + DARs

# In[ ]:


data=snap.read_dataset("combined.h5ads")


# In[ ]:


data = data.subset(obs_indices = (data.obs["celltype_final"]!="doublet"),out="h5ads_final",)


# In[ ]:


data.close()


# In[ ]:


data=snap.read_dataset("h5ads_final/combined.h5ads")


# In[ ]:


snap.pl.umap(data, color='celltype_final', interactive=True, 
             width =700, height = 550)


# In[ ]:


os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/h5ads_final')
data=snap.read_dataset("combined.h5ads")


# In[ ]:


snap.pl.umap(data, color="sample", interactive=True, out_file = "combined_harmony.html", width =700, height = 550)
snap.pl.umap(data, color="batch", interactive=True, out_file = "combined_harmony_batch.html", width =700, height = 550)
snap.pl.umap(data, color='celltype_final', interactive=True, 
             width =700, height = 550,out_file = "celltype_final.html",)


# In[ ]:


print("in")


# In[ ]:


snap.pl.umap(data, color='celltype_final', interactive=True, 
             width =700, height = 550)


# In[ ]:


os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/h5ads_final')
data=snap.read_dataset("combined.h5ads")
snap.pp.harmony(data, batch = "batch")
snap.tl.umap(data, use_rep="X_spectral_harmony")
snap.pl.umap(data, color="sample", interactive=True, out_file = "combined_harmony.html", width =700, height = 550)
snap.pl.umap(data, color="batch", interactive=True, out_file = "combined_harmony_batch.html", width =700, height = 550)
data.close()


# In[ ]:


data=snap.read_dataset("combined.h5ads")


# In[ ]:


data


# In[ ]:





# In[ ]:


print("in")


# In[ ]:


get_ipython().run_cell_magic('time', '', 'snap.tl.macs3(data, groupby=\'celltype_final\',key_added=\'celltype_final\', blacklist="/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed")\ndata.close()\n')


# In[ ]:


data


# In[ ]:


data.close()


# In[ ]:


data=snap.read_dataset("combined.h5ads")


# In[ ]:


peaks = snap.tl.merge_peaks(data.uns['macs3'], chrom_sizes=snap.genome.mm10)


# In[ ]:


peaks.head()


# In[ ]:


get_ipython().run_cell_magic('time', '', "peak_mat = snap.pp.make_peak_matrix(data, use_rep=peaks['Peaks'])\npeak_mat\n")


# In[ ]:


data.close()


# In[ ]:


print("in")


# In[ ]:





# In[ ]:


#meta = pd.read_csv("meta_bestcelltype.csv", sep = ",")
meta.head()


# In[ ]:


data.obs['age'] = meta['age']
data.obs['best_celltype_batch_age'] = data.obs['best_celltype_batch']+":"+data.obs['age']


# In[ ]:


get_ipython().run_cell_magic('time', '', 'peak_mat = snap.pp.make_peak_matrix(data, chunk_size=5000 , \n                                    file = "peak_mat_iter.h5ad", \n                                    peak_file= "filteredNfixed.union.peakSet.bed")\npeak_mat\n')


# In[ ]:


peak.obs['sample_ct'] = peak.obs['sample'] + peak.obs['best_celltype']


# In[ ]:


peaks_agg=snap.tl.aggregate_X(peak, groupby="sample_ct" , file = "sample_ct_peak_agg_count.h5ad")


# In[ ]:


peaks_agg=snap.tl.aggregate_X(peak, groupby="best_celltype_batch_age" , file = "best_celltype_batch_age_peak_agg2.h5ad", normalize = "RPM")


# In[ ]:


peak = snap.read("peak_mat_iter.h5ad")


# In[ ]:


peaks_agg = snap.read(filename='best_celltype_batch_age_peak_agg2.h5ad')


# In[ ]:


p = peaks_agg.to_memory()


# In[ ]:


df = pd.DataFrame(p.X, columns = p.var.index.tolist())


# In[ ]:


df.index =  p.obs.index.tolist()


# In[ ]:


df.to_csv("sample_ct_peak_agg_count.csv")


# In[ ]:


import os
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count

def process_ct(ct):
    file = "diff_peaks_" + ct + "_2vs18_iter.csv"

    # Check if file already exists
    if os.path.exists(file):
        print(f"Skipping {ct}, output file already exists.\n")
        return

    print(ct)
    group1 = ct + ":2mo"
    group2 = ct + ":18mo"
    g1 = data.obs['best_celltype_batch_age'] == group1
    g2 = data.obs['best_celltype_batch_age'] == group2
    if g1.sum() < 25 or g2.sum() < 25: 
        print("not enough cells\n")
        return

    peaks_selected = np.logical_or(
        (df.loc[group1] > 1).to_numpy(),
        (df.loc[group2] > 1).to_numpy(),
    )
    selected_columns = df.columns[peaks_selected]
    selected_columns_df = pd.DataFrame(selected_columns, columns=['Selected Peaks'])
    selected_columns_df.to_csv(ct + "_peaks_selected.csv", index=False)

    diff_peaks = snap.tl.diff_test(
        peak,
        cell_group1=g1,
        cell_group2=g2,
        min_log_fc=0.15, min_pct=0.01,
        features=peaks_selected,
    )
    column_names = diff_peaks.columns
    diff_peaks_df = pd.DataFrame(diff_peaks, columns=column_names)
    diff_peaks_df.head()
    diff_peaks_df.to_csv(file)
    print("DONE")

if __name__ == "__main__":
    uniq_cts = data.obs['best_celltype_batch'].unique()
    # Filter out the cts that already have an output file
    existing_cts = []
    for ct in uniq_cts:
        file = "diff_peaks_" + ct + "_2vs18_iter.csv"
        if os.path.exists(file):
            print(f"Skipping {ct}, output file already exists.")
            existing_cts.append(ct)

    uniq_cts = np.setdiff1d(uniq_cts, existing_cts)
    # Determine the number of processes to use (up to the number of available CPU cores)
    num_processes = 10

    # Use multiprocessing Pool to parallelize the loop
    with Pool(num_processes) as pool:
        pool.map(process_ct, uniq_cts)


# In[ ]:


uniq_cts = data.obs['best_celltype_batch'].unique()
uniq_cts


# In[ ]:


ct = "HPFGL:Female"
file = "diff_peaks_" + ct + "_2vs18_iter_.05.csv"
    
# Check if file already exists
if os.path.exists(file):
    print(f"Skipping {ct}, output file already exists.")
print(ct)
print("\n")
group1 = ct + ":2mo"
group2 = ct + ":18mo"
g1 = data.obs['best_celltype_batch_age'] == group1
g2 = data.obs['best_celltype_batch_age'] == group2
if g1.sum() < 25 or g2.sum() < 25 : 
    print("not enough cells")
peaks_selected = np.logical_or(
    (df.loc[group1]>1).to_numpy(),
    (df.loc[group2]>1).to_numpy(),
)
selected_columns = df.columns[peaks_selected]
selected_columns_df = pd.DataFrame(selected_columns, columns=['Selected Peaks'])
selected_columns_df.to_csv(ct + "_peaks_selected.csv", index=False)
diff_peaks = snap.tl.diff_test(
    peak,
    cell_group1=g1,
    cell_group2=g2,
    min_log_fc = 0.15,min_pct = 0.05,
    features=peaks_selected
)
column_names = diff_peaks.columns
diff_peaks_df = pd.DataFrame(diff_peaks, columns=column_names)
diff_peaks_df.to_csv(file)
print("DONE")


# In[ ]:


df.loc["HPFGL:Female:9mo", "chr15:74551850-74552349"]


# In[ ]:


uniq_cts = data.obs['best_celltype_batch'].unique()
for ct in uniq_cts:
    file = "diff_peaks_" + ct + "_2vs18_iter.csv"
    
    # Check if file already exists
    if os.path.exists(file):
        print(f"Skipping {ct}, output file already exists.")
        continue
    print(ct)
    print("\n")
    group1 = ct + ":2mo"
    group2 = ct + ":18mo"
    g1 = data.obs['best_celltype_batch_age'] == group1
    g2 = data.obs['best_celltype_batch_age'] == group2
    if g1.sum() < 25 or g2.sum() < 25 : 
        print("not enough cells")
        continue
    peaks_selected = np.logical_or(
        (df.loc[group1]>1).to_numpy(),
        (df.loc[group2]>1).to_numpy(),
    )
    selected_columns = df.columns[peaks_selected]
    selected_columns_df = pd.DataFrame(selected_columns, columns=['Selected Peaks'])
    selected_columns_df.to_csv(ct + "_peaks_selected.csv", index=False)
    diff_peaks = snap.tl.diff_test(
        peak,
        cell_group1=g1,
        cell_group2=g2,
        min_log_fc = 0.15,min_pct = 0.01,
        features=peaks_selected
    )
    column_names = diff_peaks.columns
    diff_peaks_df = pd.DataFrame(diff_peaks, columns=column_names)
    diff_peaks_df.to_csv(file)
    print("DONE")


# In[ ]:


data.close()
peak.close()
peaks_agg.close()


# In[ ]:


print("out")


# In[ ]:


peak.close()


# In[ ]:


peak = snap.read("peak_mat_iter.h5ad")


# In[ ]:


data=snap.read_dataset("combined.h5ads")


# In[ ]:


last_characters = [string[-1] for string in data.obs['sample']]


# In[ ]:


data.obs['rep'] = last_characters


# In[ ]:


data.obs['best_celltype_batch_age_rep'] = data.obs['best_celltype_batch_age']+"_"+data.obs['rep']


# In[ ]:


peak.obs['best_celltype_batch_age_rep'] = data.obs['best_celltype_batch_age_rep']


# In[ ]:


ag = snap.tl.aggregate_X(peak, groupby= 'best_celltype_batch_age_rep', normalize = "RPM", file = "best_celltype_batch_age_rep_None.h5ad")


# In[ ]:


data.close()


# In[ ]:


p = ag.to_memory()


# In[ ]:


p


# In[ ]:


df = pd.DataFrame(p.X, index=p.obs_names.tolist(), columns=p.var.index.tolist())

# Define the path to save the CSV file
csv_file_path = "best_celltype_batch_age_rep_None.csv"

# Save the DataFrame to a CSV file
df.to_csv(csv_file_path)


# In[ ]:


del p


# In[ ]:


peak.close()


# In[ ]:


data.close()


# In[ ]:


print("out")


# In[ ]:


gmat=snap.read("gmat.hdf5")


# In[ ]:


data=snap.read_dataset("combined_save3.h5ads")


# In[ ]:


data


# In[ ]:


snap.pp.make_gene_matrix(data, "/projects/ps-renlab/share/Pipelines/rna-seq/annotation/gtf/gencode.vM10.annotation.gtf", inplace=False, file="gmat.hdf5", backend='hdf5', chunk_size=2500, use_x=False, id_type='gene')


# In[ ]:


data


# In[ ]:


snap.genome.mm10


# In[ ]:


gmat.var


# In[ ]:


gmat


# In[ ]:


gmat.obs['sample']


# In[ ]:


meta


# In[ ]:


gmat.obs['best_celltype'] = meta['best_celltype']


# In[ ]:


data.close()


# In[ ]:


gmat.X


# In[ ]:


g = gmat.to_memory()


# In[ ]:


g.X


# In[ ]:


g.var.index.tolist()


# In[ ]:


g.X


# In[ ]:


df = pd.DataFrame(g.X.toarray(), index=g.obs_names.tolist(), columns=g.var.index.tolist())


# In[ ]:


df.to_csv("cell_by_gene.csv",compression="gzip")


# In[ ]:


g.X.shape


# In[ ]:


g.X.toarray()


# In[ ]:


df.head()


# In[ ]:


gmat.close()


# In[ ]:




