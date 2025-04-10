#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from multiprocessing import Pool
import os
import pandas as pd
import numpy as np
import pysam

#os.chdir('/home/lamaral/projects/combined/')
meta = pd.read_csv("meta_bestcelltype.csv", sep = ",")
#female_bams = pd.read_csv("/projects/ps-renlab2/lamaral/projects/fmouse_multiome/tracksheet_with_paths3.csv", sep = ",")
female_bams =  pd.read_csv("tracksheet_with_paths_oct1723.csv", sep = ",")
threads = 8


# In[ ]:


meta = pd.read_csv("female_RNA/meta_final.csv", sep = ",")


# In[ ]:


meta.head()


# In[ ]:


female_bams.head()


# In[ ]:


female_bams["Sample"].head()


# In[ ]:


female_bams["Sample"] = "Female:" + female_bams["Sample_name"]


# In[ ]:


meta["Sample"] = "Female:" + meta["orig.ident"]


# In[ ]:


all(sample in pd.unique(meta["Sample"]) for sample in female_bams["Sample"])


# In[ ]:


meta['barcode'] = meta.apply(lambda row: row['barcode'].split(':')[1], axis=1)
#meta = meta[meta["batch"] == "Female"]


# In[ ]:


meta['barcode'] = meta['barcode']+"-1"
meta.head()


# In[ ]:


meta['celltype_final'].unique()
meta['celltype_final'] = meta['celltype_final'].str.replace('/', '-')
meta['celltype_final'] = meta['celltype_final'].str.replace(' ', '_')

meta['celltype_final'].unique()


# In[ ]:


female_bams.head()


# In[ ]:


NCPU = 16

def run():
    """Entry point of the program"""
    outPrefix = "./female_RNA/female_split_bams/"

    print("Filtering out bam files")
    generate_bams(female_bams, meta, outPrefix)

def generate_bams(female_bams, meta, outPrefix):
    """Generate separate bam files based on cell types and samples"""
    p = Pool(NCPU)
    for idx in np.unique(meta['celltype_final'].astype(str)):
        for sample in np.unique(female_bams['Sample'].astype(str)):
            bamf = female_bams.loc[female_bams["Sample"] == sample, "Path-luisa"].values[0]
            bamf = "/" + bamf + "/outs/gex_possorted_bam.bam"
            p.apply_async(generate_bam_worker, (bamf, meta, idx, sample, outPrefix))
    p.close()
    p.join()

def generate_bam_worker(bamf, meta, cluster, sample, prefix):
    """Worker function to generate filtered bam files"""
    print(cluster)
    print(sample)
    # Check if the output BAM file already exists
    bam_fname = f"{prefix}{cluster}.{sample}.bam"
    if os.path.exists(bam_fname):
        print("For metaCell =", cluster, "Sample =", sample, "BAM file already exists. Skipping.")
        bamF.close()
        return
    name = bamf
    bamF = pysam.AlignmentFile(name)
    qnames =  list(meta[(meta['celltype_final'].astype(str) == cluster) & (meta['Sample'] == sample)]['barcode'].astype(str))
    qnames_set = set(qnames)
    print(len(qnames_set))
    if len(qnames_set) > 0:
        bam_fname = f"{prefix}{cluster}.{sample}.bam"
        print("For metaCell =", cluster, "The filtered bam is writing to:", bam_fname)
        obam = pysam.AlignmentFile(bam_fname, "wb", template=bamF)
        for b in bamF.fetch(until_eof=True):
            if not b.has_tag("CB"):
                continue
            if b.get_tag('CB') in qnames_set:
                obam.write(b)
        obam.close()
        bamF.close()
        print("metaCell =", cluster, "Sample =", sample, "writing finished.")


# In[ ]:


if __name__ == "__main__":
    run()


# In[ ]:


print("out")


# In[ ]:




