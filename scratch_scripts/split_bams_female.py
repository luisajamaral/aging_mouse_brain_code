#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from multiprocessing import Pool
import os
import pandas as pd
import numpy as np
import pysam

os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/')
meta = pd.read_csv("female_RNA/meta_celltypes.csv", sep = ",")
female_bams = pd.read_csv("tracksheet_with_paths_oct1723.csv", sep = ",")
threads = 14


# In[ ]:


female_bams["Sample"] = "Female:" + female_bams["Sample_name"]


# In[ ]:


all(sample in pd.unique(meta["sample"]) for sample in female_bams["Sample"])


# In[ ]:


meta['barcode'] = meta.apply(lambda row: row['cell_id'].split(':')[2], axis=1)


# In[ ]:


#!/bin/python
#  python split_bams.py --tissue all --bam-prefix
#  ../../data/snATAC/bam.filter.sort/ --bam-suffix .filter.csort.bam --statH
#  meta_after.txt --outPrefix ~/scratch/brain_aging_mouse/after_integration/ --cores 12
NCPU = 14

def run():
    """Entry point of the program"""
    outPrefix = "/home/lamaral/scratch/combined_all/female_split_bams/"

    print("Filtering out bam files")
    generate_bams(female_bams, meta, outPrefix)

def generate_bams(female_bams, meta, outPrefix):
    """Generate separate bam files based on cell types and samples"""
    p = Pool(NCPU)
    for idx in np.unique(meta['celltype_final'].astype(str)):
        for sample in np.unique(female_bams['Sample'].astype(str)):
            bamf = female_bams.loc[female_bams["Sample"] == sample, "Path-luisa"].values[0]
            bamf = "/" + bamf + "/outs/atac_possorted_bam.bam"
            p.apply_async(generate_bam_worker, (bamf, meta, idx, sample, outPrefix))
    p.close()
    p.join()

def generate_bam_worker(bamf, meta, cluster, sample, prefix):
    """Worker function to generate filtered bam files"""
    print(cluster)
    print(sample)
    name = bamf
    bamF = pysam.AlignmentFile(name)
    qnames =  list(meta[(meta['celltype_final'].astype(str) == cluster) & (meta['sample'] == sample)]['barcode'].astype(str))
    qnames_set = set(qnames)
    print(len(qnames_set))
    if len(qnames_set) > 0:
        clust = cluster.replace(" ", "_").replace("/", "_")
        bam_fname = f"{prefix}metacell_{clust}.{sample}.bam"
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





# In[ ]:





# In[ ]:




