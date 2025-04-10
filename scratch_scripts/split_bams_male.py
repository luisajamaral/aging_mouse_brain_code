#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from multiprocessing import Pool
import os
import pandas as pd
import numpy as np
import pysam

#os.chdir('/home/lamaral/ps-renlab2/projects/combined_all')

#os.chdir('/tscc/projects/ps-renlab2/lamaral/projects/combined_all')


meta = pd.read_csv("final_meta.csv", sep = ",")
#male_bams = pd.read_csv("/home/lamaral/ps-renlab2/projects/combined_all/male_files.csv", sep = ",")

male_bams = pd.read_csv("male_files.csv", sep = ",")

threads = 14


# In[ ]:


male_bams["Sample"] = "Male:" + male_bams["Sample_id"]


# In[ ]:


all(sample in pd.unique(meta["sample"]) for sample in male_bams["Sample"])


# In[ ]:


male_bams
male_bams['path'] = male_bams['path'].str.replace('/projects/ps-renlab/lamaral/','/tscc/projects/ps-renlab/lamaral/')


# In[ ]:


male_bams


# In[ ]:


meta['barcode'] = meta.apply(lambda row: row['cell_id'].split(':')[2], axis=1)


# In[ ]:


meta.head()


# In[ ]:


np.unique(meta['celltype_final'].astype(str))[0:]


# In[ ]:


#os.chdir('/home/lamaral/scratch/combined_all/')
#os.chdir('/home/lamaral/scratch/combined_all/')
#os.chdir('../')


# In[ ]:


os.chdir('../')


# In[ ]:


os.chdir('combined_all/male_split_bams/')


# In[ ]:


np.unique(meta['celltype_final'].astype(str))[0:-1]


# In[ ]:


meta.head()


# In[ ]:


meta_filtered = meta[meta['celltype_final'] != 'doublet']
np.unique(meta_filtered['celltype_final'].astype(str))[0:]


# In[ ]:


meta = meta_filtered


# In[ ]:


np.unique(meta['celltype_final'].astype(str))[::-1]


# In[ ]:


meta['celltype_age'] = meta['celltype_final']+meta['age']


# In[ ]:


pd.set_option('display.max_rows', len(meta['celltype_age'].value_counts()))

print(meta['celltype_final'].value_counts())


# In[ ]:


meta['celltype_final'] = meta['celltype_final'].astype(str).str.replace(" ", "_").str.replace("/", "_")


# In[ ]:


meta_filtered = meta[meta['celltype_final'] == 'L2_3_IT_CTX_Glut']
np.unique(meta_filtered['celltype_final'].astype(str))[0:]


# In[ ]:


np.unique(meta['celltype_final'].astype(str))[0:]


# In[ ]:


import os
import numpy as np
import pandas as pd
import pysam
from multiprocessing import Pool, cpu_count

NCPU = 14

def run():
    """Entry point of the program"""
    outPrefix = "./male_split_bams/"
    print("Filtering out bam files")
    generate_bams(male_bams, meta, outPrefix)

def generate_bams(male_bams, meta, outPrefix):
    """Generate separate bam files based on cell types and samples"""
    p = Pool(NCPU)
    for idx in np.unique(meta['celltype_final'].astype(str))[::-1]:
        for sample in np.unique(male_bams['Sample'].astype(str)):
            bamf = male_bams.loc[male_bams["Sample"] == sample, "path"].values[0]
            idx = idx.replace(" ", "_").replace("/", "_")
            output_bam = f"{outPrefix}{idx}.{sample}.bam"#.replace(" ", "_").replace("/", "_")
            
            # Check if the output BAM file already exists
            if os.path.exists(output_bam):
                print(f"Skipping generation for metaCell {idx} and sample {sample}, output BAM file already exists.")
                continue
            #print(output_bam)
            output_bam = f"{outPrefix}metacell_{idx}.{sample}.bam"#.replace(" ", "_").replace("/", "_")

            p.apply_async(generate_bam_worker, (bamf, meta, idx, sample, outPrefix))
    p.close()
    p.join()

def generate_bam_worker(bamf, meta, cluster, sample, prefix):
    """Worker function to generate filtered bam files"""
    print(cluster)
    name = bamf
    print(name)
    print(str.encode(sample))
    bamF = pysam.AlignmentFile(name)
    qnames = list(meta[(meta['celltype_final'].astype(str).replace(" ", "_").replace("/", "_") == cluster) & (meta['sample'] == sample)]['barcode'].astype(str))
    qnames_set = set(qnames)
    print(bamF)
    print(qnames_set)
    if len(qnames_set) > 0:
        cluster = cluster.replace(" ", "_").replace("/", "_")
        bam_fname = f"{prefix}metacell_{cluster}.{sample}.bam"#.replace(" ", "_").replace("/", "_")
        print("For metaCell =", cluster, "The filtered bam is writing to:", bam_fname)
        obam = pysam.AlignmentFile(bam_fname, "wb", template=bamF)
        for b in bamF.fetch(until_eof=True):
            if b.get_tag('BX') in qnames_set:
                obam.write(b)
        obam.close()
        bamF.close()
        print("metaCell =", cluster, "Sample =", sample, "writing finished.")

if __name__ == "__main__":
    run()


# In[ ]:


meta['celltype_final'].astype(str).replace(" ", "_").replace("/", "_")[1]


# In[ ]:


cluster
meta['celltype_final'].astype(str)[1]


# In[ ]:


(meta['celltype_final'].astype(str) == cluster)


# In[ ]:


meta['celltype_final'].astype(str).value_counts()


# In[ ]:




