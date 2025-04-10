#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Parameters
cpu = 8
group_name = "DG_Glut"
mem_gb = 10


# In[ ]:


import os
import re
import glob
import pysam
import networkx as nx
import pandas as pd
import numpy as np
from itertools import combinations
from subprocess import check_output
import xarray as xr
from pybedtools import BedTool
from collections import defaultdict, Counter
import dask
from ALLCools.plot import *
from ALLCools.mcds import MCDS, RegionDS
from ALLCools.dmr import call_dms, call_dmr
import pathlib
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from ALLCools.mcds.correlation import region_correlation, get_corr_table
from ALLCools.dmr.parse_methylpy import methylpy_to_region_ds
import seaborn as sns
from multiprocessing import Pool, Manager, Process
from functools import partial
import argparse
from pingouin import partial_corr
from ALLCools.mcds import MCDS
import cooler


# In[ ]:


#group_name = 'DG_Glut'


# In[ ]:


ct = group_name


# ## identify deg

# In[ ]:


if ct == 'STR_D1_Gaba':
    deg_ct  = 'STR_D12_Gaba'
    deg_dir = f'/data/female-amb/Diff.Result/DEG.stats/{deg_ct}'
    expr = pd.read_hdf(f'{deg_dir}/expr.hdf').T
    stats = np.load(f"{deg_dir}/{deg_ct}.2mo-{deg_ct}.18mo.npz") #fc:(2mo/18mo)
    result = pd.DataFrame({'fc': (expr[f'{deg_ct}.18mo'] /expr[f'{deg_ct}.2mo']).values,'fdr':stats['fdr']}, index = expr.index)
else:
    deg_dir = f'/data/female-amb/Diff.Result/DEG.stats/{ct}'
    expr = pd.read_hdf(f'{deg_dir}/expr.hdf').T
    stats = np.load(f"{deg_dir}/{ct}.2mo-{ct}.18mo.npz") #fc:(2mo/18mo)
    result = pd.DataFrame({'fc': (expr[f'{ct}.18mo'] /expr[f'{ct}.2mo']).values,'fdr':stats['fdr']}, index = expr.index)


# In[ ]:


result['log2(18mo/2mo)'] = np.log2(result['fc'])
result['-log10(padj)'] = -np.log10(result['fdr'])
result.head()


# In[ ]:


sig_result = result[(result['fdr'] < 0.05) & (abs(result['log2(18mo/2mo)']) > 0.1)]
sig_result.shape


# In[ ]:


fig, ax = plt.subplots(figsize = (3,3), dpi = 150)
sns.scatterplot(data = result,
                ax = ax,
                color = 'lightgrey',
                s = 3,
                x = 'log2(18mo/2mo)',
                y = '-log10(padj)')

sns.scatterplot(data = sig_result[sig_result['log2(18mo/2mo)'] > 0],
                ax = ax,
                color = 'red',
                s = 3,
                x = 'log2(18mo/2mo)',
                y = '-log10(padj)')

sns.scatterplot(data = sig_result[sig_result['log2(18mo/2mo)'] < 0],
                ax = ax,
                color = 'blue',
                s = 3,
                x = 'log2(18mo/2mo)',
                y = '-log10(padj)')

for i, row in sig_result.iterrows():
    if row['-log10(padj)'] > 100 and (row['log2(18mo/2mo)'] > 3 or row['log2(18mo/2mo)'] < -3):
        ax.text(row['log2(18mo/2mo)'], row['-log10(padj)'], 
                row.name, 
                color='black', fontsize=6)


# ## cal abc score

# In[ ]:


chrom_size_path = '/ref/m3C/mm10.main.nochrM.nochrY.chrom.sizes'
gene_meta_path = '/data/metadata/gencode.vM22.basic.annotation.gene.flat.tsv.gz'
dmr_zarr_path = f"{ct}.AllDMR.mcds"
cool_path = '/data/female-amb/AMB.CoolFiles/CellType.Age.Raw.5kb.mcool'
leg = [f"{ct}.{age}" for age in ['8wk','9mo','18mo']]
age_order = ['2mo','9mo','18mo']


# In[ ]:


mouse_chrs = ['chr' + str(x) for x in range(1,20)] + ['chrX']
mouse_size = pd.read_csv(chrom_size_path, sep="\t", index_col=0, names=['length']).loc[mouse_chrs]


# In[ ]:


mouse_genes = pd.read_csv(gene_meta_path,sep="\t")[['chrom', 'start', 'end', 'gene_id','gene_name']]
mouse_genes.columns = ['chrom', 'start', 'end', 'geneID','gene_name']
mouse_genes['geneID'] = [re.sub("\.[0-9]+$", "", x) for x in mouse_genes['geneID']]
mouse_genes.index = mouse_genes['geneID']
mouse_genes = mouse_genes[~mouse_genes['chrom'].isin(['chrY', 'chrM'])]


# In[ ]:


atac_rpm_path = f"/home/qzeng_salk_edu/project/240429_abc_atac/celltype_age_RPM_files/{ct}_RPM.txt"
atac_rpm =pd.read_csv(atac_rpm_path, sep = '\t')


# In[ ]:


peak_bed = pd.DataFrame({'chrom':[_id.split(':')[0] for _id in atac_rpm.index],
                         'start':[int(_id.split(':')[1].split('-')[0]) for _id in atac_rpm.index],
                         'end':[int(_id.split(':')[1].split('-')[-1]) for _id in atac_rpm.index]}, index = atac_rpm.index)
peak_bed.head()


# In[ ]:


atac_rpm =pd.read_csv(atac_rpm_path, sep = '\t')
atac_rpm.columns = [_.split(':')[-1] for _ in atac_rpm.columns] 
dmr_activity_dict = atac_rpm[age_order].T
dmr_activity_dict.index = leg
dmr_activity_dict = dmr_activity_dict.to_dict()


# In[ ]:


def get_gene_abc_score(group, gene_id):
    ct_age_cool = cooler.Cooler(f'{cool_path}/{group}.raw.mcool::resolutions/10000')

    max_distance = 5000000
    ABC_score = defaultdict(dict)
    total_ABC = dict()
    
    gene_coords = mouse_genes.loc[gene_id].to_dict()
    chrom = gene_coords['chrom']
    gene_start = int(gene_coords['start']) - 2000
    gene_end = int(gene_coords['end']) + 2000
    
    start = int(gene_coords['start']) - max_distance
    end = int(gene_coords['start']) + max_distance
    start = 1 if start < 0 else start
    end = mouse_size.loc[chrom, 'length'] if end > mouse_size.loc[chrom, 'length'] else end

    ct_age_dmr = peak_bed
    gene_dmr = ct_age_dmr[(ct_age_dmr['chrom'] == gene_coords['chrom']) & (ct_age_dmr['start'] > start)
                                       & (ct_age_dmr['end'] < end)]
    gene_dmr = gene_dmr[gene_dmr['end'] - gene_dmr['start']>=10]


    # calculate interactions of each DMR to target gene
    dmr_contacts = defaultdict(dict)
    contacts = ct_age_cool.matrix(balance=False, as_pixels=True, join=True).fetch(f'{chrom}:{start}-{end}')

    gene_contacts_upper = contacts[(contacts['start1'] >= gene_start) & (contacts['start1'] <= gene_end)]
    gene_contacts_down = contacts[(contacts['start2'] >= gene_start) & (contacts['start2'] <= gene_end)]
    gene_contacts_upper = gene_contacts_upper[(gene_contacts_upper['start2'] >= gene_start) & (gene_contacts_upper['start2'] <= end)]
    gene_contacts_down = gene_contacts_down[(gene_contacts_down['start1'] >= start) & (gene_contacts_down['start1'] <= gene_start)]

    for z, row in gene_contacts_upper.iterrows():
        z_bin_dmrs = gene_dmr[(gene_dmr['start'] >= row.start2) & (gene_dmr['end'] <= row.end2)]         
        for dmr in z_bin_dmrs.index:
            dmr_contacts[dmr] = row['count']
    
    for y, row in gene_contacts_down.iterrows():
        y_bin_dmrs = gene_dmr[(gene_dmr['start'] >= row.start1) & (gene_dmr['end'] <= row.end1)]         
        for dmr in y_bin_dmrs.index:
            dmr_contacts[dmr] = row['count']
    
    for dmr in dmr_contacts:
        total_ABC = dmr_activity_dict[dmr][group] * dmr_contacts[dmr]

    for dmr in dmr_contacts:
        EG = f'{dmr}-{gene_id}'
        try:
            activity = dmr_activity_dict[dmr][group]
            contact = dmr_contacts[dmr]
            ABC_score[EG] = activity, contact, (activity *  contact/ total_ABC)
        except:
            ABC_score[EG][celltype] = np.nan
    ABC_score_df = pd.DataFrame.from_dict(ABC_score, orient='index')
    #ABC_score_df.columns = [group]
    return ABC_score_df


# In[ ]:


all_genes = mouse_genes[mouse_genes['gene_name'].isin(sig_result.index)].index
len(all_genes)


# In[ ]:


import ray 
ray.init(ignore_reinit_error=True)

@ray.remote(num_cpus = 2)
def get_gene_abc_score(group, gene_id):
    ct_age_cool = cooler.Cooler(f'{cool_path}/{group}.raw.mcool::resolutions/10000')

    max_distance = 5000000
    ABC_score = defaultdict(dict)
    total_ABC = dict()
    
    gene_coords = mouse_genes.loc[gene_id].to_dict()
    chrom = gene_coords['chrom']
    gene_start = int(gene_coords['start']) - 2000
    gene_end = int(gene_coords['end']) + 2000
    
    start = int(gene_coords['start']) - max_distance
    end = int(gene_coords['start']) + max_distance
    start = 1 if start < 0 else start
    end = mouse_size.loc[chrom, 'length'] if end > mouse_size.loc[chrom, 'length'] else end

    ct_age_dmr = peak_bed
    gene_dmr = ct_age_dmr[(ct_age_dmr['chrom'] == gene_coords['chrom']) & (ct_age_dmr['start'] > start)
                                       & (ct_age_dmr['end'] < end)]
    gene_dmr = gene_dmr[gene_dmr['end'] - gene_dmr['start']>=10]


    # calculate interactions of each DMR to target gene
    dmr_contacts = defaultdict(dict)
    contacts = ct_age_cool.matrix(balance=False, as_pixels=True, join=True).fetch(f'{chrom}:{start}-{end}')

    gene_contacts_upper = contacts[(contacts['start1'] >= gene_start) & (contacts['start1'] <= gene_end)]
    gene_contacts_down = contacts[(contacts['start2'] >= gene_start) & (contacts['start2'] <= gene_end)]
    gene_contacts_upper = gene_contacts_upper[(gene_contacts_upper['start2'] >= gene_start) & (gene_contacts_upper['start2'] <= end)]
    gene_contacts_down = gene_contacts_down[(gene_contacts_down['start1'] >= start) & (gene_contacts_down['start1'] <= gene_start)]

    for z, row in gene_contacts_upper.iterrows():
        z_bin_dmrs = gene_dmr[(gene_dmr['start'] >= row.start2) & (gene_dmr['end'] <= row.end2)]         
        for dmr in z_bin_dmrs.index:
            dmr_contacts[dmr] = row['count']
    
    for y, row in gene_contacts_down.iterrows():
        y_bin_dmrs = gene_dmr[(gene_dmr['start'] >= row.start1) & (gene_dmr['end'] <= row.end1)]         
        for dmr in y_bin_dmrs.index:
            dmr_contacts[dmr] = row['count']
    
    for dmr in dmr_contacts:
        total_ABC = dmr_activity_dict[dmr][group] * dmr_contacts[dmr]

    for dmr in dmr_contacts:
        EG = f'{dmr}-{gene_id}'
        try:
            activity = dmr_activity_dict[dmr][group]
            contact = dmr_contacts[dmr]
            ABC_score[EG] = activity, contact, (activity *  contact/ total_ABC)
        except:
            ABC_score[EG][celltype] = np.nan
    ABC_score_df = pd.DataFrame.from_dict(ABC_score, orient='index')
    #ABC_score_df.columns = [group]
    return ABC_score_df


# In[ ]:


for group in [leg[0],leg[1],leg[2]]:
    ct, age = group.split('.')
    results = [get_gene_abc_score.remote(group, gene) for gene in all_genes]
    all_results = ray.get(results)
    all_results = pd.concat(all_results)
    all_results.columns = ['activity', 'contact','abc_score']
    all_results.to_csv(f"{group}.abc_score.csv")
    print(f"{group} done")    


# In[ ]:





# In[ ]:


# import ray 
# ray.init(ignore_reinit_error=True)

# @ray.remote(num_cpus = 2)
# def get_gene_abc_score_by_chunk(group, gene_list):
#     tmp = []
#     for gene_id in gene_list:
#         abc_df = get_gene_abc_score(group, gene_id)
#         tmp.append(abc_df)
#     tmp = pd.concat(tmp)
#     return tmp


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




