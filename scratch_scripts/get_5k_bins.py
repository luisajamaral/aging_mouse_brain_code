#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import snapatac2 as snap
import scanpy as scs
import pandas as pd
import numpy as np
import os
import sys

os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/h5ads_filter/anndatas/')


# In[ ]:


"../anndatas_5k/"+filename


# In[ ]:


directory = "." 
adatas = []

for filename in os.listdir(directory):
    if filename.endswith(".h5ad"):
        print(filename)
        data=snap.read(filename)
        snap.pp.add_tile_matrix(data, bin_size=5000, inplace=False, chunk_size=50000, file = "../anndatas_5k/"+filename)
        data
        data.close()

