#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from wmb import brain, cemba, mm10

from ALLCools.clustering import *
from ALLCools.mcds import MCDS
from ALLCools.plot import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy
os.chdir('/home/lamaral/ps-renlab2/projects/combined_all/methylation_data/MCDS/')


# In[ ]:


# Parameters
chrom_to_remove = ["chrX", "chrY", "chrM", "chrL"]
cpu = 1
group_name = "All"
mem_gb = 1
n_cell = 5
remove_lower_features = 0.2
var_dim = "chrom5k"
zscore_abs_cutoff = 3


# In[ ]:


cells = pd.read_csv('mC_META_230814.csv', index_col=0, header=None).index
cells.name = 'cell'


# In[ ]:


cells


# In[ ]:


cemba.CEMBA_SNM3C_3C_CHROM100K_RAW_ZARR_PATH


# In[ ]:


mcds = MCDS.open("/home/lamaral/ps-renlab2/projects/combined_all/methylation_data/MCDS/M.MCDS/mcds_5kb/chunk_9.mcds/*",
                 var_dim=var_dim,
                 use_obs=cells, engine = "scipy"
                )


# In[ ]:





# In[ ]:


cemba.CEMBA_SNMC_MCDS_PATH


# In[ ]:




