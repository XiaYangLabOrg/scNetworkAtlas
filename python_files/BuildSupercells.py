#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Set number of threads to use
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.sparse import load_npz 
import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.insert(1,'../../../src/')

from supercellHelpers import *
from buildGRNHelpers import *
from MergeNetworksHelpers import *


dataset_name = sys.argv[1]
humanMATRIX = load_npz(dataset_name)

#create outfile as a variable to store the cell type name
outfile = dataset_name.split("/")[-1].split(".npz")[0]

#re-assemble the path to the txt file that contains the names of the cell type's expressed genes from the HCA pre-processing
genes = "/".join(dataset_name.split("/")[:-1]) + "/" + str(outfile)+ "_genes.txt"

adata = sc.AnnData(np.asarray(humanMATRIX.todense())).T

#match the gene names to their corresponding expression data
adata.var.index = pd.read_csv(genes, header = None).to_numpy().ravel()
adata.var_names_make_unique()

#make the supercell
adata_merged = supercell_pipeline(adata,
                                  ngenes=2000,
                                  npcs=20,
                                  ncell=500,
                                  verbose=True)
os.makedirs('../merged_adata',exist_ok=True)
adata_merged.write('../merged_adata/'+str(outfile)+'.h5ad')
