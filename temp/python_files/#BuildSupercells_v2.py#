#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Set number of threads to use
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scing import supercellHelpers as pb
import scipy.sparse
from scipy.sparse import load_npz 
import warnings
warnings.filterwarnings("ignore")

import sys

dataset_name = sys.argv[1]
out_dir = sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
celltype_col = sys.argv[3]
sample_col = sys.argv[4]
if sys.argv[5]=='True':
    recluster=True
else:
    recluster=False
pref = dataset_name.split("/")[-1].split(".h5ad")[0]

adata = sc.read_h5ad(dataset_name)
# pseudobulk
pb.pseudobulk_pipeline(adata=adata,
                       stratify_by=[celltype_col,sample_col],
                       save_by=[celltype_col],
                       n_hvgs=2000, n_pcs=40, n_neighbors=10,
                       pb_n_neighbors=10, pb_max_overlap=5, max_pb_pergroup=500, min_pb_pergroup=30, recluster=recluster,
                       out_dir=out_dir, pref=pref,random_state=0)
