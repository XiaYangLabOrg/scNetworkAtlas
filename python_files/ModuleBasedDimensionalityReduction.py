#!/usr/bin/env python
# coding: utf-8

# In[4]:


import leidenalg        # Version 0.8.4
import igraph as ig
import numpy as np
import pandas as pd
from collections import Counter
import scanpy as sc
#import matplotlib.pyplot as plt

#import seaborn as sns
from ctxcore.genesig import GeneSignature 
# UNCOMMENT BELOW AND COMMENT OUT ABOVE IF USING AN OLDER VERSION OF PYSCENIC
#from pyscenic.genesig import GeneSignature
from pyscenic.aucell import aucell, derive_auc_threshold, create_rankings
from sklearn.decomposition import PCA

import sys
sys.path.insert(1,'../../src/')
# In[5]:


adata = sc.read(sys.argv[1])
network = pd.read_csv(sys.argv[2], index_col=0)
resolution_parameter = sys.argv[4]


# In[8]:


#network = network[["source", "target"]]


# In[9]:


#adata


# In[10]:


# create and partition graph
G = ig.Graph.TupleList([tuple(x) for x in network.values],
                       directed = True)

partition = leidenalg.find_partition(G, leidenalg.CPMVertexPartition, resolution_parameter = float(resolution_parameter));
groups = np.array(partition.membership)


# In[11]:


# organizing gene membership dataframe
genes = G.vs['name']
gene_membership = pd.DataFrame(np.array([genes,groups]).T, columns = ['genes','cluster_membership'])
gene_membership = gene_membership.sort_values('genes')


# In[12]:


# remove modules with 2 or fewer genes
group_count = Counter(gene_membership.cluster_membership)
groups_to_keep = []
for g in group_count:
    if group_count[g] > 2:
        groups_to_keep.append(g)
groups_to_keep = np.array(groups_to_keep)
gene_membership = gene_membership.loc[np.isin(gene_membership.cluster_membership,
                                              groups_to_keep)]


# In[13]:


gene_membership.cluster_membership = gene_membership.cluster_membership.astype(str)


# In[14]:


# save gene memberships for pathway enrichment analysis
cell_type = sys.argv[3]
gene_membership.to_csv('./gene_memberships/'+str(sys.argv[3])+'.gene_membership.'+str(sys.argv[4])+'.csv.gz', index=False)



