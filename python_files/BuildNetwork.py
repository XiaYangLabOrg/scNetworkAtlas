#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Set number of threads to use
import os
import scanpy as sc
from scing import build
import warnings
warnings.filterwarnings("ignore")


# In[10]:


import sys
#print(f"Name of the script      : {sys.argv[0]=}")
#print(f"Arguments of the script : {sys.argv[1:]=}")


# In[9]:


#args = sys.argv[1:]
#def main(args=None):
    #if args is None:
        #args = sys.argv[1:]
#print(args)

#arg_line = " ".join(sys.argv[1:])
#print(arg_line)


# In[4]:

adata_merged = sys.argv[1]
iteration = sys.argv[2]

outfile = adata_merged.split("/")[-1].split(".h5ad")[0]
adata = sc.read(adata_merged)
# In[ ]:






# In[5]:


#CHANGE RANGE + GRN BUILDER 
all_edges = []

np.random.seed()

adata_saved = adata.copy()

########################################GENES FILTER##############################
# sc.pp.filter_genes(adata_saved, min_cells = (0.1*adata_saved.shape[0]), inplace=True, copy=False)
    
    # -1 all genes
    # 100 neighbors for each gene
    # 10 pcs
    # 0.7 subsample per run
scing = build.grnBuilder(adata_saved, -1, 100, 10,0.7,
                      'test','test',1,int(4e9),True, int(iteration))
scing.subsample_cells()

scing.filter_genes()
scing.filter_gene_connectivities()
scing.build_grn()
  

df_edges = scing.edges

# In[12]:
dirname = f'saved_networks/intermediate_data/{outfile}/'
if (not os.path.exists(dirname)):
	os.makedirs(dirname)

print(dirname+str(outfile)+'.network.'+str(iteration)+'.csv.gz')   
df_edges.to_csv(dirname+str(outfile)+'.network.'+str(iteration)+'.csv.gz', index = False)



# In[ ]:





