#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import os
import sys

#make dictionary to match broader cell types
if sys.argv[1].endswith('tsv'):
    delim = '\t'
else:
    delim=','

dataframe = pd.read_csv(sys.argv[1], delimiter=delim) 
dictionary = dict(dataframe.values)
print(dataframe)
#get all file names from the folder with cell atlas data
files = os.listdir(sys.argv[2])
files = [f for f in files if ".h5ad" in f] #makes sure all of them are h5ad


# In[13]:


#function to get gene expression matrix as an npz file and get the gene names
def GENEEXPdata(tissuefile):
    celltype_col = sys.argv[3]
    #makes adata object from h5ad file and extracts tissue name
    adata = sc.read(f"{sys.argv[2]}/{tissuefile}")
    if len(tissuefile.split("_"))>1:
        tissueNAME = "_".join(tissuefile.split("_")[1:]).split(".h5ad")[0]
    else:
        tissueNAME = tissuefile.split(".h5ad")[0]
    
    #makes folder with tissue name
    os.makedirs(f"{tissueNAME}", exist_ok=True)
    
    #get the broader cell types
    if not all(adata.obs[celltype_col].isin(dictionary.keys())):
        for i in range(adata.shape[0]):
            OLD = adata.obs[celltype_col][i]
            NEW = "-".join((adata.obs[celltype_col][i]).split(" ")) #hyphenates the original cell types so it matches the dictionary
            adata.obs[celltype_col].replace(to_replace = OLD, value = NEW, inplace=True) #edits the adata with the hyphenated version

    adata.obs["broader_cell_type"] = adata.obs[celltype_col].map(lambda x:dictionary[x]) #writes a new column that has the broader cell type corresponding with the cell

    #extract all the unique broader cell types
    UNIQUEcelltypes = set()
    for i in range(adata.shape[0]):
        UNIQUEcelltypes.add(adata.obs["broader_cell_type"][i]) #adds to set if the broader cell type of that cell is a unique one

    adata = adata.T
    
    #save to npz file
    for i in UNIQUEcelltypes:
        NEWadata = adata[:,adata.var.broader_cell_type == i] #find the subset/portion of the adata that has cells that are that specific broader cell type
        if NEWadata.shape[1] >= 100:                         #makes sure that there are more than 100 cells of that broader cell type for PCA in SCING
            
            #eliminate genes expressed in less than 1% of cells, subset adata to genes expressed in more than 1% of cells
            print(NEWadata.shape)
            NEWadata = NEWadata[np.mean(NEWadata.X.todense() > 0, axis=1) > 0.01,:]
            fileNAME = (tissueNAME)+"_"+("_".join(i.split(" ")))
            
            #extract gene names into a txt file within the folder
            if 'gene_symbol' not in NEWadata.obs.columns:
                NEWadata.obs['gene_symbol'] = NEWadata.obs.index
            with open(f'./{tissueNAME}/{fileNAME}_genes.txt','w') as f:
                for g in NEWadata.obs.gene_symbol:
                    f.write(g+"\n")
            
            SPARSEmatrix = scipy.sparse.csr_matrix(np.asarray(NEWadata.X.todense().astype(np.float32))) #turns that portion into npz file we can work with
            scipy.sparse.save_npz(f'./{tissueNAME}/{fileNAME}', SPARSEmatrix) #saves the npz file
        
            #prevent memory leaks
            del SPARSEmatrix
        
        # del adata
        del NEWadata

    return
# In[14]:


#make the npz files
for i in files:
    GENEEXPdata(i)


# In[ ]:





# In[ ]:




