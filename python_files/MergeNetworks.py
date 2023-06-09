#!/usr/bin/env python
# coding: utf-8

# In[1]:
if __name__ == "__main__":

# Set number of threads to use
	import os
	import warnings
	warnings.filterwarnings("ignore")
	import sys
	sys.path.insert(1,'../../src/')
	from supercellHelpers import *
	from buildGRNHelpers import *
	from MergeNetworksHelpers import *
	import csv

	import scanpy as sc
# In[ ]:


	adata_merged = sc.read(sys.argv[1])
	all_edge_files = sys.argv[2]

#turn all edge files into one mega file
	files = os.listdir(all_edge_files)
	files = [f for f in files if ".csv.gz" in f]

	all_edges = []
	for i in files:
		all_edges.append(pd.read_csv(all_edge_files + i))


# In[3]:


	outfile = sys.argv[1].split("/")[-1].split(".h5ad")[0]
#outfile


# In[13]:


	merger = NetworkMerger(adata_merged,
                    all_edges,
                       0.2,
                    'saved_networks/final_edges',
                    outfile,
                    12,
                    int(2e9),
                    True)


# In[14]:


	merger.preprocess_network_files()
	merger.remove_reversed_edges()
	merger.remove_cycles()
	merger.get_triads()
	merger.remove_redundant_edges()


# In[15]:


	all_edges = merger.edge_df.sort_values(by='importance',
                          ascending=False)


# In[16]:


	all_edges
	all_edges.to_csv("saved_networks/final_edges/"+outfile+'.csv.gz')
	#print(all_edges.source)
	#print(all_edges.target)
	#print(all_edges.importance)

# In[ ]:





# In[ ]:




