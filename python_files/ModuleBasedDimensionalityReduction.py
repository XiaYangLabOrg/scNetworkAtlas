#!/usr/bin/env python
# coding: utf-8

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



def remove_small_modules(gene_membership, min_module_size):
    min_module_size
    group_count = gene_membership.cluster_membership.value_counts()
    groups_to_keep = group_count[group_count > min_module_size].index
    return gene_membership[gene_membership.cluster_membership.isin(groups_to_keep)]

def binary_search(l_idx, r_idx, target_Q1_size, precision, verbose=False):
    verboseprint = print if verbose else lambda *a, **k: None # function prints only if verbose=True

    while r_idx - l_idx >= 0.0001:
        resolution_parameter = (l_idx + r_idx) / 2
        
        partition = leidenalg.find_partition(G, leidenalg.CPMVertexPartition, resolution_parameter = resolution_parameter,seed=0)
        groups = np.array(partition.membership)
        
        # organizing gene membership dataframe
        genes = G.vs['name']
        gene_membership = pd.DataFrame(np.array([genes,groups]).T, columns = ['genes','cluster_membership'])
        gene_membership = gene_membership.sort_values('genes')

        # remove small modules with 10 or fewer genes
        verboseprint('Filtering small modules out')
        gene_membership = remove_small_modules(gene_membership, min_module_size=10)

        if gene_membership.shape[0]==0:
            verboseprint('All modules filtered out')
            # decrease resolution
            r_idx = resolution_parameter
            continue
        
        # record the size of each module
        verboseprint('Checking 25th percentile module size')
        module_sizes = gene_membership['cluster_membership'].value_counts().reset_index(name='size')
        module_sizes.columns = ['Module','Module_Size']
        
        q1 = module_sizes['Module_Size'].quantile(0.25)
        
        verboseprint(l_idx, resolution_parameter, r_idx)

        if abs(q1 - target_Q1_size) <= precision:
            print('resolution acquired')
            break

        elif q1 < target_Q1_size - precision: ## Q1 is too small
            r_idx = resolution_parameter
        
        else: ## Q1 is too large
            l_idx = resolution_parameter
        
        

    if r_idx - l_idx >= 0.0001:
        return (resolution_parameter, gene_membership)
    else:
        return resolution_parameter, None



if __name__ == '__main__':
    network = pd.read_csv(sys.argv[1], index_col=0)
    outfile = sys.argv[2]
    q1 = int(sys.argv[3])
    # Q1s = sys.argv[3]

    # network = pd.read_csv('../saved_networks/final_edges/Bladder_bladder_cell.csv.gz',index_col=0)
    # outfile = './gene_memberships/Bladder_bladder_cell.gene_membership.20.csv.gz'
    # q1 = 20
    # create and partition graph
    G = ig.Graph.TupleList([tuple(x) for x in network.values], directed = True)

    precision = 2
    max_res=1
    min_res=0.00001

    while True:
        resolution_parameter, gene_membership = binary_search(min_res, max_res, q1, precision, verbose=False)
        if type(gene_membership) == pd.DataFrame:
            print(q1,resolution_parameter)
            print('module size stats: ')
            print(gene_membership.cluster_membership.value_counts().describe())
            gene_membership.to_csv(outfile, index=False)
            break
        else:
            print("Binary Search failed at " + str(resolution_parameter) + ". Increase search space")
            max_res+=1
