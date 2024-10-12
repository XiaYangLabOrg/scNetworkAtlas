#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import networkx as nx
import leidenalg
import igraph as ig
import argparse




### module-> module

def remove_small_modules(gene_membership, min_module_size):
    # min_module_size
    group_count = gene_membership.module.value_counts()
    groups_to_keep = group_count[group_count > min_module_size].index
    return gene_membership[gene_membership.module.isin(groups_to_keep)]

def binary_search(l_idx, r_idx, target_Q1_size, precision, verbose=False):
    verboseprint = print if verbose else lambda *a, **k: None # function prints only if verbose=True

    while r_idx - l_idx >= 0.0001:
        resolution_parameter = (l_idx + r_idx) / 2
        
        partition = leidenalg.find_partition(G, leidenalg.CPMVertexPartition, resolution_parameter = resolution_parameter,seed=0)
        groups = np.array(partition.membership)
        
        # organizing gene membership dataframe
        genes = G.vs['name']
        gene_membership = pd.DataFrame(np.array([genes,groups]).T, columns = ['genes','module'])
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
        module_sizes = gene_membership['module'].value_counts().reset_index(name='size')
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

def compute_modularity(G_nx, r):
    community = nx.community.louvain_communities(G_nx, resolution=r)
    return nx.community.modularity(G_nx, community, resolution=r), community

def GeneMemByMod_optimizer(G, min_module_size, max_module_size, partition_type=leidenalg.CPMVertexPartition, resolution_range=(0,0.01), return_stats=False):
    optimiser = leidenalg.Optimiser()
    optimiser.set_rng_seed(0)
    # run partitioning across range of resolutions defined by optimizer
    profile = optimiser.resolution_profile(G, partition_type=partition_type,
                                        resolution_range=resolution_range)
    # summarize results
    if return_stats:
        partition_stats = pd.DataFrame(columns=['res','modularity', 'total_modules', 'total_genes','remaining_modules','remaining_genes', 'max_module_size'])
    
    # identify resolution with highest modularity and with modules smaller than max module size.
    best_partition = profile[0]
    for p in profile:
        res = p.resolution_parameter
        max_size = max(p.sizes())
        modularity = p.modularity
        # add to summary results
        if return_stats:
            total_modules = len(p.sizes())
            total_genes = p.n
            mod_sizes = np.array(p.sizes())
            remaining_modules = mod_sizes[((mod_sizes >= min_module_size) & (mod_sizes<=max_module_size))]
            remaining_modules_num = remaining_modules.shape[0]
            remaining_genes = remaining_modules.sum()
            partition_stats.loc[partition_stats.shape[0]] = [res, modularity, total_modules, total_genes, remaining_modules_num, remaining_genes, max_size]
        # update best partition
        if max_size <= max_module_size and modularity > best_partition.modularity:
            best_partition = p
    
    # get gene membership of best partition
    groups = np.array(best_partition.membership)
    genes = G.vs['name']
    gene_membership = pd.DataFrame(np.array([genes,groups]).T, columns = ['genes','module'])
    gene_membership = remove_small_modules(gene_membership, min_module_size=min_module_size)
    print(f"best resolution: {best_partition.resolution_parameter} \t| best modularity: {best_partition.modularity}")
    if return_stats:
        return (gene_membership, partition_stats, best_partition.resolution_parameter)
    else:
        return gene_membership

# function for gene membership
def GeneMemByMod_cpm(G, threshold = 10, maxsize_threshold=300, resolutions=list(np.arange(0.001,0.02,0.001))):
    best, count = 0, 0
    best_r = 0
    # mod_li = []
    # n_li = []
    best_community = None
    # optimization
    for r in resolutions: # add range choice - or 
        _community = leidenalg.find_partition(G, leidenalg.CPMVertexPartition, resolution_parameter = r,seed=0)
        modularity = _community.modularity # modify later, resolution parameters seems not a good parameter to optimize
        # mod_li.append(modularity)
        # n_li.append(len(_community))
        # gene_membership, module = [],[]
        # count = 0
        # for community in _community:
        #     if len(community) < threshold:
        #         continue
        #     else:
        #         for gene in community:
        #             gene_membership.append(gene)
        #             module.append(count)
        #         count += 1
        # gene_membership = pd.DataFrame([gene_membership, module]).T
        # gene_membership.columns = ['genes','module']
        
        # module_size = gene_membership.module.value_counts().reset_index()
        # module_size.columns = ['module','size']
        if max(_community.sizes()) > maxsize_threshold:
            continue
        if best < modularity:
            best_r = r
            best = modularity
            best_community = _community
            # count = 0
    print(f"best resolution: {best_r} \t| best modularity: {best}")
    
    groups = np.array(best_community.membership)
    # organizing gene membership dataframe
    genes = G.vs['name']
    gene_membership = pd.DataFrame(np.array([genes,groups]).T, columns = ['genes','module'])

    # remove small modules with 10 or fewer genes
    gene_membership = remove_small_modules(gene_membership, min_module_size=threshold)
    return gene_membership


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Module Detection")
    parser.add_argument('input', type=str, help='path to network file')
    parser.add_argument('output', type=str, help='output filename')
    parser.add_argument('min_module_size', type=int, default = 10, help='threshold fot minimum module size')
    parser.add_argument('max_module_size', type=int, default = 300, help='threshold fot maximum module size')
    
    args = parser.parse_args()
    input = args.input
    output = args.output
    min_module_size = args.min_module_size
    max_module_size = args.max_module_size
    
    if input.endswith('txt'):
        network = pd.read_csv(input, sep='\t')
    elif input.endswith('csv') or input.endswith('csv.gz'):
        network = pd.read_csv(input)
    else:
        print('File must be txt or csv')
        exit()
    #G_nx = nx.DiGraph(G.get_edgelist())
    #G_nx = nx.DiGraph(network)
    G = ig.Graph.TupleList([tuple(x) for x in network.values], directed = True)
    
    # gene_membership = GeneMemByMod_cpm(G, threshold=min_module_size, maxsize_threshold=max_module_size)
    gene_membership = GeneMemByMod_optimizer(G, min_module_size=min_module_size, max_module_size=max_module_size, resolution_range=(0,0.01))
    
    print('module size stats: ')
    print(gene_membership.module.value_counts().describe())
    
    if output.endswith('.csv'):
        gene_membership.to_csv(output, index=False)
    else:
        gene_membership.to_csv(output, index=False, sep='\t')
    