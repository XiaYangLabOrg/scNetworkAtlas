#!/usr/bin/env python
# coding: utf-8

# Set number of threads to use
import os
import scanpy as sc
import sys
from scing import build
import warnings
warnings.filterwarnings("ignore")

adata_file = sys.argv[1]
outfile = sys.argv[2]
seed = int(sys.argv[3])

adata = sc.read(adata_file)

########################################GENES FILTER##############################
# sc.pp.filter_genes(adata_saved, min_cells = (0.1*adata_saved.shape[0]), inplace=True, copy=False)
outdir = '/'.join(outfile.split('/')[:-1])
os.makedirs(outdir, exist_ok=True)

scing = build.grnBuilder(adata=adata, 
                         ngenes=-1, 
                         nneighbors=100, 
                         npcs=10,
                         subsample_perc=0.7,
                         prefix=outfile.split("/")[-1],
                         outdir=outdir,
                         ncore=1,
                         mem_per_core="auto",
                         verbose=True,
                         random_state=seed)
scing.subsample_cells()
scing.filter_genes()
scing.filter_gene_connectivities()
scing.build_grn()

scing.edges.to_csv(outfile, index = False)





