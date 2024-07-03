#!/usr/bin/env python
# coding: utf-8

# Set number of threads to use
import os
import scanpy as sc
import sys
from scing import build
import warnings
warnings.filterwarnings("ignore")
import argparse
import psutil

parser = argparse.ArgumentParser(description="Build Intermediate GRNs")
parser.add_argument('adata_file', type=str, help='path to pseudobulk h5ad file')
parser.add_argument('outfile', type=str, help='output filename')
parser.add_argument('-s','--seed', type=int, help='random seed', default=0)

args = parser.parse_args()
adata_file = args.adata_file
outfile = args.outfile
seed = args.seed

if mem_per_core == 0:
    mem_per_core = "auto"
print(f"adata_file: {adata_file}\n\
outfile: {outfile}\n\
seed: {seed}"
)

adata = sc.read_h5ad(adata_file)

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
                         ncore=psutil.cpu_count(),
                         mem_per_core='auto',
                         verbose=True,
                         random_state=seed)
scing.subsample_cells()
scing.filter_genes()
scing.filter_gene_connectivities()
scing.build_grn()

scing.edges.to_csv(outfile, index = False)





