#!/usr/bin/env python
# coding: utf-8


# Set number of threads to use
import os
import numpy as np
import pandas as pd
import scipy.sparse
from scipy.sparse import load_npz 
import warnings
warnings.filterwarnings("ignore")

import sys
from scing import supercells as pb
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description="Supercells")
parser.add_argument('dataset_name', type=str, help='path to data h5ad file')
parser.add_argument('out_dir', type=str, help='output directory')
parser.add_argument('--stratify_by', nargs='+', type=str, help='a list of column names to stratify by (e.g. cell_type_column sample_column). Separated by whitespace')

args = parser.parse_args()
dataset_name = args.dataset_name
out_dir = args.out_dir
stratify_by = args.stratify_by
# Remove any "None" strings
if stratify_by:
    stratify_by = [i for i in stratify_by if i!='None']
    if len(stratify_by) == 0:
        stratify_by = None 
print(f"dataset_name: {dataset_name}\n\
out_dir: {out_dir}\n\
stratify_by: {stratify_by}")

if dataset_name.endswith('.npz'):
    humanMATRIX = load_npz(dataset_name)

    #create outfile as a variable to store the cell type name
    outfile = dataset_name.split("/")[-1].split(".npz")[0]

    #re-assemble the path to the txt file that contains the names of the cell type's expressed genes from the HCA pre-processing
    genes = "/".join(dataset_name.split("/")[:-1]) + "/" + str(outfile)+ "_genes.txt"

    adata = sc.AnnData(np.asarray(humanMATRIX.todense())).T

    #match the gene names to their corresponding expression data
    adata.var.index = pd.read_csv(genes, header = None).to_numpy().ravel()
    adata.var_names_make_unique()
else:
    adata = sc.read_h5ad(dataset_name)
    #create outfile as a variable to store the cell type name
    outfile = dataset_name.split("/")[-1].split(".h5ad")[0]

#make the supercell
if stratify_by:
    for colnames,df in adata.obs.groupby(stratify_by):
        adata_ct = adata[df.index]
        adata_merged = pb.supercell_pipeline(adata_ct,
                                        ngenes=2000,
                                        npcs=20,
                                        ncell=500,
                                        verbose=True)
        os.makedirs(out_dir,exist_ok=True)
        if isinstance(colnames, str):
            adata_merged.write(f'{out_dir}/{outfile}_{colnames}.h5ad')
        else:
            adata_merged.write(f'{out_dir}/{outfile}_{"_".join(colnames)}.h5ad')
else:
    adata_merged = pb.supercell_pipeline(adata,
                                    ngenes=2000,
                                    npcs=20,
                                    ncell=500,
                                    verbose=True)
    os.makedirs(out_dir,exist_ok=True)
    adata_merged.write(f'{out_dir}/{outfile}.h5ad')