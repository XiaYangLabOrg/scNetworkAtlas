#!/usr/bin/env python
# coding: utf-8


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
import argparse

parser = argparse.ArgumentParser(description="Supercells")
parser.add_argument('dataset_name', type=str, help='path to data h5ad file')
parser.add_argument('out_dir', type=str, help='output directory')
parser.add_argument('--stratify_by', nargs='+', type=str, help='a list of column names to stratify by (e.g. cell_type_column sample_column). Separated by whitespace')
parser.add_argument('--save_by', nargs='+', type=str, help='a list of column names to save by (e.g. cell_type_column sample_column). Separated by whitespace')

args = parser.parse_args()
dataset_name = args.dataset_name
out_dir = args.out_dir
stratify_by = args.stratify_by
save_by = args.save_by


# Remove any "None" strings
if stratify_by:
    stratify_by = [i for i in stratify_by if i!='None']
    if len(stratify_by) == 0:
        stratify_by = None 
if save_by:
    # save by columns must also be in stratify
    if stratify_by:
        save_by = list(set(save_by).intersection(stratify_by))
        if len(save_by) == 0:
            save_by = None
    else:
        save_by = None
    print(f"Viable save_by columns that are also in stratify_by: {save_by}")
print(f"dataset_name: {dataset_name}\n\
out_dir: {out_dir}\n\
stratify_by: {stratify_by}\n\
save_by: {save_by}")

os.makedirs(out_dir, exist_ok=True)
pref = dataset_name.split("/")[-1].split(".h5ad")[0]

adata = sc.read_h5ad(dataset_name)
# pseudobulk
pb.pseudobulk_pipeline(adata=adata,
                       stratify_by=stratify_by,
                       save_by=save_by,
                       n_hvgs=2000, n_pcs=40, n_neighbors=10,
                       pb_n_neighbors=10, pb_max_overlap=5, max_pb_pergroup=500, min_pb_pergroup=30, recluster=True,
                       out_dir=out_dir, pref=pref,random_state=0)
