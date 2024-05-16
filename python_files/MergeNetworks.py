#!/usr/bin/env python
# coding: utf-8

import os
import warnings
warnings.filterwarnings("ignore")
import sys
import pandas as pd
import scanpy as sc
from scing import merge
import psutil
import argparse

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Merge Intermediate GRNs")
	parser.add_argument('adata_file', type=str, help='path to pseudobulk h5ad file')
	parser.add_argument('all_edge_files', type=str, help='directory to intermediate networks')
	parser.add_argument('outfile', type=str, help='output filename')
	parser.add_argument('-c','--consensus', type=float, help='consensus threshold', default=0.5)

	args = parser.parse_args()
	adata_file = args.adata_file
	all_edge_files = args.all_edge_files
	consensus = args.consensus
	outfile = args.outfile

	print(f"adata_file: {adata_file}\n\
	all_edge_files: {all_edge_files}\n\
	consensus: {consensus}\n\
	outfile: {outfile}")

	adata_merged = sc.read_h5ad(adata_file)
	print(f"number of cpus: {psutil.cpu_count(logical=False)}")
	outdir = '/'.join(outfile.split('/')[:-1])
	os.makedirs(outdir, exist_ok=True)
	# turn all edge files into one mega file
	files = os.listdir(all_edge_files)
	files = [f for f in files if ".csv.gz" in f]

	all_edges = []
	for i in files:
		all_edges.append(pd.read_csv(all_edge_files + i))

	remove_cycles = False
	merger = merge.NetworkMerger(adata_merged,
					all_edges,
					consensus,
					remove_cycles,
					outdir,
					outfile.split("/")[-1],
					ncore=psutil.cpu_count(),
					mem_per_core="auto",
					verbose=True)

	merger.preprocess_network_files()
	merger.remove_reversed_edges()
	merger.remove_cycles()
	merger.get_triads()
	merger.remove_redundant_edges()

	all_edges = merger.edge_df.sort_values(by='importance',
							ascending=False)

	all_edges.to_csv(outfile, index = False)
