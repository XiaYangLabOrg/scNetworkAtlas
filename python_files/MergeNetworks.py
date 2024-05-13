#!/usr/bin/env python
# coding: utf-8

if __name__ == "__main__":

	import os
	import warnings
	warnings.filterwarnings("ignore")
	import sys
	import pandas as pd
	import scanpy as sc
	from scing import merge
	import psutil


	adata_merged = sc.read(sys.argv[1])
	all_edge_files = sys.argv[2]
	consensus = float(sys.argv[3])
	outfile = sys.argv[4]
	remove_cycles = False
	print(f"number of cpus: {psutil.cpu_count(logical=False)}")
	outdir = '/'.join(outfile.split('/')[:-1])
	os.makedirs(outdir, exist_ok=True)
	# turn all edge files into one mega file
	files = os.listdir(all_edge_files)
	files = [f for f in files if ".csv.gz" in f]

	all_edges = []
	for i in files:
		all_edges.append(pd.read_csv(all_edge_files + i))

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
