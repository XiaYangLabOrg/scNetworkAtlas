#!/bin/bash

# Run from root directory /path/to/Allen_10X/
##### INPUTS #####
base_dir=$1
mapping_file=$2
adata_dir=$3
celltype_column=$4
##### INPUTS #####

tissue_dir="tissue_adata/"
mkdir $tissue_dir

qsub ../shell_scripts/run_cellmapping.sh ${mapping_file} ${adata_dir} ${celltype_column} ${tissue_dir}
