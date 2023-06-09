#!/bin/bash

# Run from root directory /path/to/Allen_10X/
##### INPUTS #####
base_dir="/u/project/xyang123/shared/reference/single_cell_databases/"
mapping_file="${base_dir}all_celltypes/human_Allen_10X.cleaned.tsv"
adata_dir="${base_dir}human/Allen_10X/adatas/"
celltype_column="celltypes"
##### INPUTS #####

tissue_dir="tissue_adata/"
mkdir $tissue_dir

qsub shell_scripts/run_cellmapping.sh ${mapping_file} ${adata_dir} ${celltype_column} ${tissue_dir}
