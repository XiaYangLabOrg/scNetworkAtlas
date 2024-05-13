#!/bin/bash

# Run from root directory /path/to/Allen_10X/
tissue_dir="tissue_adata"
supercell_dir="supercells"
celltype_col=$1
sample_col=$2
recluster=$( [ "$3" = "True" ] && echo true || echo false )
tissue_celltype_file=tissue_celltype_file.txt


if [[ ! -f ${tissue_dir}/${tissue_celltype_file} ]];
then
    find $tissue_dir | grep train.*h5ad | awk -F'.h5ad' '{print $1}' | sort > ${tissue_dir}/${tissue_celltype_file}
fi

num_lines=$(cat ${tissue_dir}/${tissue_celltype_file} | wc -l)
echo $num_lines


qsub -t 1-${num_lines}:1 ../shell_scripts/run_supercells_v2.sh ${tissue_dir}/${tissue_celltype_file} ${supercell_dir} ${celltype_col} ${sample_col} ${recluster}
