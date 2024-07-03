#!/bin/bash

# Run from root directory /path/to/Allen_10X/
tissue_dir=$1
supercell_dir=$2
celltype_col=$3
sample_col=$4
tissue_celltype_file=$5


if [[ ! -f ${tissue_dir}/${tissue_celltype_file} ]];
then
    find $tissue_dir | grep h5ad | awk -F'.h5ad' '{print $1}' | sort > ${tissue_dir}/${tissue_celltype_file}
fi

num_lines=$(cat ${tissue_dir}/${tissue_celltype_file} | wc -l)
echo $num_lines


qsub -t 1-${num_lines}:1 temp/shell_scripts/run_supercells_v2.sh ${tissue_dir}/${tissue_celltype_file} ${supercell_dir} ${celltype_col} ${sample_col}
