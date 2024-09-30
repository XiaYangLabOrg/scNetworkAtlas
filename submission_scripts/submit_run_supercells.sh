#!/bin/bash

# Run from root directory /path/to/Allen_10X/

tissue_dir=$1
supercell_dir=$2
filetype=$3
celltype_col=$4
tissue_celltype_file=$5


echo $filetype

if [[ ! -f ${tissue_dir}/${tissue_celltype_file} ]];
then
    find $tissue_dir | grep ${filetype} | awk -F ${filetype} '{print $1}' | sort > ${tissue_dir}/${tissue_celltype_file}
fi

num_lines=$(cat ${tissue_dir}/${tissue_celltype_file} | wc -l)
echo $num_lines

qsub -t 1-${num_lines}:1 temp/shell_scripts/run_supercells.sh ${tissue_dir}/${tissue_celltype_file} ${supercell_dir} ${celltype_col} ${filetype}
