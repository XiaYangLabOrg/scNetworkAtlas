#!/bin/bash

# Run from root directory /path/to/Allen_10X/

tissue_dir="tissue_adata/"
supercell_dir="supercells"
filetype="npz"


echo $filetype

cd ${tissue_dir}
tissue_celltype_file=tissue_celltype_file.txt

if [[ "${filetype}" == "npz" ]];
then
    find ./ | grep genes | awk -F'_genes.txt' '{print $1}' > ${tissue_celltype_file}
elif [[ "${filetype}" == "h5ad" ]]
then
    find ./ | grep h5ad | awk -F'.h5ad' '{print $1}' > ${tissue_celltype_file}
else
    echo "filetype must be either npz or h5ad"
    exit
fi

cd ../
num_lines=$(cat ${tissue_dir}/${tissue_celltype_file} | wc -l)
echo $num_lines


qsub -t 1-${num_lines}:1 shell_scripts/run_supercells.sh ${tissue_dir} ${tissue_celltype_file} ${supercell_dir} ${filetype}