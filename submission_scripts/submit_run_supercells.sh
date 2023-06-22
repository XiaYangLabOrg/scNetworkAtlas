#!/bin/bash

# Run from root directory /path/to/Allen_10X/

tissue_dir="tissue_adata/"
supercell_dir="supercells"

cd ${tissue_dir}
tissue_celltype_file=tissue_celltype_file.txt

find ./ | grep genes | awk -F'_genes.txt' '{print $1}' > ${tissue_celltype_file}
cd ../


num_lines=$(cat ${tissue_dir}/${tissue_celltype_file} | wc -l)

echo $num_lines
qsub -t 1-${num_lines}:1 shell_scripts/run_supercells.sh ${tissue_dir} ${tissue_celltype_file} ${supercell_dir}
