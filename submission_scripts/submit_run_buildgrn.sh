#!/bin/bash

# Run from root directory /path/to/Allen_10X/
##### INPUTS #####
num_networks=100
##################

supercell_dir="supercells"
supercell_file="supercells/supercells.txt"
out_dir="saved_networks/intermediate_data"
mkdir -p $out_dir

if [ ! -f ${supercell_file} ]
then
    ls $supercell_dir | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi

num_lines=$(cat ${supercell_file} | wc -l)
echo $num_lines

while read celltype;
do
    qsub -t 1:${num_networks} shell_scripts/run_buildgrn.sh $celltype $supercell_dir $out_dir

done < $supercell_file
