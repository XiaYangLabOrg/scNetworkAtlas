#!/bin/bash

# Run from root directory /path/to/Allen_10X/
##### INPUTS #####
num_networks=100
##################

supercell_dir="supercells"
supercell_file="supercells.txt"
cd $supercell_dir

if [ ! -f ${supercell_file} ]
then
    find ./ | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi
cd ../

num_lines=$(cat ${supercell_dir}/${supercell_file} | wc -l)
echo $num_lines

mkdir saved_networks
mkdir saved_networks/intermediate_data/
while read celltype;
do
        qsub -t 1:${num_networks} shell_scripts/run_buildgrn.sh ${celltype} ${supercell_dir}

done < ${supercell_dir}/$supercell_file
