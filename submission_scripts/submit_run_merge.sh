#!/bin/bash

# Run from root directory /path/to/Allen_10X/

# make final network directory
mkdir saved_networks/final_edges
supercell_dir="supercells"
supercell_file="supercells.txt"
cd ${supercell_dir}
if [ ! -f ${supercell_file} ]
then
    find ./ | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi
cd ../

num_lines=$(cat ${supercell_dir}/${supercell_file} | wc -l)
echo $num_lines
qsub -t 1-${num_lines}:1 shell_scripts/run_merge.sh $supercell_dir $supercell_file

