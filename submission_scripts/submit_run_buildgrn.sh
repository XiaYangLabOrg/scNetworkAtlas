#!/bin/bash

# Run from root directory /path/to/Allen_10X/
##### INPUTS #####
num_networks=$1
supercell_dir=$2
supercell_file=$3
out_dir=$4
ncore=$5
mem_per_core=$6
##################
mkdir -p $out_dir

if [ ! -f ${supercell_file} ]
then
    ls $supercell_dir | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi

num_lines=$(cat ${supercell_file} | wc -l)
echo $num_lines

while read celltype;
do
    qsub -t 1:${num_networks} -l h_rt=23:00:00,h_data=${mem_per_core}G -pe shared ${ncore} temp/shell_scripts/run_buildgrn.sh ${celltype} ${supercell_dir} ${out_dir}
done < $supercell_file
