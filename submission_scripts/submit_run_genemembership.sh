#!/bin/bash

# Run from root directory /path/to/atlas/
##### INPUTS #####
network_dir=$1
network_file=$2
out_dir=$3
min_module_size=$4
max_module_size=$5
network_ext=$6
##################

if [ ! -f ${network_file} ]
then
    ls $network_dir | grep ${network_ext} | awk -F "${network_ext}" '{print $1}' > ${network_dir}/${network_file}
fi

num_networks=$(cat ${network_dir}/${network_file} | wc -l)
echo $num_networks
mkdir -p $out_dir

qsub -t 1:${num_networks} temp/shell_scripts/run_genemembership.sh ${network_dir} ${network_file} ${out_dir} ${min_module_size} ${max_module_size} ${network_ext}

