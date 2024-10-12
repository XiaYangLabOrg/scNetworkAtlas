#!/bin/bash

# Run from root directory /path/to/atlas/
##### INPUTS #####
network_dir=$1
network_file=$2
out_dir=$3
min_module_size=$4
max_module_size=$5
network_ext=$6
submit_command=$7
##################

if [ ! -f ${network_dir}/${network_file} ]
then
    ls $network_dir | grep ${network_ext} | awk -F "${network_ext}" '{print $1}' > ${network_dir}/${network_file}
fi

num_networks=$(cat ${network_dir}/${network_file} | wc -l)
echo $num_networks
mkdir -p $out_dir

if [ $submit_command = "qsub" ]
then
    qsub -t 1:${num_networks} temp/shell_scripts/run_genemembership.sh ${network_dir} ${network_file} ${out_dir} ${min_module_size} ${max_module_size} ${network_ext}
elif [ $submit_command = "bash" ]
then
    mkdir -p "jobout/"
    bash temp/shell_scripts/run_genemembership_local.sh ${network_dir} ${network_file} ${out_dir} ${min_module_size} ${max_module_size} ${network_ext} "jobout/gene_memberships.$(date "+%m%d%Y.%s")"
else
    echo "submit_command must be qsub or bash, not ${submit_command}"
fi
