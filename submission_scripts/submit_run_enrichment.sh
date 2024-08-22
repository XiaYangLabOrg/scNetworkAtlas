#!/bin/bash

# Run from root directory /path/to/atlas/
##### INPUTS #####
module_dir=$1
module_file=$2
out_dir=$3
pathway=$4
pathway_size_min=$5
pathway_size_max=$6
pathway_col=$7
module_col=$8
##################


cd $module_dir
if [ ! -f ${module_file} ]
then
    find ./ | grep txt | awk -F'.txt' '{print $1}' > ${module_file}
fi
cd ../
num_modules=$(cat ${module_dir}/${module_file} | wc -l)
echo $num_modules
mkdir $out_dir

qsub -t 1:${num_modules} temp/shell_scripts/run_enrichment.sh ${module_dir} ${module_file} ${out_dir} ${pathway} ${pathway_size_min} ${pathway_size_max} ${pathway_col} ${module_col}

