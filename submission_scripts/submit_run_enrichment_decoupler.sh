#!/bin/bash

# Run from root directory /path/to/atlas/
##### INPUTS #####
module_dir=$1
module_file=$2
module_name_col=$3
module_gene_col=$4
pathway=$5
min_overlap=$6
pathway_size_min=$7
pathway_size_max=$8
out_dir=$9
submit_command=${10}
##################

IFS=',' read -r -a pathway_array <<< $pathway
if [ ! -f ${module_dir}/${module_file} ]
then
    ls ${module_dir}| grep txt | awk -F'.txt' '{print $1}' > ${module_dir}/${module_file}
fi

num_modules=$(cat ${module_dir}/${module_file} | wc -l)
echo $num_modules
mkdir -p $out_dir

for p in ${pathway_array[@]}
do
    if [ $submit_command = "qsub" ]
    then
        qsub -t 1:${num_modules} temp/shell_scripts/run_enrichment_decoupler.sh $module_dir $module_file $module_name_col $module_gene_col $p $min_overlap $pathway_size_min $pathway_size_max $out_dir
    elif [ $submit_command = "bash" ]
    then
        mkdir -p "jobout/"
        bash temp/shell_scripts/run_enrichment_decoupler_local.sh $module_dir $module_file $module_name_col $module_gene_col $p $min_overlap $pathway_size_min $pathway_size_max $out_dir > "jobout/enrichment.${p}.$(date "+%m%d%Y.%s")" 2>&1 &
    fi
done
