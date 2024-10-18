#!/bin/bash

# Run from root directory /path/to/atlas/
##### INPUTS #####
module_dir=$1
module_file=$2
module_name_col=$3
module_gene_col=$4
pathway_file=$5
pathway_db=$6
pathway_name_col=$7
pathway_gene_col=$8
min_overlap=$9
pathway_size_min=${10}
pathway_size_max=${11}
out_dir=${12}
submit_command=${13}
##################

IFS=',' read -r -a pathway_file_array <<< $pathway_file
IFS=',' read -r -a pathway_db_array <<< $pathway_db

if [ ! -f ${module_dir}/${module_file} ]
then
    ls ${module_dir}| grep txt | awk -F'.txt' '{print $1}' > ${module_dir}/${module_file}
fi

num_modules=$(cat ${module_dir}/${module_file} | wc -l)
echo $num_modules
mkdir -p $out_dir

num_pathways=${#pathway_file_array[@]}
for (( i=0; i<$num_pathways; i++ ))
do
    if [ $submit_command = "qsub" ]
    then
        qsub -t 1:${num_modules} temp/shell_scripts/run_enrichment.sh $module_dir $module_file $module_name_col $module_gene_col ${pathway_file_array[$i]} ${pathway_db_array[$i]} $pathway_name_col $pathway_gene_col $min_overlap $pathway_size_min $pathway_size_max $out_dir
    elif [ $submit_command = "bash" ]
    then
        mkdir -p "jobout/"
        bash temp/shell_scripts/run_enrichment_local.sh $module_dir $module_file $module_name_col $module_gene_col ${pathway_file_array[$i]} ${pathway_db_array[$i]} $pathway_name_col $pathway_gene_col $min_overlap $pathway_size_min $pathway_size_max $out_dir > "jobout/enrichment.${pathway_db_array[$i]}.$(date "+%m%d%Y.%s")" 2>&1 &  
    fi
done

