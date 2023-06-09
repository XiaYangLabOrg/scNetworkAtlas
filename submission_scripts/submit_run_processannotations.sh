#!/bin/bash

# Run from root directory /path/to/atlas/ 

##### INPUTS #####

#Set mode
mode="test" # or default
###################

# Extract all cell types
supercell_dir="supercells/"
supercell_file="supercells.txt"
cd ${supercell_dir}
if [ ! -f ${supercell_file} ]
then
    find ./ | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi
cd ../


if [[ {$mode} == "test" ]]
then
    num_supercells=1
elif [[ {$mode} == "default" ]]
    num_supercells=$(cat ${supercell_dir}/${supercell_file} | wc -l)
else
    echo "mode argument must be test or default"
    exit
fi

echo $num_supercells


# Set user inputs
intermediate_dir='./pathway_annotations/intermediate_annotations'
final_dir='./pathway_annotations/final_annotations'
mkdir pathway_annotations/final_annotations/

qsub -t 1:${num_supercells} shell_scripts/run_processannotations.sh ${supercell_dir}/${supercell_file} ${intermediate_dir} ${final_dir}


