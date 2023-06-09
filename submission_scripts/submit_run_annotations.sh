#!/bin/bash

# Run from root directory /path/to/atlas/ 


##### INPUTS #####

#Set mode
mode="test" # or default
dbs="/u/scratch/m/mikechen/pathway_databases/genesets/human" # or mouse

###################
# Create folders
mkdir pathway_annotations
mkdir pathway_annotations/intermediate_annotations/

modules_dir='gene_memberships/'
intermediate_dir='./pathway_annotations/intermediate_annotations/'
convertToHuman="False" # "True" or "False"



# Extract all cell types
supercell_dir="supercells"
supercell_file="supercells.txt"
cd ${supercell_dir}
if [ ! -f ${supercell_file} ]
then
    find ./ | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi
cd ../
if [[ {$mode} == "test" ]]
then
    PARAMETERS=('0.01')
    num_supercells=1
elif [[ {$mode} == "default" ]]
    PARAMETERS=('0.002' '0.005' '0.007' '0.01' '0.03' '0.05')
    num_supercells=$(cat ${supercell_dir}/${supercell_file} | wc -l)
else
    echo "mode argument must be test or default"
    exit
fi


echo $num_supercells

for i in "${PARAMETERS[@]}"
do
    qsub -t 1:${num_supercells} shell_scripts/run_annotations.sh ${supercell_dir}/${supercell_file} ${modules_dir} ${i} ${dbs} ${intermediate_dir}

done

