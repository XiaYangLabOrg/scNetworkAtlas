#!/bin/bash

# Run from root directory /path/to/atlas/ 


##### INPUTS #####
supercell_dir="supercells"
supercell_file="supercells.txt"
#Set mode
rerun="False" # or True
mode="test" # test or default
dbs="/u/scratch/m/mikechen/pathway_databases/genesets/human_toy" # human or mouse
# number of cores to run job
num_cores=12
# module size parameters (should match submit_run_genemembership.sh)
q1_module_sizes=('20' '35' '50')
###################
# Create folders
if [ $rerun = "True" ]
then
    mkdir -p temp_files
fi

modules_dir='gene_memberships/'
intermediate_dir='./pathway_annotations/intermediate_annotations'
if [[ "${mode}" == "test" ]]
then
    intermediate_dir="${intermediate_dir}_test"
fi
convertToHuman="False" # "True" or "False"
mkdir -p pathway_annotations
mkdir -p $intermediate_dir



# Extract all cell types
cd ${supercell_dir}
if [ ! -f ${supercell_file} ]
then
    find ./ | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi
cd ../
if [[ "${mode}" == "test" ]]
then
    q1_module_sizes=(${q1_module_sizes[0]})
    num_supercells=1
    num_cores=4
elif [[ "${mode}" == "default" ]]
then
    num_supercells=$(cat ${supercell_dir}/${supercell_file} | wc -l)
else
    echo "mode argument must be test or default"
    exit
fi

echo $num_supercells

for i in "${q1_module_sizes[@]}"
do
    # rerun on supercells without pathway results
    if [ $rerun = "True" ]
    then
        param_subset_file="temp_files/supercells_to_run_${i}.txt"
        rm -f $param_subset_file # remove if exists
        while read supercell
        do
            if [ ! -f ${intermediate_dir}/${supercell}/${i}_full.txt ]
            then
                echo "$supercell" >> ${param_subset_file}
            else
                # remove extraneous files
                ls ${intermediate_dir}/${supercell}/${i}_* | grep -v "full.txt" | xargs rm
            fi
        done < ${supercell_dir}/${supercell_file}
                
        if [ -f ${param_subset_file} ]
        then
            # number of supercells without results
            num_rerun_supercells=$(cat ${param_subset_file} | wc -l)
            echo "Rerunning ${i} module size on ${num_rerun_supercells} celltypes"
            cat $param_subset_file
            qsub -t 1:${num_rerun_supercells} -pe shared ${num_cores} ../shell_scripts/run_annotations.sh ${param_subset_file} ${modules_dir} ${i} ${dbs} ${intermediate_dir} ${convertToHuman} ${num_cores}
            
        else
            echo "No need to rerun ${i} module size"
            
        fi  
    else
        qsub -t 1:${num_supercells} -pe shared ${num_cores} ../shell_scripts/run_annotations.sh ${supercell_dir}/${supercell_file} ${modules_dir} ${i} ${dbs} ${intermediate_dir} ${convertToHuman} ${num_cores}
    fi
    echo -e "\n"
done
