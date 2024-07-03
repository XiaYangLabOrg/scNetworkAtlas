#!/bin/bash

# Run from root directory /path/to/atlas/
##### INPUTS #####
q1_module_sizes=('20' '35' '50')
supercell_dir="supercells"
supercell_file="supercells.txt"
##################

cd $supercell_dir
if [ ! -f ${supercell_file} ]
then
    find ./ | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi
cd ../
num_supercells=$(cat ${supercell_dir}/${supercell_file} | wc -l)
echo $num_supercells
mkdir gene_memberships
for i in "${q1_module_sizes[@]}"
do
	qsub -t 1:${num_supercells} temp/shell_scripts/run_genemembership.sh ${supercell_dir} ${supercell_file} ${i}
done

