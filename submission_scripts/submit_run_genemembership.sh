#!/bin/bash

# Run from root directory /path/to/atlas/
##### INPUTS #####
PARAMETERS=('0.002' '0.005' '0.007' '0.01' '0.03' '0.05')
##################

supercell_dir="supercells"
supercell_file="supercells.txt"
cd $supercell_dir

if [ ! -f ${supercell_file} ]
then
    find ./ | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi
cd ../

num_supercells=$(cat ${supercell_dir}/${supercell_file} | wc -l)
echo $num_supercells

mkdir gene_memberships

for i in "${PARAMETERS[@]}"
do
	parameter="$i"
	qsub -t 1:${num_supercells} shell_scripts/run_genemembership.sh ${supercell_dir} ${supercell_file} ${parameter}
done

