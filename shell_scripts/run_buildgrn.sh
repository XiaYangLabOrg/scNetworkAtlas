#!/bin/bash
#$ -cwd -V -N RUNNETWORK -l h_data=16G,h_rt=23:00:00 -j y -o jobout/BUILDGRN.$JOB_ID
. /u/local/Modules/default/init/modules.sh

# module load anaconda3
# conda activate scing

iteration=$SGE_TASK_ID
supercell=$1

echo "${supercell} ${iteration}"
supercell_dir="supercells"
file=saved_networks/intermediate_data/${supercell}/${supercell}.network.${iteration}.csv.gz
if [ -f "$file" ]
then
	sleep 10m
	exit
else 
	python3 python_files/BuildNetwork.py ${supercell_dir}/${supercell}.h5ad ${iteration}
	sleep 5m
	exit
fi
