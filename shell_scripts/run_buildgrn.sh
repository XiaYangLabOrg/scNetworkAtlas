#!/bin/bash
#$ -cwd -V -N RUNNETWORK -j y -o jobout/BUILDGRN.$JOB_ID



iteration=$SGE_TASK_ID
supercell=$1
supercell_dir=$2
out_dir=$3
ncore=$4
mem_per_core=$5

echo "${supercell} ${iteration}"
out_file=${out_dir}/${supercell}/${supercell}.network.${iteration}.csv.gz
if [ -f "$out_file" ]
then
	sleep 10m
	exit
else 
	echo "BuildNetwork.py started on: " `date `
	echo " "
	python3 temp/python_files/BuildNetwork.py ${supercell_dir}/${supercell}.h5ad ${out_file} --seed ${iteration}
	echo "Job $JOB_ID end on: " `date `
	echo " "
	sleep 5m
	exit
fi
