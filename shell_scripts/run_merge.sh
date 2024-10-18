#!/bin/bash

#$ -cwd -V -N MERGENETWORKS -j y -o jobout/MERGE.$JOB_ID

supercell_dir=$1
supercell_file=$2
consensus=$3
intermediate_dir=$4
out_dir=$5

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		out_file=$out_dir/${line}.${consensus}.txt
		if [ -f "$out_file" ]; then
			sleep 10m
			exit
		else
			echo "MergeNetworks.py started on: " `date `
			echo " "
			python temp/python_files/MergeNetworks.py ${supercell_dir}/${line}.h5ad ${intermediate_dir}/${line}/ $out_file --consensus $consensus
			echo "Job $JOB_ID end on: " `date `
			echo " "
			sleep 5m
			exit
		fi
	fi
	counter=$((${counter}+1)) 
done < ${supercell_file}
	
