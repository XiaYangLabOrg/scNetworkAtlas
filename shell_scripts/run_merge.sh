#!/bin/bash

#$ -cwd -V -N MERGENETWORKS -l h_data=4G,h_rt=23:59:59 -pe shared 12 -M eplau -m a -j y -o jobout/MERGE.$JOB_ID
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
		out_file=$out_dir/${line}.${consensus}.csv.gz
		if [ -f "$out_file" ]; then
			sleep 10m
			exit
		else
			echo "MergeNetworks.py started on: " `date `
			echo " "
			python3 ../python_files/MergeNetworks.py ${supercell_dir}/${line}.h5ad ${intermediate_dir}/${line}/ $consensus $out_file
			echo "Job $JOB_ID end on: " `date `
			echo " "
			sleep 5m
			exit
		fi
	fi
	counter=$((${counter}+1)) 
done < ${supercell_file}
	
