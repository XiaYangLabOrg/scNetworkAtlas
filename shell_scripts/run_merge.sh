#!/bin/bash

#$ -cwd -V -N MERGENETWORKS -l h_data=16G,h_rt=23:00:00 -pe shared 4 -M eplau -m a -j y -o jobout/jobout.$JOB_ID
counter=1
supercell_dir=$1
supercell_file=$2
# . /u/local/Modules/default/init/modules.sh

# module load anaconda3
# conda activate scing

while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		file=saved_networks/final_edges/${line}.csv.gz

		# file=/u/home/j/juliemt/juliemt/SCING/human_atlas/saved_networks/final_edges/${line}.csv.gz
		if [ -f "$file" ]; then
			sleep 10m
			exit
		else
			python3 python_files/MergeNetworks.py ${supercell_dir}/${line}.h5ad saved_networks/intermediate_data/${line}/
			sleep 5m
			exit
		fi
		
	fi
	counter=$((${counter}+1)) 
done < ${supercell_dir}/${supercell_file}
	
