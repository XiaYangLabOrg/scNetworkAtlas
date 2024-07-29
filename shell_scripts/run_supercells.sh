#!/bin/bash

#$ -cwd -V -N SUPERCELLS -l h_data=32G,h_rt=23:00:00 -j y -o jobout/supercells.$JOB_ID

tissue_celltype_file=$1
supercell_dir=$2
celltype_col=$3
filetype=$4
counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $line
		file=../${supercell_dir}/${line}.${filetype}
		if [ -f "$file" ]; then
			sleep 10m
			exit
		else
			python3 temp/python_files/BuildSupercells.py ${line}${filetype} ${supercell_dir} --stratify_by $celltype_col
			# python3 temp/python_files/BuildSupercells.py ${line}${filetype} ${supercell_dir}
			echo "sleeping"
			sleep 5m 
			exit
		fi
		break
	fi
	counter=$((${counter}+1)) 
done < ${tissue_celltype_file}

