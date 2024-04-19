#!/bin/bash

#$ -cwd -V -N SUPERCELLS -l h_data=32G,h_rt=23:00:00 -M eplau -m ea -j y -o jobout/supercells.$JOB_ID

# . /u/local/Modules/default/init/modules.sh
# module load anaconda3
# conda activate scing

tissue_celltype_file=$1
supercell_dir=$2
celltype_col=$3
sample_col=$4
recluster=$5
counter=1

while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $line
		file=${supercell_dir}/${line}.h5ad
		if [ -f "$file" ]; then
			sleep 10m
			exit
		else
			python3 python_files/BuildSupercells_v2.py ${line}.h5ad ${supercell_dir} ${celltype_col} ${sample_col} ${recluster}
			echo "sleeping"
			sleep 5m 
			exit
		fi
		break
	fi
	counter=$((${counter}+1)) 
done < ${tissue_celltype_file}


