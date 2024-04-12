#!/bin/bash

#$ -cwd -V -N SUPERCELLS -l h_data=32G,h_rt=23:00:00 -M eplau -m ea -j y -o jobout/supercells.$JOB_ID

# . /u/local/Modules/default/init/modules.sh
# module load anaconda3
# conda activate scing

tissue_dir=$1
tissue_celltype_file=$2
supercell_dir=$3
filetype=$4
counter=1
cd $tissue_dir
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $line
        start_time=$(date +%s) # Record start time
		file=../${supercell_dir}/${line}.h5ad
		if [ -f "$file" ]; then
			sleep 10m
			exit
		else
			python3 ../python_files/BuildSupercells.py ${line}.${filetype}
            echo "BuildSupercells.py completed on: " `date`
            echo " "
            end_time=$(date +%s) # Record end time
            duration=$((end_time - start_time)) # Calculate duration
            echo "Total time taken: $duration seconds"
            echo "Job $JOB_ID end on: " `date`
            echo " "
            mkdir -p timing_info
            timing_file="../timing_info/supercells.txt"
            echo "Start Time: $start_time" > $timing_file
            echo "End Time: $end_time" >> $timing_file
            echo "Duration: $duration seconds" >> $timing_file
            echo "sleeping"
            sleep 5m
            exit
		fi
		break
	fi
	counter=$((${counter}+1)) 
done < ${tissue_celltype_file}


