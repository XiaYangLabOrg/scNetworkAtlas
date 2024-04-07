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

start_time=$(date +%s) # Record start time

if [ -f "$file" ]
then
	sleep 10m
	exit
else 
	python3 python_files/BuildNetwork.py ${supercell_dir}/${supercell}.h5ad ${iteration}
    echo " "
    end_time=$(date +%s) # Record end time
    duration=$((end_time - start_time)) # Calculate duration
    echo "Total time taken: $duration seconds"
    echo "Job $JOB_ID end on: " `date`
    echo " "
    
    mkdir -p timing_info
    timing_file="timing_info/buildgrn.txt"
    echo "Start Time: $start_time" > $timing_file
    echo "End Time: $end_time" >> $timing_file
    echo "Duration: $duration seconds" >> $timing_file
    sleep 5m
    exit
fi
