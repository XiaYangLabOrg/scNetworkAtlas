#!/bin/bash

#$ -cwd -V -N MERGENETWORKS -l h_data=4G,h_rt=23:59:59 -pe shared 12 -M eplau -m a -j y -o jobout/MERGE.$JOB_ID
supercell_dir=$1
supercell_file=$2
consensus=$3
intermediate_dir=$4
out_dir=$5

start_time=$(date +%s) # Record start time

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
			python3 python_files/MergeNetworks.py ${supercell_dir}/${line}.h5ad ${intermediate_dir}/${line}/ $consensus $out_file
            
            python_exit_status=$? # Capture exit status of Python script
            if [ $python_exit_status -eq 0 ]; then
                echo "MergeNetworks.py completed successfully on: " `date`
            else
                echo "MergeNetworks.py failed with exit code $python_exit_status on: " `date`
            fi
            
            echo " "
            end_time=$(date +%s) # Record end time
            duration=$((end_time - start_time)) # Calculate duration
            echo "Total time taken: $duration seconds"
            echo "Job $JOB_ID end on: " `date`
            echo " "
            
            mkdir -p ../timing_info
            timing_file="../timing_info/merge.txt"
            echo "Start Time: $start_time" > $timing_file
            echo "End Time: $end_time" >> $timing_file
            echo "Duration: $duration seconds" >> $timing_file
            sleep 5m
            exit
		fi
	fi
	counter=$((${counter}+1)) 
done < ${supercell_file}
	
