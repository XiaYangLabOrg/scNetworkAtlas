#$ -cwd -V -N ANNOTATIONS -l h_data=4G,h_rt=23:00:00 -pe shared 2 -M eplau -m a -j y -o jobout/annotations.$JOB_ID

# . /u/local/Modules/default/init/modules.sh
#module load anaconda3
#conda activate scing

supercell_file=$1
intermediate_dir=$2
final_dir=$3

echo $supercell_file
echo $intermediate_dir
echo $final_dir

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
        start_time=$(date +%s) # Record start time
		file=${final_dir}/${line}/FILTERED.txt
		if [ -f "$file" ]; then
			sleep 120
			exit
		else
			Rscript --vanilla ./r_files/ProcessAnnotations.r ${line} ${intermediate_dir} ${final_dir}
            echo "ProcessAnnotations.r completed on: " `date`
            echo " "
            end_time=$(date +%s) # Record end time
            duration=$((end_time - start_time)) # Calculate duration
            echo "Total time taken: $duration seconds"
            echo "Job $JOB_ID end on: " `date`
            echo " "
            mkdir -p timing_info
            timing_file="../timing_info/process_annotations.txt"
            echo "Start Time: $start_time" > $timing_file
            echo "End Time: $end_time" >> $timing_file
            echo "Duration: $duration seconds" >> $timing_file
            sleep 5m
            exit
		fi
		
	fi
	counter=$((${counter}+1)) 
done < ./${supercell_file}

