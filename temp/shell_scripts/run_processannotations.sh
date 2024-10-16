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
		file=${final_dir}/${line}/FILTERED.txt
		if [ -f "$file" ]; then
			sleep 120
			exit
		else
			Rscript --vanilla temp/r_files/ProcessAnnotations.r ${line} ${intermediate_dir} ${final_dir}
			sleep 5m
			exit
		fi
		
	fi
	counter=$((${counter}+1)) 
done < ./${supercell_file}

