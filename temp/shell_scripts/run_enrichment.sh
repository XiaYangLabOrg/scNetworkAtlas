
#$ -cwd -V -N ENRICHMENT -l h_data=16G,h_rt=2:00:00 -j y -o jobout/enrichment.$JOB_ID

# . /u/local/Modules/default/init/modules.sh

#module load anaconda3
#conda activate scing

module_dir=$1
module_file=$2
out_dir=$3
pathway=$4
pathway_size_min=$5
pathway_size_max=$6
pathway_col=$7
module_col=$8

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		outfile=./${out_dir}/${line}.enrichment.txt
		if [ -f "$outfile" ]; then
			sleep 5m
			exit
		else
			python3 temp/python_files/enrichment.py ./${module_dir}/${line}.txt ${outfile} ${pathway} ${pathway_size_min} ${pathway_size_max} ${pathway_col} ${module_col}
			echo "sleeping"
			sleep 5m
			exit
		fi
		
	fi
	counter=$((${counter}+1)) 
done < $module_dir/${module_file}
	
