
#$ -cwd -V -N GENEMEMBERSHIPS -l h_data=16G,h_rt=2:00:00 -j y -o jobout/gene_memberships.$JOB_ID

# . /u/local/Modules/default/init/modules.sh

#module load anaconda3
#conda activate scing

network_dir=$1
network_file=$2
out_dir=$3
min_module_size=$4
max_module_size=$5
network_ext=$6

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		outfile=./${out_dir}/${line}gene_membership.txt
		if [ -f "$outfile" ]; then
			sleep 5m
			exit
		else
			python3 temp/python_files/ModuleBasedDimensionalityReduction.py ${network_dir}/${line}${network_ext} $outfile $min_module_size $max_module_size
			echo "sleeping"
			sleep 5m
			exit
		fi
		
	fi
	counter=$((${counter}+1)) 
done < $network_dir/${network_file}
	
