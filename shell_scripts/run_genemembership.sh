
#$ -cwd -V -N GENEMEMBERSHIPS -l h_data=16G,h_rt=2:00:00 -M eplau -m a -j y -o jobout/gene_memberships.$JOB_ID

# . /u/local/Modules/default/init/modules.sh

#module load anaconda3
#conda activate scing

supercell_dir=$1
supercells=$2
q1_size=$3

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		outfile=./gene_memberships/${line}.gene_membership.${q1_size}.csv.gz
		if [ -f "$outfile" ]; then
			sleep 5m
			exit
		else
			python3 temp/python_files/ModuleBasedDimensionalityReduction.py ./saved_networks/final_edges/${line}.csv.gz ${outfile} ${q1_size}
			echo "sleeping"
			sleep 5m
			exit
		fi
		
	fi
	counter=$((${counter}+1)) 
done < $supercell_dir/${supercells}
	
