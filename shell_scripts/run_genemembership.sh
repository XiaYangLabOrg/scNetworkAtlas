
#$ -cwd -V -N GENEMEMBERSHIPS -l h_data=16G,h_rt=2:00:00 -M eplau -m a -j y -o jobout/gene_memberships.$JOB_ID

# . /u/local/Modules/default/init/modules.sh

#module load anaconda3
#conda activate scing

supercell_dir=$1
supercells=$2
parameter=$3

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		file=./gene_memberships/${line}.gene_membership.${parameter}.csv.gz
		if [ -f "$file" ]; then
			sleep 5m
			exit
		else
			python3 ./python_files/ModuleBasedDimensionalityReduction.py $supercell_dir/${line}.h5ad ./saved_networks/final_edges/${line}.csv.gz ${line} ${parameter}
			sleep 5m
			exit
		fi
		
	fi
	counter=$((${counter}+1)) 
done < $supercell_dir/${supercells}
	
