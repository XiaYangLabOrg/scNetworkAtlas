
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
jobout=$7

while read line;
do
	outfile=./${out_dir}/${line}gene_membership.txt
	if [ -f "$outfile" ]; then
		continue
	else
		echo "ModuleBasedDimensionalityReduction.py started on: " `date ` >> $jobout
		echo " "  >> $jobout
		python3 temp/python_files/ModuleBasedDimensionalityReduction.py ${network_dir}/${line}${network_ext} $outfile $min_module_size $max_module_size | tee -a $jobout
		echo "ModuleBasedDimensionalityReduction.py ended on: " `date `  >> $jobout
		echo " "  >> $jobout
		
	fi
done < $network_dir/${network_file}
