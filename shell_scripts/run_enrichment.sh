
#$ -cwd -V -N ENRICHMENT -l h_data=16G,h_rt=2:00:00 -j y -o jobout/enrichment.$JOB_ID

# . /u/local/Modules/default/init/modules.sh

#module load anaconda3
#conda activate scing

module_dir=$1
module_file=$2
module_name_col=$3
module_gene_col=$4
pathway_file=$5
pathway_db=$6
pathway_name_col=$7
pathway_gene_col=$8
min_overlap=$9
pathway_size_min=${10}
pathway_size_max=${11}
n_background_genes=20000
out_dir=${12}

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		outfile=./${out_dir}/${line}.${pathway_db}.enrichment.txt
		if [ -f "$outfile" ]; then
			sleep 5m
			exit
		else
			echo "Job $JOB_ID enrichment.py started on: " `date `
			echo " "
			python3 temp/python_files/enrichment.py \
				./${module_dir}/${line}.txt \
				--module_name_col $module_name_col \
				--module_gene_col $module_gene_col \
				--pathway_file ${pathway_file} \
				--pathway_name_col $pathway_name_col \
				--pathway_gene_col $pathway_gene_col \
				--min_overlap $min_overlap \
				--pathway_size_min ${pathway_size_min} \
				--pathway_size_max ${pathway_size_max} \
				--n_background_genes $n_background_genes \
				--out_file ${outfile}
			
			echo "Job $JOB_ID enrichment.py ended on: " `date `
			echo " "
			echo "sleeping"
			sleep 5m
			exit
		fi
		
	fi
	counter=$((${counter}+1)) 
done < $module_dir/${module_file}
	
