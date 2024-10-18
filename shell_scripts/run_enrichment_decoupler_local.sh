
module_dir=$1
module_file=$2
module_name_col=$3
module_gene_col=$4
pathway=$5
min_overlap=$6
pathway_size_min=$7
pathway_size_max=$8
n_background_genes=20000
out_dir=$9
jobout=${10}

while read line;
do
	outfile=./${out_dir}/${line}.${pathway}.enrichment.txt
	if [ -f "$outfile" ]; then
		continue
	else
		echo "Job $JOB_ID enrichment.py started on: " `date `
		echo " "
		python3 temp/python_files/enrichment.py \
			./${module_dir}/${line}.txt \
			--module_name_col $module_name_col \
			--module_gene_col $module_gene_col \
			--use_msigdb \
			--pathway $pathway \
			--min_overlap $min_overlap \
			--pathway_size_min ${pathway_size_min} \
			--pathway_size_max ${pathway_size_max} \
			--n_background_genes $n_background_genes \
			--out_file ${outfile}
		echo "Job $JOB_ID enrichment.py ended on: " `date `
		echo " "
	fi
	
done < $module_dir/${module_file}
	
