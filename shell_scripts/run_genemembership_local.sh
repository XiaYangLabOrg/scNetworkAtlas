
network_dir=$1
network_file=$2
out_dir=$3
min_module_size=$4
max_module_size=$5
network_ext=$6

while read line;
do
	outfile=./${out_dir}/${line}gene_membership.txt
	if [ -f "$outfile" ]; then
		continue
	else
		echo "ModuleBasedDimensionalityReduction.py started on: " `date `
		echo " "
		python3 temp/python_files/ModuleBasedDimensionalityReduction.py ${network_dir}/${line}${network_ext} $outfile $min_module_size $max_module_size
		echo "ModuleBasedDimensionalityReduction.py ended on: " `date `
		echo " "
		
	fi
done < $network_dir/${network_file}
