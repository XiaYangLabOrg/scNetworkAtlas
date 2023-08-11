#$ -cwd -V -N ANNOTATIONS -l h_data=12G,h_rt=23:00:00 -M eplau -m a -j y -o jobout/annotations.$JOB_ID

# . /u/local/Modules/default/init/modules.sh

#module load anaconda3
#conda activate enrichr

supercell_file=$1
modules_dir=$2
q1=$3
db=$4
intermediate_dir=$5
convertToHuman=$6
num_cores=$7
# final_dir=$6

echo $supercell_file
echo $modules_dir
echo $q1
echo $db
echo $intermediate_dir
# echo $final_dir

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		file="./pathway_annotations/intermediate_annotations/${line}/${q1}_full.txt"
		if [ -f "$file" ];
		then
			echo "$file exists, sleeping"
			sleep 5m
			exit
		fi
		Rscript --vanilla ./r_files/PathwayEnrichmentAnnotation.r ${modules_dir} ${line} ${q1} ${db} ${intermediate_dir} ${convertToHuman} ${num_cores}
		echo "sleeping"
		sleep 5m
		exit
	fi
	counter=$((${counter}+1)) 
done < ./${supercell_file}

