#$ -cwd -V -N ANNOTATIONS -l h_data=16G,h_rt=23:00:00 -pe shared 2 -M eplau -m a -j y -o jobout/annotations.$JOB_ID

# . /u/local/Modules/default/init/modules.sh

#module load anaconda3
#conda activate enrichr

supercell_file=$1
modules_dir=$2
res=$3
db=$4
intermediate_dir=$5
convertToHuman=$6
# final_dir=$6

echo $supercell_file
echo $modules_dir
echo $res_param
echo $db
echo $intermediate_dir
# echo $final_dir

counter=1
while read line;
do
	if [ $SGE_TASK_ID = $counter ]
	then 
		echo $counter
		file="./pathway_annotations/intermediate_annotations/${line}/${res_param}_full.txt"
		if [ -f "$file" ];
		then
			echo "$file exists, sleeping"
			sleep 5m
			exit
		fi
		Rscript --vanilla ./r_files/PathwayEnrichmentAnnotation.r ${modules_dir} ${line} ${res} ${db} ${intermediate_dir} ${convertToHuman} # ${final_dir}
		echo "sleeping"
		sleep 5m
		exit
	fi
	counter=$((${counter}+1)) 
done < ./${supercell_file}

