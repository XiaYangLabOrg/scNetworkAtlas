#!/bin/bash
#$ -cwd -V -N RUNMAPPING -l h_data=8G,h_rt=23:00:00 -pe shared 1 -M eplau -m ea -j y -o jobout/cell_mapping.$JOB_ID -t 1:1

# . /u/local/Modules/default/init/modules.sh
# module load anaconda3
# conda activate scing

mapping_file=$1
adata_dir=$2
celltype_column=$3
tissue_dir=$4
cd $tissue_dir

start_time=$(date +%s) # Record start time

python3 ../python_files/CellMapping.py $mapping_file $adata_dir $celltype_column

end_time=$(date +%s) # Record end time
duration=$((end_time - start_time)) # Calculate duration

timing_file="timing_info/cell_mapping.txt"
echo "Start Time: $start_time" > $timing_file
echo "End Time: $end_time" >> $timing_file
echo "Duration: $duration seconds" >> $timing_file

echo "sleeping"
sleep 5m

