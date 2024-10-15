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

python3 ../temp/python_files/CellMapping.py $mapping_file $adata_dir $celltype_column

echo "sleeping"
sleep 5m

