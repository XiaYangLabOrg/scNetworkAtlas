#!/bin/bash

# Run from root directory /path/to/Allen_10X/

# make final network directory
supercell_dir="supercells"
supercell_file="supercells/supercells.txt"
intermediate_dir="saved_networks/intermediate_data"
# if testing multiple consensus thresholds
# consensus_thresholds=(0.2 0.5 0.8 1.0)
# if running one consensus threshold
consensus_thresholds=0.5
out_dir="saved_networks/final_edges"

mkdir -p $out_dir
if [ ! -f ${supercell_file} ]
then
    ls $supercell_dir | grep h5ad | awk -F'.h5ad' '{print $1}' > ${supercell_file}
fi

num_lines=$(cat ${supercell_file} | wc -l)
echo $num_lines

for consensus in ${consensus_thresholds[@]}
do
    qsub -t 1-${num_lines}:1 ../shell_scripts/run_merge.sh $supercell_dir $supercell_file $consensus $intermediate_dir $out_dir
done

