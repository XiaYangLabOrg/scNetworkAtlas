# SCING network building pipeline

## Create a directory for each atlas
From SCING directory, make human/mouse directory and the atlas directory inside that.
> mkdir human_atlas
> mkdir human_atlas/Allen_10X
> cd human_atlas/Allen_10X
> cp -r /path/to/python/shell/submission/scripts ./

## Initial Directories
- python_scripts
- shell_scripts
- submission_scripts: requires changing input files
- ../../src: SCING scripts

Run all submission scripts from the root directory /path/to/atlas/


## 01. submit_run_cell_mapping.sh
inputs: 
- base_dir: path to where all single cell database files are
- mapping_file: path to where the cell type mapping file is for your atlas
- adata_dir: path to where the adata file is for your atlas
- celltype_column: name of the adata column containing the original cell type labels (explore h5ad file in jupyter notebook beforehand)

> base_dir="/u/project/xyang123/shared/reference/single_cell_databases/"
> mapping_file="${base_dir}all_celltypes/human_Allen_10X.cleaned.tsv"
> adata_dir="${base_dir}human/Allen_10X/adatas/"
> celltype_column="celltypes"

After setting these inputs in the script, run it from the atlas directory with: <br>
> bash submission_scripts/submit_run_cell_mapping.sh

## 02. submit_run_supercells.sh
No inputs required. Must run submit_run_cell_mapping.sh before this.
> bash submission_scripts/submit_run_supercells.sh

## 03. submit_run_buildgrn.sh
inputs:
- num_networks: number of networks (set to 100 for official run, set to 2 for test run)
- merged_adata_dir: directory where supercells are

After setting this input in the script, run it from the atlas directory with: <br>
> bash submission_scripts/submit_run_buildgrn.sh

## 04. submit_run_merge.sh
No inputs required.
> bash submission_scripts/submit_run_merge.sh

## 05. submit_run_genemembership.sh
inputs:
- PARAMETERS: different resolutions at which to create the modules
- supercell_dir: supercell directory relative to atlas root directory (e.g. merged_adata)
- supercell_file: name of supercell file (e.g. supercells.txt)
- num_supercells: generated from supercell file
> bash submission_scripts/submit_run_genemembership.sh

## 06. submit_run_annotations.sh
inputs:
- mode: set to "test" to run on one cell type at one parameter; set to "default" to run on all modules
- modules_dir: directory of gene memberships created in previous script relative to atlas root directory (e.g. gene_memberships)
- dbs: path to all pathway enrichment databases
- intermediate_dir: directory where files of annotations for individual modules for each resolution will be stored
- final dir: directory where files of annotations of all modules within each resolution will be stored
- PARAMETERS: different resolutions at which to create the modules
- convertToHuman: whether the atlas is from mouse data (put True if yes, False if no)
> bash submission_scripts/submit_run_annotations.sh

## 07. submit_run_processannotations.sh
inputs:
- mode: set to "test" to run on one cell type at one parameter; set to "default" to run on all modules
- intermediate_dir: directory where files of annotations for individual modules for each resolution will be stored
- final dir: directory where files of annotations of all modules within each resolution will be stored
> bash submission_scripts/submit_run_processannotations.sh
