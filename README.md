# SCING network building pipeline

## Installation
```
git clone https://github.com/XiaYangLabOrg/scNetworkAtlas.git
cd scNetworkAtlas
conda env create -n scing --file install/scing.environment.yml 
conda activate scing
pip install pyitlib  
```

If you want to use the AUCell from SCENIC for graph based dimensionality reduction you must install pyscenic  
```
pip install pyscenic
```
If you want to perform pathway enrichment analysis with enrichr you must install R and enrichr OR use any R environment you have. Be sure that 'foreach' and 'doParallel' are installed in your R environment.
```  
conda env create -n enrichr --file install/pathway.environment.yml
conda activate enrichr
# install inside R
R
> install.packages('foreach')
> install.packages('doParallel')
> quit()
```

## Create a directory for each atlas
From scNetworkAtlas directory, make human/mouse directory and the atlas directory inside that.
```
mkdir human_atlas
mkdir human_atlas/Allen_10X
cd human_atlas/Allen_10X
cp -r /path/to/python/shell/submission/scripts ./
```

## Initial Directories
- python_scripts
- shell_scripts
- submission_scripts: requires changing input files
- ../../src: SCING scripts

Run all submission scripts from the root directory /path/to/atlas/


## 01. submit_run_cellmapping.sh
inputs: 
- base_dir: path to where all single cell database files are
- mapping_file: path to where the cell type mapping file is for your atlas
- adata_dir: path to where the adata file is for your atlas
- celltype_column: name of the adata column containing the original cell type labels (explore h5ad file in jupyter notebook beforehand)

```
base_dir="/u/project/xyang123/shared/reference/single_cell_databases/"
mapping_file="${base_dir}all_celltypes/human_Allen_10X.cleaned.tsv"
adata_dir="${base_dir}human/Allen_10X/adatas/"
celltype_column="celltypes"
```
After setting these inputs in the script, run it from the atlas directory with: <br>
```
bash submission_scripts/submit_run_cellmapping.sh
```

## 02. submit_run_supercells.sh
inputs:
- filetype: gene expression matrix input file type ("npz" or "h5ad")

```
bash submission_scripts/submit_run_supercells.sh
```

## 03. submit_run_buildgrn.sh
inputs:
- num_networks: number of networks (set to 100 for official run, set to 2 for test run)

After setting this input in the script, run it from the atlas directory with: <br>
```
bash submission_scripts/submit_run_buildgrn.sh
```

## 04. submit_run_merge.sh
No inputs required.
```
bash submission_scripts/submit_run_merge.sh
```

## 05. submit_run_genemembership.sh
inputs:
- q1_module_sizes: desired module sizes (set as the 25th percentile)
```
bash submission_scripts/submit_run_genemembership.sh
```

## 06. submit_run_annotations.sh
Must have R in your environment.
```
conda activate enrichr
```

inputs:
- rerun: set to "True" or "False" depending on whether you are rerunning jobs that did not finish in time.
- mode: set to "test" to run on one cell type at one parameter; set to "default" to run on all modules
- modules_dir: directory of gene memberships created in previous script relative to atlas root directory (e.g. gene_memberships)
- dbs: path to all pathway enrichment databases
- intermediate_dir: directory where files of annotations for individual modules for each resolution will be stored
- q1_module_sizes: desired module sizes (MUST MATCH from step 05)
```
bash submission_scripts/submit_run_annotations.sh
```

## 07. submit_run_processannotations.sh
Must have R in your environment.
```
conda activate enrichr
```
inputs:
- mode: set to "test" to run on one cell type at one parameter; set to "default" to run on all modules
- intermediate_dir: directory where files of annotations for individual modules for each resolution will be stored
- final_dir: directory where files of annotations of all modules within each resolution will be stored
```
bash submission_scripts/submit_run_processannotations.sh
```
