# SCING Pipeline

## Installing the Pipeline

If you have not installed SCING, follow the installation instructions and conda environment setup in the SCING repo [here](https://github.com/XiaYangLabOrg/SCING). Be sure that SCING version is up to date with the following.

```bash
cd /path/to/SCING/repo
git pull origin
cd ../
```


Then, clone this repository.

```bash
git clone https://github.com/XiaYangLabOrg/scNetworkAtlas.git
cd scNetworkAtlas
git checkout --track origin/scGRNdb.v1
```

## Creating Your Project Directory

Create a directory for your project.

```bash
mkdir <your_project_name>
cd <your_project_name>
```

Each project will be configured differently, so copy in a config and run_pipeline file for this project.

```bash
cp /path/to/scNetworkAtlas/config.py .
cp /path/to/scNetworkAtlas/run_pipeline.py .
```

Note to those who have used previous versions of this repository:

## Configuring Project Settings

Fill out `config.py` with your project and development environment details.

- You may ignore (leave as is) config settings for steps you won't run

## Running the SCING Pipeline

To run any step, `<step>` of the SCING pipeline, run

```bash
python3 run_pipeline.py <step>
```

where `<step>` must be one of 

- [`setup`](#setup)
- [`cell_mapping`](#cell_mapping)
- [`pseudobulking`](#pseudobulking-supercells)
- [`build_grn`](#build_grn)
- [`merge_networks`](#merge-networks)

*In progress*

- `gene_membership`
- `annotations`
- `process_annotations`

### Setup

This step copies all submission, shell, and python scripts to your project directory.

Inputs:

- `main_branch_path`: path to scNetworkAtlas
- `base_dir`: for cell atlas projects, this is the cell atlas data directory. For other projects, this is your project directory
- `conda_init_script`: path to conda.sh, normally /path/to/miniconda3/etc/profile.d/conda.sh

### Cell_mapping

This step updates cell type labels in your single cell data. If you already have the cell type labels, you can skip this step.

Inputs:

- `base_dir`: same as setup `base_dir`
- `mapping_file`: tab-separated file with columns for `Original Cell Type` and - - `Broader Cell Type`
- `adata_dir`: path to adata directory
- `celltype_column`: cell type column in the single cell object

### Pseudobulking (Supercells)

This step performs the cell pseudobulking using leiden clustering.

Inputs:

- `tissue_dir`: path to adata directory
- `supercell_dir`: output pseudobuolk adata directory
- `filetype`: file type for counts data (h5ad or npz)
- `celltype_col`: cell type column in the single cell object
- `tissue_celltype_file`: name of txt file to store all existing adata paths

### Build_grn

This step builds intermediate GRNs by bootstrapping the pseudobulk cells.

Inputs:

- `num_networks`: number of intermediate networks
- `supercell_dir`: pseudobulk adata directory
- `supercell_file`: name of txt file to store all existing supercell file paths
- `out_dir`: output directory
- `ncore`: number of cores used to build each network (default is 1)
- `mem_per_core`: memory per core in GB (default is 16)

### Merge Networks

This step merges intermediate GRNs.

Inputs:

- `supercell_dir`: pseudobulk adata directory
- `supercell_file`: name of txt file to store all existing supercell file paths
- `intermediate_dir`: directory to intermediate networks
- `consensus`: list of consensus thresholds to test (default is [0.5])
- `out_dir`: output directory
- `ncore`: number of cores used to build each network (default is 12)
- `mem_per_core`: memory per core in GB (default is 4)

### Gene_membership (In Progress)

This step performs module detection in the final network.


### Pathway Enrichment (In Progress)

This step performs pathway enrichment on modules

<!-- 
## Updates in Progress

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
``` -->
