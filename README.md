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
# Run only if you want to use the development version of the pipeline
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
- [`gene_membership`](#gene-membership)
- [`pathway_enrichment`](#pathway-enrichment)

### Setup

This step copies all submission, shell, and python scripts to your project directory.

Inputs:

- `main_branch_path`: path to scNetworkAtlas
- `base_dir`: for cell atlas projects, this is the cell atlas data directory. For other projects, this is your project directory

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

### Gene Membership

This step performs module detection in the final network.

Inputs:

- `network_dir`: directory of networks
- `network_file`: name of file to store all existing network paths. (If `network_ext` is txt, `network_file` must not end in txt)
- `network_ext`: file extension of network files (Does not start with a period. For example, use `txt` not `.txt`)
- `out_dir`: output directory
- `min_module_size`: minimum module size
- `max_module_size`: maximum module size
- `submit_command`: submit command (either "qsub" to run on cluster or "bash" to run locally)

### Pathway Enrichment

This step performs pathway enrichment on modules. It requires another conda environment called decoupler.

```bash
# ensure you are in base env
cd /path/to/scNetworkAtlas
conda env create --name decoupler install/decoupler_environment.yml
conda activate decoupler
cd /path/to/project_folder
```

If you want to use your own pathway database file, run the `enrichment` step and edit the `enrichment` parameters in the config file. **NOTE**: Make sure that the module gene symbols match the pathway gene symbols. If your pathway file does not match your module's gene symbol format, the decoupler API offers gene symbol conversion here: https://decoupler-py.readthedocs.io/en/latest/generated/decoupler.translate_net.html. You can convert your modules, and proceed with this step.

Inputs:

- `module_dir`: directory of modules
- `module_file`: name of file to store all existing module file paths (If module files are txt files, `module_file` must not end in txt)
- `module_name_col`: module column name in module file
- `module_gene_col`: gene column name in module file
- `pathway_file`: pathway database file. Must be 2 columns. If you want to include multiple pathway files, add them as a list (e.g. `['file1.txt','file2.txt']`)
- `pathway_db`: name of pathway database. Each pathway_file must have a pathway_db name. (e.g. `['pathway_db1','pathway_db2']`)
- `pathway_name_col`: pathway column name in pathway file
- `pathway_gene_col`: gene column name in pathway file
- `min_overlap`: minimum pathway-module overlap required for enrichment analysis
- `pathway_size_min`: minimum pathway size
- `pathway_size_max`: maximum pathway size
- `out_dir`: output directory
- `submit_command`: submit command (either "qsub" to run on cluster or "bash" to run locally)


If you do not have a pathway database file, the decoupler python package provides an API for msigdb pathway databases. Run the `enrichment_decoupler` step and edit the `enrichment_decoupler` parameters in the config file

Inputs:

- `module_dir`: directory of modules
- `module_file`: name of file to store all existing module file paths (If module files are txt files, `module_file` must not end in txt)
- `module_name_col`: module column name in module file
- `module_gene_col`: gene column name in module file
- `pathway`: pathway name in decoupler API. Possible names are listed beneath this argument in `config.py`. If you want to include multiple pathway files, add them as a list (e.g. `['pathway1','pathway2']`)
- `min_overlap`: minimum pathway-module overlap required for enrichment analysis
- `pathway_size_min`: minimum pathway size
- `pathway_size_max`: maximum pathway size
- `out_dir`: output directory
- `submit_command`: submit command (either "qsub" to run on cluster or "bash" to run locally)

