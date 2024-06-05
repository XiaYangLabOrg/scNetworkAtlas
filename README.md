## To Do

# SCING Pipeline

## Installing the Pipeline

First, clone this repository.

```
git clone https://github.com/XiaYangLabOrg/scNetworkAtlas.git
```

Now, checkout a new local branch that tracks the scGRNdb.v1 branch.

```
cd scNetworkAtlas
git checkout --track origin/scGRNdb.v1
```

Next, create a conda environment for SCING.

```
conda env create -n scing --file install/scing.environment.yml 
conda activate scing
```

TODO: I don't know why we need this package lol

```
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

## Creating Your Project Directory

Create a directory for your project.

```
mkdir <your_project_name>
cd <your_project_name>
```

Each project will be configured differently, so copy in a config file for this project.

```
cp ../config.py .
```

Also copy in the pipeline runner.

```
cp ../run_pipeline.py .
```

Note to those who have used previous versions of this repository:

- Each project directory no longer needs its own copy of submission scripts
- `run_pipeline.py` handles proper submission as long as your `config.py` is properly filled out

## Configuring Project Settings

Fill out `config.py` with your project and development environment details.

- You may ignore (leave as is) config settings for steps you won't run

## Running the SCING Pipeline

To run any step, `<step>` of the SCING pipeline, run

```
python3 run_pipeline.py <step>
```

where `<step>` must be one of 

- `cell_mapping`
- `pseudobulking_v2`
- `pseudobulking`
- `build_grn`
- `merge_networks`
- `gene_membership`
- `annotations`
- `process_annotations`

## Note to Those Peeking Into `submission_scripts`

You'll notice we run the shell scripts via a call to `../shell_scripts` even though `submission_scripts` and `shell_scripts` are in the same directory. This is because we intend for users to run the pipeline script from within their data directories (which are subdirectories of the cloned scNetworkAtlas repo). So, when a submission script is run, `shell_scripts` is in the parent directory of the current working directory.

Analogous reasoning for why `shell_scripts` call Python scripts via `../python_files`

## Old README Stuff

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
