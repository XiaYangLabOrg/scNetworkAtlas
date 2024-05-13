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
