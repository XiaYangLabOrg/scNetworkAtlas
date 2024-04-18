### HOW TO USE ###

# PREREQUISITES: scing miniconda environment exists (build from SCING repo) #

# BEFORE YOU BEGIN: fill out scing_config below with the appropriate details for your project #

# Usage: python run_pipeline.py <setup|cell_mapping|pseudobulking|build_grn|merge_networks|gene_membership|annotations|process_annotations> #

# Run all steps, including setup, in the directory of your new project #

# This will set up all project files and activate the environment necessary to run SCING #

### END OF HOW TO USE ###

import sys

# User-configurable inputs
# -----------------------------------

# Configuration settings (these settings are read into submission scripts)
base_dir = "/u/project/xyang123/shared/reference/single_cell_databases/"
scing_config = {
    'cell_mapping': {
        'base_dir': base_dir,
        'mapping_file': f"${base_dir}all_celltypes/mouse_TabMurSenis_FACS.cleaned.tsv",
        'adata_dir': f"${base_dir}mouse/Tab_Mur_Senis_FACS/adatas/",
        'celltype_column': "celltypes"
    },
    'pseudobulking': {
        'filetype': 'npz' # npz or h5ad
    },
    'build_grn': {
        'num_networks': 100
    },
    'merge_networks': {
        'config': None # no inputs needed for merge networks step
    },
    'gene_membership': {
        'q1_module_sizes': ('20' '35' '50')
    },
    'annotations': {
        'rerun': "False",
        'mode': "test", # test or default
        'dbs': "/u/scratch/m/mikechen/pathway_databases/genesets/human_toy", # human or mouse
        'num_cores': 12,
        'q1_module_sizes': ('20' '35' '50'), # should match gene_membership config
        'modules_dir': 'gene_memberships/',
        'intermediate_dir': './pathway_annotations/intermediate_annotations'
    },
    'process_annotations': {
        'mode': 'test', # test or default
        'intermediate_dir': './pathway_annotations/intermediate_annotations',
        'final_dir': './pathway_annotations/final_annotations'
    },
}
    
# -----------------------------------
# End of User-configurable inputs



# Set up
# -----------------------------------

def setup():
    # Set up project structure (copy over files)
    # -----------------------------------
    import shutil

    # TODO: clone from Git repo once stable version is released

    # copy in files needed to run SCING
    source_base_dir = '~/scratch/dry_run_pipeline/'
    destination_dir = '.' # current directory assumed to be new project location

    # these names should reflect directory names from the scNetworkAtlas GitHub repo
    dirs = ['python_files', 'shell_scripts', 'submission_scripts']
    try:
        for directory in dirs:
            shutil.copytree(f'source_base_dir/{directory}', f'{destination_dir}/{directory}', dirs_exist_ok=False) # don't allow overwriting an existing directory
    except FileExistsError:
        print(f"Directory '{destination_dir}' already exists.")
    except OSError as e:
        print(f"Error: {e.strerror}")
    # -----------------------------------
    # End of Set up (copy over files)
    
    
    # Make bookkeeping directories
    # -----------------------------------

    import os

    # Make directories as needed
    os.makedirs("timing_info") # store timing metrics
    os.makedirs("jobout") # store job logs

    # -----------------------------------
    # End of Make directories

    
    # Activate SCING conda env
    # -----------------------------------

    import subprocess

    def activate_conda_environment(conda_env):
        try:
            # Activate the specified conda environment
            activate_cmd = f"conda activate {conda_env}"
            subprocess.run(activate_cmd, shell=True, check=True)
            print(f"Activated conda environment: {conda_env}")
        except subprocess.CalledProcessError:
            print(f"Error: Failed to activate conda environment {conda_env}")

    # if activation is not working, ensure the correct anaconda/miniconda installation is being used on your system
    activate_conda_environment('scing')

    # -----------------------------------
    # End of Set up conda env

    
# -----------------------------------
# End of Set up



# Run pipeline steps
# -----------------------------------

def cell_mapping():
    pass

def pseudobulking():
    pass

def build_grn():
    pass

def merge_networks():
    pass

def gene_membership():
    pass

def annotations():
    pass

def process_annotations():
    pass

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_pipeline.py <setup|cell_mapping|pseudobulking|build_grn|merge_networks|gene_membership|annotations|process_annotations>")
        sys.exit(1)
    
    step = sys.argv[1]

    if step == "setup":
        setup()
    elif step == "cell_mapping":
        cell_mapping()
    elif step == "pseudobulking":
        pseudobulking()
    elif step == "build_grn":
        build_grn()
    elif step == "merge_networks":
        merge_networks()
    elif step == "gene_membership":
        gene_membership()
    elif step == "annotations":
        annotations()
    elif step == "process_annotations":
        process_annotations()
    else:
        print("Invalid step:", step)
        print("Usage: python run_pipeline.py <setup|cell_mapping|pseudobulking|build_grn|merge_networks|gene_membership|annotations|process_annotations>")
        sys.exit(1)

# -----------------------------------
# End of Run pipeline steps

## TODO: create a clean system to delete all pipeline steps (like the actual python_files dir) once the pipeline is done running and user only wants outputs. Basically delete all non-generated files
