### HOW TO USE ###

# PREREQUISITES: scing miniconda environment exists (build from SCING repo) #

# BEFORE YOU BEGIN: fill out scing_config below with the appropriate details for your project #

# Usage: python run_pipeline.py <setup|cell_mapping|pseudobulking|build_grn|merge_networks|gene_membership|annotations|process_annotations> #

# Run all steps, including setup, in the directory of your new project #

# This will set up all project files and activate the environment necessary to run SCING #

### END OF HOW TO USE ###

import sys

# Fetch config from current dir
from config import base_dir, scing_config

# Set up
# -----------------------------------

def setup():
    # Set up project structure (copy over files)
    # -----------------------------------
    print("\nCopying in project files...")
    import shutil
    import os

    # TODO: clone from Git repo once stable version is released

    # copy in files needed to run SCING
    #source_base_dir = '~/scratch/dry_run_pipeline'
    source_base_dir = '..'
    source_base_dir = os.path.expanduser(source_base_dir)  # resolve the tilde character into the actual home directory (python's shutil doesn't do this automatically)
    destination_dir = '.' # current directory assumed to be new project location

    # these names should reflect directory names from the scNetworkAtlas GitHub repo
    dirs = ['python_files', 'shell_scripts', 'submission_scripts']
    try:
        for directory in dirs:
            shutil.copytree(f'{source_base_dir}/{directory}', f'{destination_dir}/{directory}', dirs_exist_ok=False) # don't allow overwriting an existing directory
    except FileExistsError:
        print(f"Directory '{destination_dir}' already exists.")
    except OSError as e:
        print(f"Error: {e.strerror}")
    # -----------------------------------
    # End of Set up (copy over files)
    
    
    # Make bookkeeping directories
    # -----------------------------------
    print("\nCreating bookkeeping directories...")

    # Make directories as needed
    os.makedirs("timing_info", exist_ok=True) # store timing metrics
    os.makedirs("jobout", exist_ok=True) # store job logs

    # -----------------------------------
    # End of Make directories

    
    # Activate SCING conda env
    # -----------------------------------
    print("\nActivating conda env...")

    def activate_conda_environment(conda_env):
        try:
            combined_command = f"source {scing_config['conda_init_script']} && conda activate {conda_env}"
            print(f"\nATTENTION: PLEASE RUN \n`{combined_command}`\n to activate conda environment {conda_env}")
        except Exception as e:
            print(f"Error: {e}")

    # if activation is not working, ensure the correct anaconda/miniconda installation is being used on your system
    activate_conda_environment('scing')

    # -----------------------------------
    # End of Set up conda env

    
# -----------------------------------
# End of Set up



# Run pipeline steps
# -----------------------------------

import os

# TODO: must ensure correct conda env before each step (prompt user)
def cell_mapping():
    config = scing_config['cell_mapping']
    cmd = f"bash submission_scripts/submit_run_cellmapping.sh {base_dir} {config['mapping_file']} {config['adata_dir']} {config['celltype_column']}"
    os.system(cmd)

def pseudobulking():
    config = scing_config['pseudobulking_v2']
    cmd = f"bash submission_scripts/submit_run_supercells_v2.sh {config['celltype_col']} {config['sample_col']} {config['recluster']}"
    os.system(cmd)

def build_grn():
    cmd = f"bash submission_scripts/submit_run_buildgrn.sh {scing_config['build_grn']['num_networks']}"
    os.system(cmd)

def merge_networks():
    cmd = 'bash submission_scripts/submit_run_merge.sh'
    os.system(cmd)

def gene_membership():
    pass

def annotations():
    pass

def process_annotations():
    pass

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 run_pipeline.py <setup|cell_mapping|pseudobulking|build_grn|merge_networks|gene_membership|annotations|process_annotations>")
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
        print("Usage: python3 run_pipeline.py <setup|cell_mapping|pseudobulking|build_grn|merge_networks|gene_membership|annotations|process_annotations>")
        sys.exit(1)

# -----------------------------------
# End of Run pipeline steps

## TODO: create a clean system to delete all pipeline steps (like the actual python_files dir) once the pipeline is done running and user only wants outputs. Basically delete all non-generated files
