### HOW TO USE ###

# PREREQUISITES: scing miniconda environment exists (build from SCING repo) #

# BEFORE YOU BEGIN: fill out config.py with the appropriate details for your project #

# Run all steps in the directory of your project directory #

### END OF HOW TO USE ###

USAGE = "python run_pipeline.py <cell_mapping|pseudobulking_v2|pseudobulking|build_grn|merge_networks|gene_membership|annotations|process_annotations|clean>"

import sys

# Fetch config from current dir
from config import base_dir, scing_config

def activate_conda_environment(conda_env):
    try:
        combined_command = f"source {scing_config['conda_init_script']} && conda activate {conda_env}"
        print(f"\nATTENTION: PLEASE RUN \n`{combined_command}`\n to activate conda environment {conda_env}")
    except Exception as e:
        print(f"Error: {e}")

# Set up
# -----------------------------------

# these names should reflect directory names from the scNetworkAtlas GitHub repo
dirs = ['python_files', 'shell_scripts', 'submission_scripts', 'r_files']

def setup():
    # Set up project structure (copy over files)
    # -----------------------------------
#     print("\nCopying in project files...")
#     import shutil
#     import os

#     # copy in files needed to run SCING
#     source_base_dir = '..'
#     source_base_dir = os.path.expanduser(source_base_dir)  # resolve the tilde character into the actual home directory (python's shutil doesn't do this automatically)
#     destination_dir = '.' # current directory assumed to be new project location
    
#     response = input(f"Are you sure you want to rerun setup? This will overwrite files in {dirs} and delete any modifications you made. (y/n): ").strip().lower()
#     if response == 'y':
#         print("Continuing setup...")
#     elif response == 'n':
#         print("Setup aborted.")
#         return
#     else:
#         print("Invalid input. Please enter 'y' or 'n'. Exiting...")
#         return
    
#     # copy in files
#     try:
#         for directory in dirs:
#             shutil.copytree(f'{source_base_dir}/{directory}', f'{destination_dir}/{directory}', dirs_exist_ok=True)
#     except FileExistsError:
#         print(f"Directory '{destination_dir}' already exists.")
#     except OSError as e:
#         print(f"Error: {e.strerror}")
#     # -----------------------------------
#     # End of Set up (copy over files)
    
    
#     # Make bookkeeping directories
#     # -----------------------------------
# #     print("\nCreating bookkeeping directories...")

# #     # Make directories as needed
# #     os.makedirs("timing_info", exist_ok=True) # store timing metrics
# #     os.makedirs("jobout", exist_ok=True) # store job logs

#     # -----------------------------------
#     # End of Make directories

    
#     # Activate SCING conda env
#     # -----------------------------------
# #     print("\nActivating conda env...")

#     # if activation is not working, ensure the correct anaconda/miniconda installation is being used on your system
#     activate_conda_environment('scing')

#     # -----------------------------------
#     # End of Set up conda env
    pass
    
# -----------------------------------
# End of Set up



# Run pipeline steps
# -----------------------------------

import os

def confirm_conda_activated():
    confirmation = input("Is your scing conda environment activated? (y/n): ").lower()
    if confirmation == 'y':
        return True
    elif confirmation == 'n':
        print("Please activate your conda environment and run the script again.")
        activate_conda_environment('scing')
        sys.exit()
    else:
        print("Invalid input. Please enter 'y' or 'n'.")
        return False

def cell_mapping():
    confirm_conda_activated()

    config = scing_config['cell_mapping']
    cmd = f"bash ../submission_scripts/submit_run_cellmapping.sh {base_dir} {config['mapping_file']} {config['adata_dir']} {config['celltype_column']}"
    os.system(cmd)

def pseudobulking_v2():
    confirm_conda_activated()

    config = scing_config['pseudobulking_v2']
    cmd = f"bash ../submission_scripts/submit_run_supercells_v2.sh {config['tissue_dir']} {config['supercell_dir']} {config['celltype_col']} {config['sample_col']} {config['tissue_celltype_file']}"
    os.system(cmd)

def pseudobulking():
    confirm_conda_activated()

    config = scing_config['pseudobulking']
    cmd = f"bash ../submission_scripts/submit_run_supercells.sh {config['tissue_dir']} {config['supercell_dir']} {config['filetype']} {config['celltype_col']} {config['tissue_celltype_file']}"
    os.system(cmd)
    
def build_grn():
    confirm_conda_activated()
    config = scing_config['build_grn']
    cmd = f"bash ../submission_scripts/submit_run_buildgrn.sh {config['num_networks']} {config['supercell_dir']} {config['supercell_file']} {config['out_dir']}"
    os.system(cmd)

def merge_networks():
    confirm_conda_activated()
    config = scing_config['merge_networks']

    # pass as string to submission script
    consensus_str = ','.join(map(str, scing_config['merge_networks']['consensus_thresholds']))

    cmd = f"bash ../submission_scripts/submit_run_merge.sh {config['supercell_dir']} {config['supercell_file']} {config['intermediate_dir']} {consensus_str} {config['out_dir']}"
    os.system(cmd)

def gene_membership():
    confirm_conda_activated()

    pass

def annotations():
    confirm_conda_activated()

    pass

def process_annotations():
    confirm_conda_activated()

    pass

def clean():
    pass
    # import shutil
    
    # confirmation = input(f"Please confirm you'd like to delete all pipeline scripts (does not delete any outputs produced by SCING, just scripts in {dirs} used to run scripts). (y/n): ").lower()
    
    # if confirmation == 'y':
    #     try:
    #         for d in dirs:
    #             shutil.rmtree(d)
    #             print(f"Directory '{d}' successfully deleted.")
    #     except FileNotFoundError:
    #         print(f"Directory '{d}' not found.")
    #     except OSError as e:
    #         print(f"Error: {e}")
    # elif confirmation == 'n':
    #     sys.exit()
    # else:
    #     print("Invalid input. Please retry.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(USAGE)
        sys.exit(1)
    
    step = sys.argv[1]

    if step == "setup":
        setup()
    elif step == "cell_mapping":
        cell_mapping()
    elif step == "pseudobulking_v2":
        pseudobulking_v2()
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
    elif step == "clean":
        clean()
    else:
        print("Invalid step:", step)
        print(USAGE)
        sys.exit(1)

# -----------------------------------
# End of Run pipeline steps
