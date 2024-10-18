### HOW TO USE ###

# PREREQUISITES: scing miniconda environment exists (build from SCING repo) #

# BEFORE YOU BEGIN: fill out config.py with the appropriate details for your project #

# Run all steps in the directory of your project directory #

### END OF HOW TO USE ###

USAGE = "python run_pipeline.py <setup|cell_mapping|pseudobulking_v2|pseudobulking|build_grn|merge_networks|gene_membership|enrichment|enrichment_decoupler|clean>"
env_activated = False

import os
import sys
import argparse

# Fetch config from current dir
from config import base_dir, scing_config

def activate_conda_environment(conda_env):
    try:
        combined_command = f"conda activate {conda_env}"
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
def setup():
    for tmp in ["submission_scripts", "python_files", "r_files", "shell_scripts"]:
        os.makedirs(f"temp/{tmp}", exist_ok=True)
    print("copying scripts")
    cmds = scing_config['copy_cmds']
    for cmd in cmds:
        os.system(cmd)
    print("done copying scripts")
    
def confirm_conda_activated(conda_env='scing'):
    if not env_activated:
        confirmation = input(f"Is your scing {conda_env} environment activated? (y/n): ").lower()
        if confirmation == 'y':
            return True
        elif confirmation == 'n':
            print("Please activate your conda environment and run the script again.")
            activate_conda_environment(conda_env)
            sys.exit()
        else:
            print("Invalid input. Please enter 'y' or 'n'.")
            return False
    else:
        print(f'{conda_env} environment is activated')
        return True

def join_list_arg(l):
    # combine list argument to string
    if isinstance(l, list):
        return ','.join(map(str,l))
    else:
        return str(l)

def cell_mapping():
    confirm_conda_activated()

    config = scing_config['cell_mapping']
    cmd = f"bash temp/submission_scripts/submit_run_cellmapping.sh {base_dir} {config['mapping_file']} {config['adata_dir']} {config['celltype_column']}"
    os.system(cmd)

def pseudobulking_v2():
    confirm_conda_activated()

    config = scing_config['pseudobulking_v2']
    cmd = f"bash temp/submission_scripts/submit_run_supercells_v2.sh {config['tissue_dir']} {config['supercell_dir']} {config['celltype_col']} {config['sample_col']} {config['tissue_celltype_file']}"
    os.system(cmd)

def pseudobulking():
    confirm_conda_activated()

    config = scing_config['pseudobulking']
    cmd = f"bash temp/submission_scripts/submit_run_supercells.sh {config['tissue_dir']} {config['supercell_dir']} {config['filetype']} {config['celltype_col']} {config['tissue_celltype_file']}"
    os.system(cmd)
    
def build_grn():
    confirm_conda_activated()

    config = scing_config['build_grn']
    cmd = f"bash temp/submission_scripts/submit_run_buildgrn.sh {config['num_networks']} {config['supercell_dir']} {config['supercell_file']} {config['out_dir']} {config['ncore']} {config['mem_per_core']}"
    os.system(cmd)

def merge_networks():
    confirm_conda_activated()

    config = scing_config['merge_networks']
    consensus_str = join_list_arg(config['consensus']) # pass as string to submission script
    cmd = f"bash temp/submission_scripts/submit_run_merge.sh {config['supercell_dir']} {config['supercell_file']} {config['intermediate_dir']} {consensus_str} {config['out_dir']} {config['ncore']} {config['mem_per_core']}"
    os.system(cmd)

def gene_membership():
    confirm_conda_activated()
    
    config = scing_config['gene_membership']
    cmd = f"bash temp/submission_scripts/submit_run_genemembership.sh {config['network_dir']} {config['network_file']} {config['out_dir']} {config['min_module_size']} {config['max_module_size']} {config['network_ext']} {config['submit_command']}"
    os.system(cmd)

def enrichment():
    confirm_conda_activated(conda_env='decoupler')
    
    config = scing_config['enrichment']
    pathway_file_str = join_list_arg(config['pathway_file'])
    pathway_db_str = join_list_arg(config['pathway_db'])
    # must be same number of pathway files as pathway db names
    assert len(pathway_file_str.split(',')) == len(pathway_db_str.split(',')), 'Error. Must have same number of pathway files as pathway dbs'
    cmd = f"bash temp/submission_scripts/submit_run_enrichment.sh {config['module_dir']} {config['module_file']} {config['module_name_col']} {config['module_gene_col']} {pathway_file_str} {pathway_db_str} {config['pathway_name_col']} {config['pathway_gene_col']} {config['min_overlap']} {config['pathway_size_min']} {config['pathway_size_max']} {config['out_dir']} {config['submit_command']}"
    os.system(cmd)

def enrichment_decoupler():
    confirm_conda_activated(conda_env='decoupler')
    
    config = scing_config['enrichment_decoupler']
    pathway_str = join_list_arg(config['pathway'])
    cmd = f"bash temp/submission_scripts/submit_run_enrichment_decoupler.sh {config['module_dir']} {config['module_file']} {config['module_name_col']} {config['module_gene_col']} {pathway_str} {config['min_overlap']} {config['pathway_size_min']} {config['pathway_size_max']} {config['out_dir']} {config['submit_command']}"
    os.system(cmd)


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
    parser = argparse.ArgumentParser(description="Running SCING Pipeline")
    parser.add_argument("step", type=str, help="step to run")
    parser.add_argument("-y", "--env_activated",  "--y", dest="env_activated", action='store_true', help="indicates conda environment is activated, skips environment check", default=False)
    args = parser.parse_args()
    if len(sys.argv) < 2:
        print(USAGE)
        sys.exit(1)
    step = args.step
    env_activated = args.env_activated
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
    elif step == "enrichment":
        enrichment()
    elif step == "enrichment_decoupler":
        enrichment_decoupler()
    elif step == "clean":
        clean()
    else:
        print("Invalid step:", step)
        print(USAGE)
        sys.exit(1)

# -----------------------------------
# End of Run pipeline steps
