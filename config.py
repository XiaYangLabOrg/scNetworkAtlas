# User-configurable inputs
# -----------------------------------

# Configuration settings (these settings are read into submission scripts)
main_branch_path = "/u/home/s/skikuchi/project-xyang123/scNetworkAtlas/"
#cloned scGRNdb repo
base_dir = "/u/project/xyang123/shared/reference/single_cell_databases/"

scing_config = {
    'conda_init_script': '~/project-xyang123/miniconda3/etc/profile.d/conda.sh',
    'copy_cmds': 
            [f"cp -r {main_branch_path}shell_scripts/* temp/shell_scripts/",
             f"cp -r {main_branch_path}submission_scripts/* temp/submission_scripts/",
             f"cp -r {main_branch_path}python_files/* temp/python_files/"],
        
    'cell_mapping': {
        'base_dir': base_dir,
        'mapping_file': f"{base_dir}all_celltypes/mouse_TabMurSenis_FACS.cleaned.tsv",
        'adata_dir': f"{base_dir}mouse/Tab_Mur_Senis_FACS/adatas/",
        'celltype_column': "celltypes"
    },
    'pseudobulking': {
        'tissue_dir': "tissue_adata",
        'supercell_dir': "supercells",
        'filetype': "h5ad", # npz or h5ad
        'celltype_col': "celltypes",
        'tissue_celltype_file': 'tissue_celltype_file.txt',
    },
    'pseudobulking_v2': { # note: only use h5ad in pseudobulking v2
        'tissue_dir': 'tissue_adata',
        'supercell_dir':'supercells',
        'tissue_celltype_file':'tissue_celltype_file.txt',        
        'celltype_col': 'Cell.Type',
        'sample_col': 'SampleID',
        'recluster': True, # TODO: not used yet
    },
    'build_grn': {
        'num_networks': 100,
        'supercell_dir': 'supercells',
        'supercell_file': 'supercell_file.txt',
        'out_dir': "saved_networks/intermediate_networks",
        'ncore': 1,
        'mem_per_core': 16,
    },
    'merge_networks': {
        'supercell_dir': 'supercells',
        'supercell_file': 'supercell_file.txt',
        'intermediate_dir': 'saved_networks/intermediate_networks',
        'consensus': [0.2, 0.5, 0.8], # must be list
        'out_dir': 'saved_networks/final_networks',
        'ncore': 12,
        'mem_per_core': 4 
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
