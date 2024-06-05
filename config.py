# User-configurable inputs
# -----------------------------------

# Configuration settings (these settings are read into submission scripts)
scing_scripts_base_dir = "/u/home/k/kylepu/scratch/env_vars_playground" # directory storing cloned scGRNdb repo
base_dir = "/u/project/xyang123/shared/reference/single_cell_databases/"
scing_config = {
    'conda_init_script': '~/project-xyang123/miniconda3/etc/profile.d/conda.sh',

    'cell_mapping': {
        'base_dir': base_dir,
        'mapping_file': f"{base_dir}all_celltypes/mouse_TabMurSenis_FACS.cleaned.tsv",
        'adata_dir': f"{base_dir}mouse/Tab_Mur_Senis_FACS/adatas/",
        'celltype_column': "celltypes"
    },
    'pseudobulking': {
        'filetype': 'npz' # npz or h5ad
    },
    'pseudobulking_v2': { # note: only use h5ad in pseudobulking v2
        'celltype_col': 'Cell.Type',
        'sample_col': 'SampleID',
        'recluster': True,
    },
    'build_grn': {
        'num_networks': 100
    },
    'merge_networks': {
#        'config': None # no inputs needed for merge networks step
        'consensus': [0.2, 0.6, 0.8]
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
