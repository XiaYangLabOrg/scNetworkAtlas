# User-configurable inputs
# -----------------------------------

# Configuration settings (these settings are read into submission scripts)
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
        'tissue_dir': "tissue_adata",
        'supercell_dir': "supercells",
        'filetype': "h5ad", # npz or h5ad
        'celltype_col': "celltypes",
        'tissue_celltype_file': 'tissue_celltype_file.txt',
    },
    'pseudobulking_v2': { # note: only use h5ad in pseudobulking v2
        'tissue_dir': 'tissue_adata',
        'supercell_dir':'supercells',
        'celltype_col': 'celltypes',
        'sample_col': 'sample',
        'tissue_celltype_file':'tissue_celltype_file.txt',
    },
    'build_grn': {
        'num_networks': 100,
        'supercell_dir':"supercells",
        'supercell_file':"supercells/supercells.txt",
        'out_dir':"saved_networks/intermediate_data",
    },
    'merge_networks': {
        'supercell_dir': "supercells",
        'supercell_file': "supercells/supercells.txt",
        'intermediate_dir': "saved_networks/intermediate_data",
        'consensus_thresholds': 0.5, # need to figure out how feed a list of values to bash
        'out_dir': "saved_networks/final_edges",
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
