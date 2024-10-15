# User-configurable inputs
# -----------------------------------

# Configuration settings (these settings are read into submission scripts)
main_branch_path = "/u/scratch/m/mikechen/scNetworkAtlas/"
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
        'supercell_file': 'supercells.txt',
        'out_dir': "saved_networks/intermediate_networks",
        'ncore': 1,
        'mem_per_core': 16,
    },
    'merge_networks': {
        'supercell_dir': 'supercells',
        'supercell_file': 'supercell_file.txt',
        'intermediate_dir': "saved_networks/intermediate_networks",
        'consensus': [0.5],
        'out_dir': 'saved_networks/final_networks',
        'ncore': 12,
        'mem_per_core': 4 
    },
    'gene_membership': {
        'network_dir': 'saved_networks/final_networks',
        'network_file': 'final_networks.txt',
        #'q1_module_sizes': ('20' '35' '50')
        'out_dir': 'gene_memberships',
        "min_module_size": 10,
        "max_module_size": 300,
        
    },
    'enrichment':{
        'modules_dir': 'gene_memberships',
        'module_file': 'gene_memberships.txt',
        'out_dir': 'enrichment',
        "pathway": "go_biological_process",
        "pathway_size_min": 10,
        "pathway_size_max": 500,
        "pathway_col": "pathway",
        "module_col": "module",
        # choice:
        # chemical_and_genetic_perturbations immunesigdb mirna_targets_mirdb go_molecular_function tf_targets_gtrf tf_targets_legacy oncogenic_signatures cell_type_signatures vaccine_response go_biological_process cancer_gene_neighborhoods 
        # cancer_modules go_cellular_component wikipathways reactome_pathways hallmark mirna_targets_legacy biocarta_pathways positional human_phenotype_ontology pid_pathways kegg_pathways
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
