# User-configurable inputs
# -----------------------------------

# Configuration settings (these settings are read into submission scripts)
main_branch_path = "/u/scratch/m/mikechen/scNetworkAtlas/"
#cloned scGRNdb repo
base_dir = "/u/project/xyang123/shared/reference/single_cell_databases/"

scing_config = {
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
        'network_ext': 'csv.gz',
        'out_dir': 'gene_memberships',
        "min_module_size": 10,
        "max_module_size": 300,
        "submit_command": "bash", # qsub or bash
        
    },
    'enrichment':{
        'module_dir': 'gene_memberships',
        'module_file': 'gene_membership_files',
        'module_name_col': 'module',
        'module_gene_col': 'genes',
        'pathway_file': ['/u/project/xyang123/shared/reference/pathways/human/GO_BP_Hs.txt','/u/project/xyang123/shared/reference/pathways/human/DisGeNET.txt'],
        'pathway_db': ['GO_BP','DisGeNET'],
        'pathway_name_col': 'MODULE',
        'pathway_gene_col': 'GENE',
        'min_overlap': 10,
        "pathway_size_min": 10,
        "pathway_size_max": 500,
        'out_dir': 'enrichment',
        'submit_command': 'bash', # qsub or bash


    },
    'enrichment_decoupler':{
        'module_dir': 'gene_memberships',
        'module_file': 'gene_membership_files',
        'module_name_col': 'module',
        'module_gene_col': 'genes',
        "pathway": ["go_biological_process","go_molecular_function"],
        'min_overlap': 10,
        "pathway_size_min": 10,
        "pathway_size_max": 500,
        'out_dir': 'enrichment_decoupler',
        'submit_command': 'bash', # qsub or bash
        # choice:
        # chemical_and_genetic_perturbations immunesigdb mirna_targets_mirdb go_molecular_function tf_targets_gtrf tf_targets_legacy oncogenic_signatures cell_type_signatures vaccine_response go_biological_process cancer_gene_neighborhoods 
        # cancer_modules go_cellular_component wikipathways reactome_pathways hallmark mirna_targets_legacy biocarta_pathways positional human_phenotype_ontology pid_pathways kegg_pathways
    },
    
}
    
# -----------------------------------
# End of User-configurable inputs
