import decoupler as dc
import pandas as pd
import numpy as np
import os
import polars as pl
import scipy.stats as stats
import sys
import argparse

#argparse
parser = argparse.ArgumentParser(description="enrichment")
parser.add_argument('module_file', type=str, help='path to module file')
parser.add_argument('--module_name_col', type=str, default='module', help='module column name in module file')
parser.add_argument('--module_gene_col', type=str, default='genes', help='gene column name in module file')
# parser.add_argument('--mod_delim', type=str, default="\t", help="delimiter for module file. Default \"\\t\"")
parser.add_argument('--out_file', type=str, help='output file')
# if using msigdb pathways through decoupler api
parser.add_argument('--use_msigdb', action="store_true", default=False, help="whether to use decoupler api to access pathways. Default False")
parser.add_argument('--pathway', type=str, default="go_biological_process", help='which pathway to use, see decoupler msigDB collection')
# if using local pathway files
parser.add_argument('--pathway_file', type=str, default=None, help="pathway file")
parser.add_argument('--pathway_name_col', type=str, default='MODULE', help='pathway column name in pathway file')
parser.add_argument('--pathway_gene_col', type=str, default='GENE', help='gene column in pathway file')
# parser.add_argument('--pathway_delim', type=str, default="\t", help="delimiter for pathway file. Default \"\\t\"")
# additional arguments
parser.add_argument('--min_overlap', type=int, default=5, help="miniumum pathway-module overlap needed for enrichment analysis")
parser.add_argument('--pathway_size_min', type=int, default=10, help='pathway size limit')
parser.add_argument('--pathway_size_max', type=int, default=500, help='pathway size limit')
parser.add_argument('--n_background_genes', type=int, default=20000, help="number of background genes for hypergeometric test. Default 20000")

args = parser.parse_args()

out_file = args.out_file
pathway = args.pathway
pathway_file = args.pathway_file
pathway_name_col = args.pathway_name_col
pathway_gene_col = args.pathway_gene_col
pathway_size_min = args.pathway_size_min
pathway_size_max = args.pathway_size_max
min_overlap = args.min_overlap

module_file = args.module_file
module_name_col = args.module_name_col
module_gene_col = args.module_gene_col

n_background_genes = args.n_background_genes
# read pathway file
if pathway_file:
    if pathway_file.endswith(".csv"):
        separator = ","
    else:
        separator = "\t"
    db = pl.read_csv(pathway_file, separator="\t").rename({pathway_name_col:'PATHWAY', pathway_gene_col:"GENE"})
    db = db.group_by("PATHWAY").agg(pl.col("GENE").unique())
elif args.use_msigdb and args.pathway:
    msigdb = dc.get_resource("MSigDB")
    msigdb = pl.from_pandas(msigdb)
    msigdb = msigdb.filter(pl.col('collection') == pathway)
    msigdb = msigdb.unique(("geneset", "genesymbol"))
    db = msigdb.group_by("geneset").agg(pl.col("genesymbol")).rename({"genesymbol": "GENE", "geneset":"PATHWAY"})
elif args.use_msig and not args.pathway:
    raise argparse.ArgumentError("pathway should both be set if --use_msigdb is set")
else:
    raise argparse.ArgumentError("--use_msig and --pathway OR --pathway_file must be set")

# read module file
if module_file.endswith(".csv"):
    separator = ","
else:
    separator = "\t"
module_df = pl.read_csv(module_file, separator=separator).cast(pl.String).rename({module_name_col:'MODULE', module_gene_col:'GENE'})
score_di = {}

# calculate overlap
for mod, mod_genes in module_df.group_by("MODULE").agg(pl.col("GENE")).iter_rows(): # module -> genes
    mod_PE = db.with_columns(
        pl.col("GENE").map_elements(lambda x: list(set(mod_genes) & set(x)), return_dtype=pl.List(pl.String)).alias("overlap_genes"),
        pl.col("GENE").map_elements(lambda x: len(x), return_dtype=pl.Int32).alias("pathway_size"),
        pl.Series("module_size", [len(mod_genes)]*len(db))
    )
    mod_PE = mod_PE.with_columns(pl.col("overlap_genes").map_elements(lambda x: len(x), return_dtype=pl.Int32).alias("overlap"))

    # filter pathways
    mod_PE = mod_PE.filter(pl.col("pathway_size") >= pathway_size_min,
                                   pl.col("pathway_size") <= pathway_size_max,
                                   pl.col("overlap") > min(min_overlap, 0.2*(len(mod_genes)))) # use smaller overlap threshold for small modules
    mod_PE = mod_PE.with_columns(
        P = stats.hypergeom.sf(mod_PE["overlap"]-1, n_background_genes, mod_PE["pathway_size"], len(mod_genes)), 
        risk_ratio = (mod_PE["overlap"] *n_background_genes / (len(mod_genes) * mod_PE["pathway_size"])),
    )
    mod_PE = mod_PE.with_columns(
        FDR = stats.false_discovery_control(mod_PE["P"] , axis=0)
    )
    if len(mod_PE) == 0:
        continue
    mod_PE = mod_PE.sort("FDR")
    score_di[mod] = mod_PE

final_PE = pl.DataFrame()
if len(score_di.keys())>0:
    for mod, x in score_di.items():
        print(mod, " | ", x['FDR'].min())
        x = x.with_columns(pl.Series("MODULE", [mod]*len(x)),
                        pl.col("overlap_genes").list.join(separator=";"))
        final_PE = pl.concat([final_PE, x], how='vertical')
    print(final_PE.head())
    final_PE = final_PE[["MODULE", "PATHWAY", "P", "FDR" , "risk_ratio", "overlap", "pathway_size", "module_size", "overlap_genes"]].sort("FDR")

    if out_file.split(".")[-1] == "txt" or out_file.split(".")[-1] == "tsv":
        final_PE.write_csv(out_file, separator="\t")
    else:
        final_PE.write_csv(out_file)
else:
    if pathway_file:
        print(f'no pathway enrichment result from {pathway_file} for {module_file}')
    else:
        print(f'no pathway enrichment result from {pathway} for {module_file}')