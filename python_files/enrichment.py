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
parser.add_argument('read_file', type=str, help='path to data h5ad file')
parser.add_argument('out_file', type=str, help='output directory')
parser.add_argument('pathway', type=str, default="go_biological_process", help='which pathway to use, see decoupler msigDB collection')
parser.add_argument('pathway_size_min', type=int, default=10, help='pathway size limit')
parser.add_argument('pathway_size_max', type=int, default=500, help='pathway size limit')
parser.add_argument('pathway_col', type=str, default='pathway', help='pathway column name')
parser.add_argument('module_col', type=str, default='module', help='module column name')

args = parser.parse_args()

read_file = args.read_file
write_file = args.out_file
pathway = args.pathway
pathway_size_min = args.pathway_size_min
pathway_size_max = args.pathway_size_max
pathway_col = args.pathway_col
module_col = args.module_col

msigdb = dc.get_resource("MSigDB")
msigdb = pl.from_pandas(msigdb)
msigdb = msigdb.filter(pl.col('collection') == pathway)
msigdb = msigdb.unique(("geneset", "genesymbol"))

db = msigdb.group_by("geneset").agg(pl.col("genesymbol")).rename({"genesymbol": "db_genesymbol"})

if read_file.endswith(".csv"):
    separator = ","
else:
    separator = "\t"
df = pl.read_csv(read_file, separator=separator)
df = df.with_columns(pl.col("genes").str.to_uppercase()) # adjustment genename

score_di = {}

for cluster, geneset in df.group_by("module").agg(pl.col("genes")).iter_rows(): # cluster -> module
    cluster_db = db.with_columns(
        pl.col("db_genesymbol").map_elements(lambda x: len(set(geneset) & set(x)), return_dtype=pl.Int32).alias("overlap"),
        pl.col("db_genesymbol").map_elements(lambda x: len(x), return_dtype=pl.Int32).alias("pathway_size"),
        pl.Series("module_size", [len(geneset)]*len(db))
    )
    
    cluster_db = cluster_db.filter(pl.col("pathway_size") >= pathway_size_min,
                                   pl.col("pathway_size") <= pathway_size_max,
                                   pl.col("overlap") > 0)
    cluster_db = cluster_db.with_columns(
        p = stats.hypergeom.sf(cluster_db["overlap"]-1, 20000, cluster_db["pathway_size"], len(geneset)), 
        risk_ratio = (cluster_db["overlap"] *20000 / (len(geneset) * cluster_db["pathway_size"])),
    )
    cluster_db = cluster_db.with_columns(
        adj_p = stats.false_discovery_control(cluster_db["p"] , axis=0)
    )
    cluster_db = cluster_db.filter(pl.col("adj_p") < 0.01)
    if len(cluster_db) == 0:
        continue
    cluster_db = cluster_db.sort("adj_p")
    score_di[cluster] = cluster_db

new = pl.DataFrame()
for cluster, x in score_di.items():
    print(cluster, " | ", x['adj_p'].min())
    x = x.with_columns(pl.Series("cluster", [cluster]*len(x)))
    new = pl.concat([new, x], how='vertical')

new = new[["cluster", "p", "adj_p" , "risk_ratio", "overlap", "pathway_size", "module_size", "geneset"]].sort("adj_p")
new = new.rename({
    "geneset": pathway_col,
    "cluster": module_col,})
new.sort("p").head(10)
if write_file.split(".")[-1] == "txt" or write_file.split(".")[-1] == "tsv":
    new.write_csv(write_file, separator="\t")
else:
    new.write_csv(write_file)
