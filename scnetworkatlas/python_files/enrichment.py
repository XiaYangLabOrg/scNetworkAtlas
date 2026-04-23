"""Pathway enrichment for one (module-file, pathway-file) pair.

Public callable: ``run(...)``. The ``__main__`` block below is a thin
argparse shim that delegates to ``run(...)``; ``python enrichment.py
...`` continues to behave exactly as before for legacy callers
(including ``shell_scripts/run_enrichment_local.sh``).

Algorithm: per-module hypergeometric test of (module ∩ pathway) over a
fixed background, then per-module FDR (BH) correction. Sort by FDR.
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional

import numpy as np  # noqa: F401  (kept for parity with original imports)
import os  # noqa: F401  (kept for parity with original imports)
import pandas as pd
import polars as pl
import scipy.stats as stats


def run(
    module_file: str,
    out_file: str,
    *,
    pathway_file: Optional[str] = None,
    use_msigdb: bool = False,
    pathway: str = "go_biological_process",
    module_name_col: str = "module",
    module_gene_col: str = "genes",
    pathway_name_col: str = "MODULE",
    pathway_gene_col: str = "GENE",
    min_overlap: int = 5,
    pathway_size_min: int = 10,
    pathway_size_max: int = 500,
    n_background_genes: int = 20000,
) -> int:
    """Compute enrichment for one (module-file, pathway) pair.

    Either ``pathway_file`` (a local TSV pathway database) **or**
    ``use_msigdb=True`` with ``pathway`` (an msigdb collection name
    fetched via ``decoupler.get_resource``) must be provided.

    Returns the number of (module, pathway) rows written to
    ``out_file``. Zero is a valid result and produces no output file
    (parity with the legacy script).
    """
    if pathway_file:
        if pathway_file.endswith(".csv"):
            separator = ","  # noqa: F841  (parity: original script also computes but does not use)
        else:
            separator = "\t"  # noqa: F841
        db = pl.read_csv(pathway_file, separator="\t").rename(
            {pathway_name_col: "PATHWAY", pathway_gene_col: "GENE"}
        )
        db = db.group_by("PATHWAY").agg(pl.col("GENE").unique())
    elif use_msigdb and pathway:
        import decoupler as dc

        msigdb = dc.get_resource("MSigDB")
        msigdb = pl.from_pandas(msigdb)
        msigdb = msigdb.filter(pl.col("collection") == pathway)
        msigdb = msigdb.unique(("geneset", "genesymbol"))
        db = msigdb.group_by("geneset").agg(pl.col("genesymbol")).rename(
            {"genesymbol": "GENE", "geneset": "PATHWAY"}
        )
    else:
        raise ValueError(
            "Either pathway_file or (use_msigdb=True with pathway=<collection>) "
            "must be provided."
        )

    if module_file.endswith(".csv"):
        separator = ","
    else:
        separator = "\t"
    module_df = pl.read_csv(module_file, separator=separator).cast(pl.String).rename(
        {module_name_col: "MODULE", module_gene_col: "GENE"}
    )
    score_di = {}

    for mod, mod_genes in module_df.group_by("MODULE").agg(pl.col("GENE")).iter_rows():
        mod_PE = db.with_columns(
            pl.col("GENE")
            .map_elements(
                lambda x: list(set(mod_genes) & set(x)),
                return_dtype=pl.List(pl.String),
            )
            .alias("overlap_genes"),
            pl.col("GENE")
            .map_elements(lambda x: len(x), return_dtype=pl.Int32)
            .alias("pathway_size"),
            pl.Series("module_size", [len(mod_genes)] * len(db)),
        )
        mod_PE = mod_PE.with_columns(
            pl.col("overlap_genes")
            .map_elements(lambda x: len(x), return_dtype=pl.Int32)
            .alias("overlap")
        )

        mod_PE = mod_PE.filter(
            pl.col("pathway_size") >= pathway_size_min,
            pl.col("pathway_size") <= pathway_size_max,
            pl.col("overlap") > min(min_overlap, 0.2 * (len(mod_genes))),
        )
        mod_PE = mod_PE.with_columns(
            P=stats.hypergeom.sf(
                mod_PE["overlap"] - 1,
                n_background_genes,
                mod_PE["pathway_size"],
                len(mod_genes),
            ),
            risk_ratio=(
                mod_PE["overlap"]
                * n_background_genes
                / (len(mod_genes) * mod_PE["pathway_size"])
            ),
        )
        mod_PE = mod_PE.with_columns(
            FDR=stats.false_discovery_control(mod_PE["P"], axis=0)
        )
        if len(mod_PE) == 0:
            continue
        mod_PE = mod_PE.sort("FDR")
        score_di[mod] = mod_PE

    final_PE = pl.DataFrame()
    if len(score_di.keys()) > 0:
        for mod, x in score_di.items():
            x = x.with_columns(
                pl.Series("MODULE", [mod] * len(x)),
                pl.col("overlap_genes").list.join(separator=";"),
            )
            final_PE = pl.concat([final_PE, x], how="vertical")
        final_PE = final_PE[
            [
                "MODULE",
                "PATHWAY",
                "P",
                "FDR",
                "risk_ratio",
                "overlap",
                "pathway_size",
                "module_size",
                "overlap_genes",
            ]
        ]
        final_PE = pd.DataFrame(final_PE)
        final_PE.columns = [
            "MODULE",
            "PATHWAY",
            "P",
            "FDR",
            "risk_ratio",
            "overlap",
            "pathway_size",
            "module_size",
            "overlap_genes",
        ]
        final_PE = final_PE.sort_values("FDR")

        if out_file.split(".")[-1] in ("txt", "tsv"):
            final_PE.to_csv(out_file, sep="\t", index=False)
        else:
            final_PE.to_csv(out_file, index=False)
        return len(final_PE)

    if pathway_file:
        print(
            f"no pathway enrichment result from {pathway_file} for {module_file}"
        )
    else:
        print(f"no pathway enrichment result from {pathway} for {module_file}")
    return 0


def _build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="enrichment")
    parser.add_argument("module_file", type=str, help="path to module file")
    parser.add_argument(
        "--module_name_col",
        type=str,
        default="module",
        help="module column name in module file",
    )
    parser.add_argument(
        "--module_gene_col",
        type=str,
        default="genes",
        help="gene column name in module file",
    )
    parser.add_argument("--out_file", type=str, help="output file")
    parser.add_argument(
        "--use_msigdb",
        action="store_true",
        default=False,
        help="whether to use decoupler api to access pathways. Default False",
    )
    parser.add_argument(
        "--pathway",
        type=str,
        default="go_biological_process",
        help="which pathway to use, see decoupler msigDB collection",
    )
    parser.add_argument(
        "--pathway_file", type=str, default=None, help="pathway file"
    )
    parser.add_argument(
        "--pathway_name_col",
        type=str,
        default="MODULE",
        help="pathway column name in pathway file",
    )
    parser.add_argument(
        "--pathway_gene_col",
        type=str,
        default="GENE",
        help="gene column in pathway file",
    )
    parser.add_argument(
        "--min_overlap",
        type=int,
        default=5,
        help="miniumum pathway-module overlap needed for enrichment analysis",
    )
    parser.add_argument(
        "--pathway_size_min", type=int, default=10, help="pathway size limit"
    )
    parser.add_argument(
        "--pathway_size_max", type=int, default=500, help="pathway size limit"
    )
    parser.add_argument(
        "--n_background_genes",
        type=int,
        default=20000,
        help="number of background genes for hypergeometric test. Default 20000",
    )
    return parser


if __name__ == "__main__":
    args = _build_argparser().parse_args()
    if not args.out_file:
        print("--out_file is required", file=sys.stderr)
        sys.exit(2)
    run(
        module_file=args.module_file,
        out_file=args.out_file,
        pathway_file=args.pathway_file,
        use_msigdb=args.use_msigdb,
        pathway=args.pathway,
        module_name_col=args.module_name_col,
        module_gene_col=args.module_gene_col,
        pathway_name_col=args.pathway_name_col,
        pathway_gene_col=args.pathway_gene_col,
        min_overlap=args.min_overlap,
        pathway_size_min=args.pathway_size_min,
        pathway_size_max=args.pathway_size_max,
        n_background_genes=args.n_background_genes,
    )
