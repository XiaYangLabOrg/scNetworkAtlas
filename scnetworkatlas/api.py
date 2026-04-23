"""Public, in-process API for the scNetworkAtlas pipeline.

This module provides the importable equivalent of two ``run_pipeline.py``
subcommands so callers can run the algorithms without spawning a
subprocess and without writing a generated ``config.py`` to disk.

The two functions cover the steps needed by the scGRNdb
``module_pathway`` flow:

    run_gene_membership(config) — module detection on a single network
                                  (≡ ``run_pipeline.py gene_membership``)
    run_enrichment(config)      — pathway enrichment for one
                                  (module, pathway) pair
                                  (≡ ``run_pipeline.py enrichment``)

Design notes
------------

* Per-pair, not per-batch. The legacy bash shims
  (``shell_scripts/run_genemembership_local.sh``,
  ``shell_scripts/run_enrichment_local.sh``) loop over ``network_file``
  and ``module_file`` listings externally, then invoke the per-pair
  Python script in each iteration. These wrappers expose only the
  per-pair Python call; the caller is expected to do its own loop. This
  matches what the algorithm modules already accept and avoids
  re-implementing the bash listing logic in Python.

* No file-based config. ``run_pipeline.py`` reads a generated
  ``config.py`` (Python source). These wrappers accept a Python ``dict``
  argument so callers do not have to write executable code to disk.

* Algorithm parity. The wrappers call the same functions in
  ``scnetworkatlas/python_files/`` that the bash shims call today
  (``ModuleBasedDimensionalityReduction.GeneMemByMod_optimizer`` and
  the body of ``enrichment.py`` factored out into a callable). No
  algorithmic change.

* Status: SKELETON. Both functions raise ``NotImplementedError``. The
  implementation will land in a follow-up commit on this branch
  (``secure-backend``) once the API shape is reviewed.
"""

from __future__ import annotations

from typing import Any, Dict


# ---------------------------------------------------------------------------
# Type aliases — kept simple here; if the package ever grows we can promote
# these to TypedDicts or pydantic models. For now ``dict[str, Any]`` is the
# narrowest contract that does not lock the wrapper into a particular
# validation library.
# ---------------------------------------------------------------------------

GeneMembershipConfig = Dict[str, Any]
GeneMembershipResult = Dict[str, Any]
EnrichmentConfig = Dict[str, Any]
EnrichmentResult = Dict[str, Any]


def run_gene_membership(config: GeneMembershipConfig) -> GeneMembershipResult:
    """Run module detection (Leiden / CPM) on a single network.

    Equivalent to ``python run_pipeline.py gene_membership`` followed by
    one iteration of ``shell_scripts/run_genemembership_local.sh`` for
    a single network file. Internally calls
    ``python_files.ModuleBasedDimensionalityReduction.GeneMemByMod_optimizer``.

    Parameters
    ----------
    config : dict
        Required keys (all values must be valid file paths or
        plain Python scalars; no shell metacharacters are interpreted):

        ``network_file`` : str
            Absolute path to the input network file. Must be readable
            and parseable by pandas as a TSV (``.txt``) or CSV
            (``.csv`` / ``.csv.gz``) edge list with columns
            (source, target[, weight]).
        ``out_file`` : str
            Absolute path where the gene-membership table will be
            written. Parent directory must exist. Output format is
            inferred from the extension (``.csv`` → comma-separated;
            anything else → tab-separated). Two columns are written:
            ``genes`` and ``module``.
        ``min_module_size`` : int
            Minimum module size (genes) — modules below this are
            dropped. Typical: 10.
        ``max_module_size`` : int
            Maximum allowed module size (genes). The Leiden resolution
            is selected so the largest module fits this cap. Typical:
            300.

    Returns
    -------
    dict
        ``{"out_file": <out_file>, "n_modules": <int>, "n_genes": <int>}``
        — the absolute output path and basic module/gene counts after
        size-filtering. Useful for callers that want to surface a
        warning when zero modules survive.

    Raises
    ------
    FileNotFoundError
        ``network_file`` does not exist or its parent directory is
        unreadable.
    ValueError
        ``min_module_size`` >= ``max_module_size``, or the network
        edge list cannot be parsed.
    NotImplementedError
        Always raised in this skeleton commit. Implementation lands in
        a follow-up.

    Notes
    -----
    Determinism: the underlying ``leidenalg.Optimiser`` is seeded with
    ``rng_seed=0`` (see ``GeneMemByMod_optimizer``); given identical
    inputs and library versions, output is reproducible across runs.

    No conda environment switching. Caller must already be in an env
    with ``numpy``, ``pandas``, ``networkx``, ``leidenalg``, and
    ``python-igraph`` installed.
    """
    raise NotImplementedError(
        "scnetworkatlas.api.run_gene_membership is a skeleton; "
        "implementation lands in a follow-up commit on the "
        "secure-backend branch."
    )


def run_enrichment(config: EnrichmentConfig) -> EnrichmentResult:
    """Run pathway enrichment for one (module-file, pathway-file) pair.

    Equivalent to ``python run_pipeline.py enrichment`` followed by
    one iteration of ``shell_scripts/run_enrichment_local.sh`` for
    a single (module, pathway) pair. Internally calls the algorithm
    body of ``python_files/enrichment.py`` after that script is
    refactored to expose a callable ``run(...)`` function (the
    ``__main__`` block keeps a thin argparse shim that calls the same
    ``run(...)`` for backwards compatibility with any external caller).

    The hypergeometric statistic, FDR correction, and risk-ratio
    calculation are unchanged from the script form.

    Parameters
    ----------
    config : dict
        Required keys:

        ``module_file`` : str
            Absolute path to a gene-membership file (the output of
            ``run_gene_membership``). Two columns: module ID and gene
            symbol; column names are configurable below.
        ``pathway_file`` : str
            Absolute path to a pathway database file (TSV). Two
            columns: pathway ID and gene symbol; column names are
            configurable below.
        ``out_file`` : str
            Absolute path where the per-pair enrichment table is
            written. Parent directory must exist. Format inferred from
            extension (``.txt`` / ``.tsv`` → tab-separated; otherwise
            comma-separated). Columns: ``MODULE``, ``PATHWAY``, ``P``,
            ``FDR``, ``risk_ratio``, ``overlap``, ``pathway_size``,
            ``module_size``, ``overlap_genes``.

        Optional keys (defaults match ``enrichment.py`` argparse):

        ``module_name_col`` : str = ``"module"``
        ``module_gene_col`` : str = ``"genes"``
        ``pathway_name_col`` : str = ``"MODULE"``
        ``pathway_gene_col`` : str = ``"GENE"``
        ``min_overlap`` : int = ``10``
            Minimum (module ∩ pathway) size required to test
            enrichment. The script applies a softer floor for small
            modules (``min(min_overlap, 0.2 * module_size)``); that
            behavior is preserved.
        ``pathway_size_min`` : int = ``10``
        ``pathway_size_max`` : int = ``500``
        ``n_background_genes`` : int = ``20000``
            Universe size for the hypergeometric test.

    Returns
    -------
    dict
        ``{"out_file": <out_file>, "n_pathways_enriched": <int>}`` —
        the absolute output path and the number of (module, pathway)
        rows in the output (post-filter, pre-FDR-correction). Zero is
        a valid result; the runner uses this to log a warning per
        pathway DB.

    Raises
    ------
    FileNotFoundError
        ``module_file`` or ``pathway_file`` does not exist.
    KeyError
        A configured column name is missing from the input file.
    NotImplementedError
        Always raised in this skeleton commit. Implementation lands in
        a follow-up.

    Notes
    -----
    The msigdb code path (``--use_msigdb`` in ``enrichment.py``) is
    intentionally **not** exposed by this wrapper. It is the only code
    path that imports ``decoupler``, and consumers of this API are
    expected to provide their own pathway files. Add a separate
    ``run_enrichment_decoupler(...)`` function later if the msigdb
    path is needed.

    No conda environment switching. Caller must already be in an env
    with ``pandas``, ``numpy``, ``polars``, and ``scipy`` installed.
    """
    raise NotImplementedError(
        "scnetworkatlas.api.run_enrichment is a skeleton; "
        "implementation lands in a follow-up commit on the "
        "secure-backend branch."
    )
