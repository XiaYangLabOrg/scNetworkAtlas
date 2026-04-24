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
  ``scnetworkatlas.python_files`` that the bash shims call today
  (``ModuleBasedDimensionalityReduction.run`` and ``enrichment.run``).
  No algorithmic change.
"""

from __future__ import annotations

import os
from typing import Any, Dict


GeneMembershipConfig = Dict[str, Any]
GeneMembershipResult = Dict[str, Any]
EnrichmentConfig = Dict[str, Any]
EnrichmentResult = Dict[str, Any]


_GENE_MEMBERSHIP_REQUIRED = ("network_file", "out_file", "min_module_size", "max_module_size")
_ENRICHMENT_REQUIRED = ("module_file", "pathway_file", "out_file")
_ENRICHMENT_OPTIONAL_DEFAULTS: Dict[str, Any] = {
    "module_name_col": "module",
    "module_gene_col": "genes",
    "pathway_name_col": "MODULE",
    "pathway_gene_col": "GENE",
    "min_overlap": 10,
    "pathway_size_min": 10,
    "pathway_size_max": 500,
    "n_background_genes": 20000,
}


def _require_keys(config: Dict[str, Any], required: tuple, fn_name: str) -> None:
    missing = [k for k in required if k not in config]
    if missing:
        raise KeyError(
            f"{fn_name}: missing required config keys: {missing}. "
            f"Got keys: {sorted(config.keys())}"
        )


def _require_existing_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{label} not found: {path}")


def _require_writable_parent(path: str, label: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(parent):
        raise FileNotFoundError(
            f"{label} parent directory does not exist: {parent}"
        )


def run_gene_membership(config: GeneMembershipConfig) -> GeneMembershipResult:
    """Run module detection (Leiden / CPM) on a single network.

    Equivalent to ``python run_pipeline.py gene_membership`` followed by
    one iteration of ``shell_scripts/run_genemembership_local.sh`` for
    a single network file. Internally calls
    ``scnetworkatlas.python_files.ModuleBasedDimensionalityReduction.run``
    which in turn calls ``GeneMemByMod_optimizer``.

    Parameters
    ----------
    config : dict
        Required keys:

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
            ``genes`` and ``module`` (both lowercase). Downstream
            ``run_enrichment`` consumes this file and defaults its
            ``module_name_col`` / ``module_gene_col`` to match these
            exact lowercase headers; the pathway-file side of that
            function uses uppercase headers — see that function's
            docstring for why.
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
        size-filtering.

    Raises
    ------
    FileNotFoundError
        ``network_file`` does not exist or ``out_file``'s parent
        directory does not exist.
    KeyError
        Required config key missing.
    ValueError
        Network file extension not recognized, or module-size bounds
        invalid.

    Notes
    -----
    Determinism: the underlying ``leidenalg.Optimiser`` is seeded with
    ``rng_seed=0``; given identical inputs and library versions, output
    is reproducible across runs.

    No conda environment switching. Caller must already be in an env
    with ``numpy``, ``pandas``, ``networkx``, ``leidenalg``, and
    ``python-igraph`` installed.
    """
    _require_keys(config, _GENE_MEMBERSHIP_REQUIRED, "run_gene_membership")

    network_file = config["network_file"]
    out_file = config["out_file"]
    min_module_size = int(config["min_module_size"])
    max_module_size = int(config["max_module_size"])

    if min_module_size <= 0:
        raise ValueError(f"min_module_size must be > 0, got {min_module_size}")
    if max_module_size <= min_module_size:
        raise ValueError(
            f"max_module_size ({max_module_size}) must be > min_module_size "
            f"({min_module_size})"
        )

    _require_existing_file(network_file, "network_file")
    _require_writable_parent(out_file, "out_file")

    from scnetworkatlas.python_files.ModuleBasedDimensionalityReduction import run as _run

    gene_membership = _run(
        input=network_file,
        output=out_file,
        min_module_size=min_module_size,
        max_module_size=max_module_size,
    )

    return {
        "out_file": os.path.abspath(out_file),
        "n_modules": int(gene_membership["module"].nunique()),
        "n_genes": int(len(gene_membership)),
    }


def run_enrichment(config: EnrichmentConfig) -> EnrichmentResult:
    """Run pathway enrichment for one (module-file, pathway-file) pair.

    Equivalent to ``python run_pipeline.py enrichment`` followed by
    one iteration of ``shell_scripts/run_enrichment_local.sh`` for
    a single (module, pathway) pair. Internally calls
    ``scnetworkatlas.python_files.enrichment.run``.

    Parameters
    ----------
    config : dict
        Required keys:

        ``module_file`` : str
            Absolute path to a gene-membership file (the output of
            ``run_gene_membership``).
        ``pathway_file`` : str
            Absolute path to a pathway database file (TSV).
        ``out_file`` : str
            Absolute path where the per-pair enrichment table is
            written. Parent directory must exist.

        Optional keys (defaults match the legacy script's argparse,
        with one exception noted below):

        ``module_name_col`` : str = ``"module"``
        ``module_gene_col`` : str = ``"genes"``
        ``pathway_name_col`` : str = ``"MODULE"``
        ``pathway_gene_col`` : str = ``"GENE"``

        Why do ``module_name_col`` and ``pathway_name_col`` share the
        letters "module" / "MODULE" in different cases? Historical
        artifact. The default ``"module"`` (lowercase) matches the
        literal column header emitted by ``run_gene_membership``. The
        default ``"MODULE"`` (uppercase) matches the literal column
        header used by every pathway-DB file shipped under
        ``scGRNdb/Data/pathways/`` — where "MODULE" names the pathway
        (e.g. ``GOBP_GLYCOLYTIC_PROCESS``), not a graph-cluster.
        Internally, both are renamed to distinct canonical names
        (``MODULE`` vs. ``PATHWAY``) before the hypergeometric test
        runs. If you supply your own pathway file with differently-
        named columns, pass ``pathway_name_col`` / ``pathway_gene_col``
        explicitly.

        ``min_overlap`` : int = ``10``
            Note: the legacy ``enrichment.py`` argparse default is
            ``5``; this wrapper uses ``10`` to match the
            ``scing_config['enrichment']['min_overlap']`` default in
            ``config.py`` (which is what callers actually inherit).
            Override explicitly if you want the lower bound.
        ``pathway_size_min`` : int = ``10``
        ``pathway_size_max`` : int = ``500``
        ``n_background_genes`` : int = ``20000``

    Returns
    -------
    dict
        ``{"out_file": <out_file>, "n_pathways_enriched": <int>}``
        — the absolute output path and the number of (module, pathway)
        rows in the output. Zero is a valid result and means no output
        file was written (parity with the legacy script).

    Raises
    ------
    FileNotFoundError
        ``module_file``, ``pathway_file``, or ``out_file``'s parent
        directory does not exist.
    KeyError
        Required config key missing or a configured column name is
        absent from the input file.

    Notes
    -----
    The msigdb code path (``--use_msigdb`` in the legacy script) is
    intentionally **not** exposed by this wrapper. It is the only code
    path that imports ``decoupler``, and consumers of this API are
    expected to provide their own pathway files. Add a separate
    ``run_enrichment_decoupler(...)`` function later if the msigdb
    path is needed.

    No conda environment switching. Caller must already be in an env
    with ``pandas``, ``numpy``, ``polars``, and ``scipy`` installed.
    """
    _require_keys(config, _ENRICHMENT_REQUIRED, "run_enrichment")

    unknown = (
        set(config.keys())
        - set(_ENRICHMENT_REQUIRED)
        - set(_ENRICHMENT_OPTIONAL_DEFAULTS.keys())
    )
    if unknown:
        raise KeyError(
            f"run_enrichment: unknown config keys: {sorted(unknown)}. "
            f"Allowed: {sorted(set(_ENRICHMENT_REQUIRED) | set(_ENRICHMENT_OPTIONAL_DEFAULTS.keys()))}"
        )

    module_file = config["module_file"]
    pathway_file = config["pathway_file"]
    out_file = config["out_file"]

    _require_existing_file(module_file, "module_file")
    _require_existing_file(pathway_file, "pathway_file")
    _require_writable_parent(out_file, "out_file")

    kwargs = dict(_ENRICHMENT_OPTIONAL_DEFAULTS)
    for key in _ENRICHMENT_OPTIONAL_DEFAULTS:
        if key in config:
            kwargs[key] = config[key]

    from scnetworkatlas.python_files.enrichment import run as _run

    n_rows = _run(
        module_file=module_file,
        out_file=out_file,
        pathway_file=pathway_file,
        **kwargs,
    )

    return {
        "out_file": os.path.abspath(out_file),
        "n_pathways_enriched": int(n_rows),
    }
