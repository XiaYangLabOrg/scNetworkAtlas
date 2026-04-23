"""scnetworkatlas — importable surface for the scNetworkAtlas pipeline.

This package exposes a small, in-process API equivalent to the relevant
``run_pipeline.py`` subcommands. It exists so that downstream consumers
(e.g. the secure mergeomicsSecure backend) can run scNetworkAtlas steps
without spawning ``python run_pipeline.py <step>`` as a subprocess and
without reading a generated ``config.py`` from disk.

Public surface (see ``scnetworkatlas.api`` for full docs):

    run_gene_membership(config: dict) -> dict   # ≡ run_pipeline.py gene_membership
    run_enrichment(config: dict) -> dict        # ≡ run_pipeline.py enrichment

The underlying algorithms are unchanged — these wrappers call the same
functions in ``scnetworkatlas/python_files/`` that ``run_pipeline.py``
calls today via its bash shims.
"""

from scnetworkatlas.api import run_gene_membership, run_enrichment

__all__ = ["run_gene_membership", "run_enrichment"]
