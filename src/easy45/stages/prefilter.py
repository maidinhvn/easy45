"""S0 — deplete organelle reads.

Plastid DNA can dominate leaf WGS libraries. We map reads against the
user-supplied plastid+mito reference (minimap2 -x map-hifi) and *remove*
reads with strong organelle hits before recruitment, so they cannot inflate
coverage or leak into downstream clustering.

If no organelle reference is provided, this stage is a pass-through and emits
a warning — the eukaryote-trained boundary stage (S2) is the second line of
defence, but explicit depletion is strongly recommended.

Output keys: {"reads": <filtered FASTA/FASTQ>}
"""

from __future__ import annotations

import logging

from ..config import Config

log = logging.getLogger("easy45")


def run(config: Config, state: dict) -> dict:
    if config.organelle_ref is None:
        log.warning("S0: no --organelle-ref given; skipping organelle depletion")
        return {"reads": config.reads}
    log.warning("S0: prefilter not implemented yet (stub)")
    return {"reads": config.reads}
