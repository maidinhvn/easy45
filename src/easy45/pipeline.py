"""Pipeline orchestration: chains the stages and handles resume.

Each stage is a callable ``fn(config, inputs) -> outputs`` living in
``easy45.stages``. The pipeline records completed stages in a manifest file
inside the work directory so a re-run can skip finished stages (``--resume``).
"""

from __future__ import annotations

import json
import logging

from .config import Config
from .stages import (
    annotate,
    boundary,
    cluster,
    consensus,
    extract,
    igs,
    prefilter,
    recruit,
)

log = logging.getLogger("easy45")

# Ordered (key, human label, callable). Each stage reads/writes files under
# config.workdir and returns a dict of named output paths.
STAGES = [
    ("prefilter", "S0 deplete organelle reads", prefilter.run),
    ("recruit",   "S1 recruit rDNA reads",      recruit.run),
    ("boundary",  "S2 locate unit boundaries",  boundary.run),
    ("extract",   "S3 cut & orient units",      extract.run),
    ("cluster",   "S4 coarse cluster ribotypes", cluster.run),
    ("consensus", "S5 consensus + variants",    consensus.run),
    ("igs",       "S6 IGS (best-effort)",       igs.run),
    ("annotate",  "S7 annotate & report",       annotate.run),
]


def _manifest_path(config: Config):
    return config.workdir / "manifest.json"


def _load_manifest(config: Config) -> dict:
    p = _manifest_path(config)
    if p.exists():
        return json.loads(p.read_text())
    return {}


def _save_manifest(config: Config, manifest: dict) -> None:
    _manifest_path(config).write_text(json.dumps(manifest, indent=2, default=str))


def run_pipeline(config: Config) -> dict:
    """Execute all stages in order, returning the accumulated outputs."""
    config.outdir.mkdir(parents=True, exist_ok=True)
    config.workdir.mkdir(parents=True, exist_ok=True)

    manifest = _load_manifest(config) if config.resume else {}
    state: dict = manifest.get("outputs", {})

    for key, label, fn in STAGES:
        if key == "igs" and not config.recover_igs:
            log.info("[skip] %s (--no-igs)", label)
            continue
        if config.resume and key in manifest.get("done", []):
            log.info("[resume] %s already complete", label)
            continue
        log.info("[run] %s", label)
        outputs = fn(config, state)
        state.update(outputs or {})
        manifest.setdefault("done", []).append(key)
        manifest["outputs"] = state
        _save_manifest(config, manifest)

    log.info("Pipeline complete. Results in %s", config.outdir)
    return state
