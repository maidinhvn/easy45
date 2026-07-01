"""Smoke tests for the easy45 skeleton — verify wiring, not stage logic."""

from easy45 import __version__
from easy45.config import Config
from easy45.cli import _build_parser
from easy45.pipeline import STAGES


def test_version():
    assert __version__ == "0.1.0"


def test_parser_builds():
    parser = _build_parser()
    args = parser.parse_args(["run", "-i", "r.fq", "-a", "anchor.fa"])
    assert args.command == "run"
    assert args.reads == "r.fq"


def test_config_defaults():
    cfg = Config(reads="r.fq", anchor_ref="anchor.fa")
    assert cfg.workdir == cfg.outdir / "work"
    assert cfg.recover_igs is True


def test_all_stages_callable():
    for key, label, fn in STAGES:
        assert callable(fn), key
