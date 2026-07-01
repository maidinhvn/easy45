"""Command-line entry point for easy45."""

from __future__ import annotations

import argparse
import logging
import sys

from . import __version__
from .config import Config
from .external import DependencyError, check_dependencies
from .pipeline import run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="easy45",
        description="Assembly-free recovery of the 45S nrDNA transcribed unit "
                    "and ribotype variants from HiFi long reads.",
    )
    p.add_argument("--version", action="version", version=f"easy45 {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # --- run ----------------------------------------------------------------
    r = sub.add_parser("run", help="run the full pipeline")
    r.add_argument("-i", "--reads", required=True, help="HiFi reads (FASTA/FASTQ[.gz])")
    r.add_argument("-a", "--anchor-ref", default=None,
                   help="45S anchor for read recruitment "
                        "(default: bundled Arabidopsis T2T 45S unit)")
    r.add_argument("-r", "--organelle-ref", default=None,
                   help="plastid+mito genomes for Stage 0 depletion (strongly recommended)")
    r.add_argument("-o", "--outdir", default="easy45_out")
    r.add_argument("-t", "--threads", type=int, default=4)
    r.add_argument("--recruit-min-matchlen", type=int, default=300,
                   help="min summed bp matched to anchor to recruit a read")
    r.add_argument("--cluster-id", type=float, default=0.97)
    r.add_argument("--min-reads", type=int, default=5)
    r.add_argument("--min-freq", type=float, default=0.05)
    r.add_argument("--no-igs", action="store_true", help="skip IGS recovery")
    r.add_argument("--no-resume", action="store_true", help="ignore previous run state")
    r.add_argument("-v", "--verbose", action="store_true")

    # --- check-deps ---------------------------------------------------------
    c = sub.add_parser("check-deps", help="verify external tools are installed")
    c.add_argument("--optional", action="store_true",
                   help="also require optional tools (infernal)")
    return p


def _setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command == "check-deps":
        try:
            found = check_dependencies(include_optional=args.optional)
        except DependencyError as e:
            print(e, file=sys.stderr)
            return 1
        for tool, path in found.items():
            print(f"  {'OK ' if path else 'MISSING'}  {tool:10s} {path or ''}")
        print("All required tools found.")
        return 0

    if args.command == "run":
        _setup_logging(args.verbose)
        try:
            check_dependencies()
        except DependencyError as e:
            print(e, file=sys.stderr)
            return 1
        from .config import DEFAULT_ANCHOR
        config = Config(
            reads=args.reads,
            anchor_ref=args.anchor_ref or DEFAULT_ANCHOR,
            organelle_ref=args.organelle_ref,
            outdir=args.outdir,
            threads=args.threads,
            recruit_min_matchlen=args.recruit_min_matchlen,
            cluster_id=args.cluster_id,
            min_reads=args.min_reads,
            min_freq=args.min_freq,
            recover_igs=not args.no_igs,
            resume=not args.no_resume,
        )
        run_pipeline(config)
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
