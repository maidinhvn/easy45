"""Run configuration and default parameters for the easy45 pipeline.

All user-tunable thresholds live here so they can be surfaced on the CLI
and overridden in one place. The :class:`Config` object is built once in
``cli.py`` and threaded through every stage.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# Default recruitment anchor shipped with the package: one intact 45S repeat
# unit (18S-ITS1-5.8S-ITS2-26S-IGS) extracted from the Arabidopsis thaliana
# T2T NOR2 assembly (GenBank OR453402). Used when --anchor-ref is not given.
DEFAULT_ANCHOR = Path(__file__).parent / "data" / "default_anchor_arabidopsis_45S.fasta"

# Rfam covariance models (SSU RF01960 + LSU RF02543) used to define the mature
# 18S 5' / 26S 3' boundary on the consensus — a general, multi-taxon,
# structure-aware boundary, independent of any single reference sequence.
DEFAULT_CM = Path(__file__).parent / "data" / "rrna_ssu_lsu.cm"


@dataclass
class Config:
    # --- inputs -----------------------------------------------------------
    reads: Path                       # HiFi reads (FASTA/FASTQ, optionally gzipped)
    anchor_ref: Path = DEFAULT_ANCHOR # 45S anchor for recruitment (defaults to bundled Arabidopsis T2T unit)
    cm_ref: Path = DEFAULT_CM         # Rfam SSU+LSU CMs for mature-boundary trimming of the consensus
    organelle_ref: Path | None = None # plastid+mito genomes for Stage 0 depletion

    # --- outputs ----------------------------------------------------------
    outdir: Path = Path("easy45_out")
    workdir: Path | None = None       # intermediate files; defaults to <outdir>/work
    threads: int = 4

    # --- Stage 1: recruitment ---------------------------------------------
    recruit_preset: str = "map-hifi"  # minimap2 preset
    recruit_min_matchlen: int = 300   # min summed residue matches to the anchor

    # --- Stage 4: clustering (coarse grouping only) -----------------------
    cluster_id: float = 0.97          # vsearch --id; intentionally loose, fine-grained
                                      # variant calling happens in variant.py

    # --- variant calling (variant.py) -------------------------------------
    min_reads: int = 5                # a ribotype must be supported by >= this many reads
    min_freq: float = 0.05            # ... and represent >= this fraction of unit reads
    mask_homopolymer: int = 4         # mask indels inside homopolymer runs >= this length
    min_variant_sites: int = 2        # a variant ribotype must differ at >= this many SNV sites

    # --- Stage 6: IGS (best-effort) ---------------------------------------
    recover_igs: bool = True
    igs_length_bin: int = 50          # bp window for grouping IGS by length class

    # --- behaviour --------------------------------------------------------
    resume: bool = True               # skip stages whose outputs already exist
    keep_intermediate: bool = True

    def __post_init__(self) -> None:
        self.reads = Path(self.reads)
        self.anchor_ref = Path(self.anchor_ref)
        self.cm_ref = Path(self.cm_ref)
        if self.organelle_ref is not None:
            self.organelle_ref = Path(self.organelle_ref)
        self.outdir = Path(self.outdir)
        if self.workdir is None:
            self.workdir = self.outdir / "work"
        self.workdir = Path(self.workdir)
