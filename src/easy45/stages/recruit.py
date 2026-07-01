"""S1 — recruit candidate rDNA reads.

Map (organelle-depleted) reads against the conserved 18S/5.8S/26S anchor with
``minimap2 -x map-hifi`` and keep reads carrying a confident conserved hit.
The anchor's coding regions are highly conserved across angiosperms, while its
ITS/ETS/IGS portions are too divergent to align across genera — so in practice
recruitment is driven by the conserved genes and is naturally nuclear-specific
(bacterial-type plastid 16S/23S do not pass map-hifi thresholds).

A read is recruited when its summed residue matches to the anchor reach
``config.recruit_min_matchlen`` (a read may have several PAF lines, e.g. one per
gene; matches are summed per read).

Input keys:  {"reads": <FASTA/FASTQ[.gz]>}
Output keys: {"recruited": <FASTA of candidate reads>}
"""

from __future__ import annotations

import logging
import subprocess

from ..config import Config
from ..external import run as sh
from ..io import detect_format

log = logging.getLogger("easy45")


def _recruited_ids(paf_path, min_matchlen: int) -> set[str]:
    """Sum PAF residue matches (col 10, 0-based index 9) per query read."""
    matches: dict[str, int] = {}
    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            matches[f[0]] = matches.get(f[0], 0) + int(f[9])
    return {q for q, m in matches.items() if m >= min_matchlen}


def run(config: Config, state: dict) -> dict:
    reads = state.get("reads", config.reads)
    paf = config.workdir / "s1_recruit.paf"
    ids_file = config.workdir / "s1_recruited.ids"
    recruited = config.workdir / "s1_recruited.fasta"

    # 1. map reads (query) against the anchor (target); PAF -> file (not memory)
    log.info("S1: mapping reads vs anchor (%s)", config.anchor_ref.name)
    with open(paf, "w") as out:
        subprocess.run(
            ["minimap2", "-x", config.recruit_preset, "-t", str(config.threads),
             str(config.anchor_ref), str(reads)],
            check=True, stdout=out, stderr=subprocess.PIPE, text=True,
        )

    # 2. select reads with enough conserved match
    ids = _recruited_ids(paf, config.recruit_min_matchlen)
    ids_file.write_text("\n".join(sorted(ids)) + ("\n" if ids else ""))
    log.info("S1: recruited %d reads (>=%d bp matched to anchor)",
             len(ids), config.recruit_min_matchlen)
    if not ids:
        raise RuntimeError("S1 recruited 0 reads — check anchor/reads or lower "
                           "--recruit-min-matchlen")

    # 3. extract recruited reads, normalised to FASTA
    #    (seqkit grep preserves input format; convert FASTQ -> FASTA via a pipe)
    if detect_format(reads) == "fastq":
        grep = subprocess.Popen(
            ["seqkit", "grep", "-f", str(ids_file), str(reads)],
            stdout=subprocess.PIPE)
        with open(recruited, "w") as out:
            subprocess.run(["seqkit", "fq2fa"], stdin=grep.stdout,
                           stdout=out, check=True)
        grep.stdout.close()
        if grep.wait() != 0:
            raise RuntimeError("seqkit grep failed during recruitment")
    else:
        sh(["seqkit", "grep", "-f", str(ids_file), str(reads), "-o", str(recruited)])

    return {"recruited": recruited}
