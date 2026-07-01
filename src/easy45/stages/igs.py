"""S6 — IGS recovery (best-effort, complete only).

The intergenic spacer (IGS/NTS+ETS) lies between the 3' end of one unit's 28S
and the 5' start of the next unit's 18S. We recover it only from reads that
span BOTH flanks cleanly (a 28S and a following 18S on the same strand, with a
plausible gap). Reads that do not fully span the spacer are DISCARDED — easy45
never emits a partial/flagged IGS.

Complete IGS sequences are strand-normalised, binned by length class
(config.igs_length_bin, since IGS length is polymorphic) and a POA consensus is
built per bin (POA is never run across length classes, which would fabricate
gaps). If no read spans an IGS, igs.fasta is empty.

Input keys:  {"recruited","barrnap_gff"}
Output keys: {"igs"}
"""

from __future__ import annotations

import logging
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq

from ..config import Config

log = logging.getLogger("easy45")

IGS_MIN, IGS_MAX = 300, 10000   # plausible bp span of one intergenic spacer
CONSENSUS_MAX_READS = 100


def _features_by_read(gff_path):
    reads: dict[str, dict] = {}
    want = {"18S_rRNA": "18S", "28S_rRNA": "28S"}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "rRNA":
                continue
            name = next((a.split("=", 1)[1] for a in c[8].split(";")
                         if a.startswith("Name=")), "")
            gene = want.get(name)
            if gene:
                reads.setdefault(c[0], {"18S": [], "28S": []})[gene].append(
                    (int(c[3]), int(c[4]), c[6]))
    return reads


def _find_igs(feats):
    """Return (lo, hi, strand) of one complete IGS (28S 3' -> next 18S 5'), or None."""
    best = None
    for s28, e28, st28 in feats["28S"]:
        for s18, e18, st18 in feats["18S"]:
            if st18 != st28:
                continue
            if st28 == "+" and s18 > e28:          # IGS between 28S.end and next 18S.start
                lo, hi, gap = e28 + 1, s18 - 1, s18 - e28
            elif st28 == "-" and e18 < s28:        # reversed: 18S.end .. 28S.start
                lo, hi, gap = e18 + 1, s28 - 1, s28 - e18
            else:
                continue
            if IGS_MIN <= gap <= IGS_MAX and (best is None or gap < best[3]):
                best = (lo, hi, st28, gap)
    if best is None:
        return None
    return best[0], best[1], best[2]


def _abpoa(seqs, tmp_fasta):
    with open(tmp_fasta, "w") as f:
        for i, s in enumerate(seqs):
            f.write(f">{i}\n{s}\n")
    proc = subprocess.run(["abpoa", str(tmp_fasta)],
                          check=True, capture_output=True, text=True)
    return "".join(ln for ln in proc.stdout.splitlines() if not ln.startswith(">"))


def run(config: Config, state: dict) -> dict:
    feats = _features_by_read(state["barrnap_gff"])
    seqs = SeqIO.index(str(state["recruited"]), "fasta")
    igs_fa = config.outdir / "igs.fasta"

    # collect complete IGS sequences (strand-normalised)
    igs_seqs = []
    for rid, f in feats.items():
        found = _find_igs(f)
        if found is None:
            continue
        lo, hi, strand = found
        sub = seqs[rid].seq[lo - 1:hi]
        if strand == "-":
            sub = sub.reverse_complement()
        igs_seqs.append(str(sub))
    seqs.close()

    # bin by length class, consensus within each bin
    bins: dict[int, list[str]] = {}
    for s in igs_seqs:
        bins.setdefault(len(s) // config.igs_length_bin, []).append(s)

    # keep only well-supported length classes (discard low-count noise tail)
    kept = {b: m for b, m in bins.items() if len(m) >= config.min_reads}

    tmp = config.workdir / "s6_abpoa_input.fasta"
    n_bins = 0
    with open(igs_fa, "w") as out:
        for b, members in sorted(kept.items(), key=lambda kv: len(kv[1]), reverse=True):
            cons = members[0] if len(members) == 1 else _abpoa(members[:CONSENSUS_MAX_READS], tmp)
            out.write(f">igs_lenclass~{b*config.igs_length_bin} reads={len(members)} "
                      f"len={len(cons)}\n{cons}\n")
            n_bins += 1

    log.info("S6: %d complete IGS from reads -> %d length class(es) "
             "with >=%d reads (%d sparse classes discarded)",
             len(igs_seqs), n_bins, config.min_reads, len(bins) - len(kept))
    return {"igs": igs_fa}
