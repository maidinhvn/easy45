"""Ribotype variant validation - the scientific core of easy45.

A cluster passing the abundance gate (S5: min_reads + min_freq) is only a
*candidate*. Comparing its consensus to the primary ribotype, three outcomes:

  noise    - differs only by homopolymer/STR indels and < min_variant_sites
             substitutions. HiFi residual error is dominated by *systematic*
             indels in homopolymer tracts that recur across reads, so an
             abundance filter alone cannot remove them; we judge on
             SUBSTITUTION sites, masking mononucleotide indels.

  foreign  - aligns to the primary at low identity (< MIN_VARIANT_IDENTITY) or
             over little of its length. A real intragenomic ribotype variant is
             kept >~95% identical to the primary by concerted evolution; a 45S
             that is only ~76% identical is a different organism (e.g. a fungal
             endophyte recruited via the conserved 18S/26S) or a pseudogene -
             NOT a ribotype variant. Calling these "variants" would falsely
             flag a hybrid, so they are excluded.

  variant  - highly similar to the primary (>= MIN_VARIANT_IDENTITY, well
             covered) yet differs at >= min_variant_sites real substitution
             sites: a genuine ribotype (hybrid / allopolyploid signal).

Candidate and primary consensuses are aligned with minimap2 --cs.
"""

from __future__ import annotations

import logging
import re
import subprocess
from dataclasses import dataclass, field

log = logging.getLogger("easy45")

MIN_VARIANT_IDENTITY = 0.90   # a real ribotype variant must be >= this identical to primary
MIN_ALIGNED_FRAC = 0.80       # ... over at least this fraction of its length


@dataclass
class VariantVerdict:
    category: str                    # "variant" | "noise" | "foreign"
    substitutions: int
    homopolymer_indels: int
    other_indels: int
    identity: float
    aligned_frac: float
    sub_positions: list = field(default_factory=list)

    @property
    def is_real(self) -> bool:
        return self.category == "variant"


def _best_alignment(primary_fa, variant_fa, tmp_paf):
    """Align variant (query) to primary; return (cs, matches, block_len, qlen, qcov)."""
    with open(tmp_paf, "w") as out:
        subprocess.run(["minimap2", "-c", "--cs=short", str(primary_fa), str(variant_fa)],
                       check=True, stdout=out, stderr=subprocess.DEVNULL, text=True)
    best = None
    best_match = -1
    with open(tmp_paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            matches, block = int(f[9]), int(f[10])
            qlen, qstart, qend = int(f[1]), int(f[2]), int(f[3])
            cs = next((t[5:] for t in f[12:] if t.startswith("cs:Z:")), None)
            if cs is not None and matches > best_match:
                best_match = matches
                best = (cs, matches, block, qlen, (qend - qstart) / qlen)
    return best


def _parse_cs(cs):
    """Walk a short cs string -> (n_subs, n_hp_indels, n_other_indels, sub_positions).

    A mononucleotide insertion/deletion is treated as homopolymer/STR-slip noise
    and masked; substitutions are always counted as real divergence.
    """
    subs = hp_indels = other_indels = 0
    positions = []
    refpos = 0
    for op, val in re.findall(r"([:*+\-])([A-Za-z0-9]+)", cs):
        if op == ":":
            refpos += int(val)
        elif op == "*":
            subs += 1
            positions.append(refpos + 1)
            refpos += 1
        elif op == "-":
            hp_indels += 1 if len(set(val.lower())) == 1 else 0
            other_indels += 0 if len(set(val.lower())) == 1 else 1
            refpos += len(val)
        elif op == "+":
            hp_indels += 1 if len(set(val.lower())) == 1 else 0
            other_indels += 0 if len(set(val.lower())) == 1 else 1
    return subs, hp_indels, other_indels, positions


def confirm_variant(primary_fa, variant_fa, tmp_paf,
                    min_variant_sites: int) -> VariantVerdict:
    """Classify a candidate consensus as variant / noise / foreign."""
    aln = _best_alignment(primary_fa, variant_fa, tmp_paf)
    if aln is None:
        # no alignment to primary at all -> a different sequence entirely
        return VariantVerdict("foreign", 0, 0, 0, 0.0, 0.0)
    cs, matches, block, qlen, qcov = aln
    identity = matches / block if block else 0.0
    subs, hp, other, pos = _parse_cs(cs)

    if identity < MIN_VARIANT_IDENTITY or qcov < MIN_ALIGNED_FRAC:
        category = "foreign"          # contaminant / pseudogene, not a ribotype variant
    elif subs >= min_variant_sites:
        category = "variant"
    else:
        category = "noise"            # homopolymer/minor only
    return VariantVerdict(category, subs, hp, other, identity, qcov, pos)
