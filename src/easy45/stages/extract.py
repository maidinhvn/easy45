"""S3 - cut and orient transcribed units.

For each spanning read we pair an 18S with its downstream 28S (the pair whose
genomic span falls in one-unit range) and excise that unit, strand-normalised
to 18S->26S orientation. Reads with more than one unit are handled (one unit
is cut, not the whole 18S..28S stretch across copies).

Two FASTAs are written, both keyed by read id so membership maps across them:
  * units (core, barrnap boundaries) -> clustering substrate (S4). No margin,
    so the variable ETS does not fragment clusters.
  * units_margin (core +/- EDGE_MARGIN bp) -> consensus substrate (S5), giving
    the CM trim (S5) room to recover the mature 18S 5'/26S 3' termini that
    barrnap clips.

Input keys:  {"recruited": <FASTA>, "barrnap_gff": <GFF>, "spanning_ids": <ids>}
Output keys: {"units": <core FASTA>, "units_margin": <margin FASTA>}
"""

from __future__ import annotations

import logging

from Bio import SeqIO

from ..config import Config

log = logging.getLogger("easy45")

UNIT_MIN, UNIT_MAX = 4000, 9000   # plausible bp span of one transcribed unit (18S start..28S end)
ITS_MAX_GAP = 1500                # max bp between 18S and 28S within one unit (ITS1+5.8S+ITS2);
                                  # excludes IGS-spanning 28S->next-18S mis-pairs (short-IGS taxa)
# barrnap (HMM) trims a few-to-tens of bp off the 18S 5' / 26S 3' termini. We
# also emit a margin-widened unit so the true mature ends are inside it; the
# consensus is trimmed back to the precise mature boundary by CMs in S5.
EDGE_MARGIN = 150


def _features_by_read(gff_path):
    """Return {read_id: {'18S': [(s,e,strand)...], '28S': [...]}}, all hits."""
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
            if gene is None:
                continue
            reads.setdefault(c[0], {"18S": [], "28S": []})[gene].append(
                (int(c[3]), int(c[4]), c[6]))
    return reads


def _pick_unit(feats):
    """Choose one (lo, hi, strand) transcribed unit from a read's 18S/28S hits.

    A valid transcribed unit is 18S--ITS1--5.8S--ITS2--26S, i.e. the 18S lies 5'
    of the 28S within ONE repeat, separated only by the ITS region. We therefore
    require the 18S to be upstream of the 28S (transcription order) with an
    ITS-sized gap (<= ITS_MAX_GAP). This rejects the IGS-spanning mis-pair
    (a 28S followed by the NEXT repeat's 18S), which otherwise sneaks in when the
    IGS is short enough that 28S->next-18S still falls in one-unit span range.
    Among valid pairs, prefer the smallest span (one copy, not two).
    """
    best = None
    for s18, e18, st18 in feats["18S"]:
        for s28, e28, st28 in feats["28S"]:
            if st18 != st28:
                continue
            if st18 == "+":            # 18S then 28S; gap = ITS region
                gap, lo, hi = s28 - e18, s18, e28
            else:                      # minus strand: 28S upstream in genomic coords
                gap, lo, hi = s18 - e28, s28, e18
            if not (-100 <= gap <= ITS_MAX_GAP):   # ITS-sized, not IGS-sized
                continue
            span = hi - lo + 1
            if UNIT_MIN <= span <= UNIT_MAX and (best is None or span < best[3]):
                best = (lo, hi, st18, span)
    if best is None:
        return None
    lo, hi, strand, _ = best
    return lo, hi, strand


def run(config: Config, state: dict) -> dict:
    recruited = state["recruited"]
    feats = _features_by_read(state["barrnap_gff"])
    spanning = set(state["spanning_ids"].read_text().split())
    units_fa = config.workdir / "s3_units.fasta"
    margin_fa = config.workdir / "s3_units_margin.fasta"

    seqs = SeqIO.index(str(recruited), "fasta")
    n_written = n_skipped = 0
    with open(units_fa, "w") as core, open(margin_fa, "w") as marg:
        for rid in spanning:
            f = feats.get(rid)
            if not f:
                n_skipped += 1
                continue
            picked = _pick_unit(f)
            if picked is None:
                n_skipped += 1
                continue
            lo, hi, strand = picked
            seq = seqs[rid].seq
            mlo = max(1, lo - EDGE_MARGIN)
            mhi = min(len(seq), hi + EDGE_MARGIN)
            core_sub = seq[lo - 1:hi]                 # GFF coords are 1-based inclusive
            marg_sub = seq[mlo - 1:mhi]
            if strand == "-":
                core_sub = core_sub.reverse_complement()
                marg_sub = marg_sub.reverse_complement()
            # id = read id (shared across both files + membership); coords in description
            core.write(f">{rid} unit={lo}-{hi} strand={strand}\n{core_sub}\n")
            marg.write(f">{rid} unit={mlo}-{mhi} strand={strand}\n{marg_sub}\n")
            n_written += 1
    seqs.close()

    log.info("S3: cut %d transcribed units (core + margin); %d spanning reads "
             "skipped (no single-unit 18S-28S pair)", n_written, n_skipped)
    if n_written == 0:
        raise RuntimeError("S3 produced 0 units - check S2 boundaries")
    return {"units": units_fa, "units_margin": margin_fa}
