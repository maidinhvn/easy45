#!/usr/bin/env python3
"""Reproducibly build the default easy45 recruitment anchor.

Fetches a slice of the Arabidopsis thaliana T2T NOR2 array (GenBank OR453402)
and extracts one complete 45S repeat unit. See
src/easy45/data/ANCHOR_PROVENANCE.md for full provenance and citation.

Usage:  python scripts/make_default_anchor.py
Output: src/easy45/data/default_anchor_arabidopsis_45S.fasta

No third-party dependencies (stdlib urllib only).
"""

from __future__ import annotations

import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NOR2_ACCESSION = "OR453402"   # A. thaliana T2T NOR2 (Fultz et al. 2023, Sci Adv)
NOR2_SLICE_BP = 60000         # enough to contain several repeat units
ANCHOR_18S = "X16077"         # A. thaliana 18S rRNA gene; gene = bases 88..1891
PROBE_LEN = 40                # 5' 18S k-mer used to find unit boundaries

OUT = Path(__file__).resolve().parent.parent / "src" / "easy45" / "data" \
    / "default_anchor_arabidopsis_45S.fasta"


def efetch_fasta(accession: str, start: int | None = None, stop: int | None = None) -> str:
    params = {"db": "nucleotide", "id": accession, "rettype": "fasta", "retmode": "text"}
    if start and stop:
        params.update(seq_start=str(start), seq_stop=str(stop))
    url = EUTILS + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url) as r:               # noqa: S310 (trusted host)
        text = r.read().decode()
    return "".join(line for line in text.splitlines() if not line.startswith(">")).upper()


def main() -> None:
    print(f"Fetching {ANCHOR_18S} (18S anchor) ...")
    g18s = efetch_fasta(ANCHOR_18S)[88 - 1:1891]         # GenBank gene coords
    probe = g18s[:PROBE_LEN]

    print(f"Fetching {NOR2_ACCESSION}:1-{NOR2_SLICE_BP} (NOR2 slice) ...")
    nor2 = efetch_fasta(NOR2_ACCESSION, 1, NOR2_SLICE_BP)

    # locate 18S 5' starts
    starts, i = [], nor2.find(probe)
    while i != -1:
        starts.append(i)
        i = nor2.find(probe, i + 1)
    print(f"18S start positions: {starts}")

    # pick first consecutive pair whose spacing matches the canonical period (~10 kb),
    # skipping anomalous array-edge hits
    unit = None
    for a, b in zip(starts, starts[1:]):
        if 9500 <= (b - a) <= 10600:
            unit = nor2[a:b]
            # 1-based inclusive GenBank coords: start = a+1, end = b (b is the
            # 0-based exclusive slice end, i.e. 1-based inclusive last base).
            coords = f"{NOR2_ACCESSION}.1:{a + 1}-{b}"
            break
    if unit is None:
        raise SystemExit("No canonical ~10 kb repeat unit found; check NOR2 slice.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        f.write(f">45S_unit Arabidopsis_thaliana T2T_NOR2 {coords} one_complete_repeat\n")
        for j in range(0, len(unit), 70):
            f.write(unit[j:j + 70] + "\n")
    print(f"Wrote {OUT} ({len(unit)} bp, {coords})")


if __name__ == "__main__":
    main()
