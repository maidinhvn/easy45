"""S2 — locate unit boundaries on each read.

Run ``barrnap --kingdom euk`` on the recruited reads to get per-read
coordinates and strand of the 18S, 5.8S and 28S rRNA genes. A read is flagged
as *spanning* the transcribed unit when it contains a complete (non-partial)
18S AND a complete 28S on the same strand — i.e. the ETS->26S unit is fully
inside the read. These coordinates drive cutting/orientation in S3.

(ITS1/ITS2 delimitation with ITSx is deferred to S7, where it runs once on the
clean consensus rather than on every noisy read.)

Input keys:  {"recruited": <FASTA>}
Output keys: {"boundaries": <TSV per-read feature coords>,
              "spanning_ids": <text file, one read id per line>}
"""

from __future__ import annotations

import logging
import subprocess

from ..config import Config

log = logging.getLogger("easy45")

# barrnap euk feature names we care about (5S is a separate locus -> ignored)
_GENES = {"18S_rRNA": "18S", "5_8S_rRNA": "5_8S", "28S_rRNA": "28S"}


def _parse_barrnap(gff_path):
    """Return {read_id: {gene: (start, end, strand, is_partial)}}.

    Keeps the longest hit per gene per read (best/most-complete copy).
    """
    reads: dict[str, dict] = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "rRNA":
                continue
            seqid, start, end, strand, attrs = c[0], int(c[3]), int(c[4]), c[6], c[8]
            name = next((a.split("=", 1)[1] for a in attrs.split(";")
                         if a.startswith("Name=")), "")
            gene = _GENES.get(name)
            if gene is None:
                continue
            is_partial = "partial" in attrs.lower()
            rec = reads.setdefault(seqid, {})
            prev = rec.get(gene)
            if prev is None or (end - start) > (prev[1] - prev[0]):
                rec[gene] = (start, end, strand, is_partial)
    return reads


def run(config: Config, state: dict) -> dict:
    recruited = state["recruited"]
    gff = config.workdir / "s2_barrnap.gff"
    tsv = config.workdir / "s2_boundaries.tsv"
    span_ids = config.workdir / "s2_spanning.ids"

    log.info("S2: locating rRNA genes with barrnap (--kingdom euk)")
    with open(gff, "w") as out:
        subprocess.run(
            ["barrnap", "--kingdom", "euk", "--threads", str(config.threads),
             str(recruited)],
            check=True, stdout=out, stderr=subprocess.PIPE, text=True,
        )

    reads = _parse_barrnap(gff)

    spanning = []
    with open(tsv, "w") as t:
        t.write("read_id\tstrand\thas_18S\thas_5_8S\thas_28S\tfull_18S\tfull_28S"
                "\tunit_start\tunit_end\tspanning\n")
        for rid, genes in reads.items():
            g18, g58, g28 = genes.get("18S"), genes.get("5_8S"), genes.get("28S")
            full18 = bool(g18 and not g18[3])
            full28 = bool(g28 and not g28[3])
            strand = g18[2] if g18 else (g28[2] if g28 else ".")
            same_strand = g18 and g28 and g18[2] == g28[2]
            is_span = bool(full18 and full28 and same_strand)
            # transcribed unit extent = 18S start .. 28S end (strand-aware)
            ustart = uend = "."
            if g18 and g28:
                lo = min(g18[0], g28[0]); hi = max(g18[1], g28[1])
                ustart, uend = lo, hi
            t.write(f"{rid}\t{strand}\t{bool(g18)}\t{bool(g58)}\t{bool(g28)}"
                    f"\t{full18}\t{full28}\t{ustart}\t{uend}\t{is_span}\n")
            if is_span:
                spanning.append(rid)

    span_ids.write_text("\n".join(spanning) + ("\n" if spanning else ""))
    log.info("S2: %d reads carry rRNA genes; %d span a complete transcribed unit",
             len(reads), len(spanning))
    if not spanning:
        raise RuntimeError("S2 found 0 reads spanning a complete unit — "
                           "check recruited reads / anchor")

    return {"boundaries": tsv, "spanning_ids": span_ids, "barrnap_gff": gff}
