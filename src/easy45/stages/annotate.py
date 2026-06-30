"""S7 — annotate final ribotype sequences and write the report.

Runs ITSx once on the final ribotype sequences (primary consensus + any
variants) to delimit SSU(18S) / ITS1 / 5.8S / ITS2 / LSU(28S) precisely, emits
those as GFF3, and writes a human-readable run report (ribotype table summary,
IGS status, and warnings such as skipped organelle depletion or a hybrid
signal).

Input keys:  {"consensus","variants","ribotypes", optional "igs"}
Output keys: {"gff","report"}
"""

from __future__ import annotations

import logging
import re
import subprocess

from Bio import SeqIO

from ..config import Config

log = logging.getLogger("easy45")

# ITSx region label -> GFF feature name
_REGION = {"SSU": "18S_rRNA", "ITS1": "ITS1", "5.8S": "5_8S_rRNA",
           "ITS2": "ITS2", "LSU": "28S_rRNA"}
_POS_RE = re.compile(r"(SSU|ITS1|5\.8S|ITS2|LSU): (\d+)-(\d+)")


def _combine_inputs(consensus_fa, variants_fa, dest):
    n = 0
    with open(dest, "w") as out:
        for src in (consensus_fa, variants_fa):
            if src and src.exists():
                txt = src.read_text()
                out.write(txt)
                n += txt.count(">")
    return n


def _parse_itsx_positions(pos_path):
    """Yield (seqname, [(feature, start, end), ...]) from ITSx .positions.txt."""
    with open(pos_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            name = parts[0]
            feats = []
            for region, s, e in _POS_RE.findall(line):
                feats.append((_REGION[region], int(s), int(e)))
            yield name, feats


def run(config: Config, state: dict) -> dict:
    combined = config.workdir / "s7_ribotypes.fasta"
    n_seqs = _combine_inputs(state["consensus"], state.get("variants"), combined)

    prefix = config.workdir / "s7_itsx"
    gff = config.outdir / "annotation.gff3"
    report = config.outdir / "report.txt"
    its = config.outdir / "its.fasta"
    its_parts = config.outdir / "its_parts.fasta"

    if n_seqs == 0:
        # No ribotype passed the abundance gate (e.g. highly heterogeneous or
        # contaminated units that did not cluster into any dominant group).
        # Emit empty outputs + a report instead of crashing on a missing ITSx run.
        log.warning("S7: no ribotype to annotate (no cluster passed the abundance "
                    "gate); writing empty annotation/ITS + report")
        gff.write_text("##gff-version 3\n")
        its.write_text("")
        its_parts.write_text("")
        _write_report(config, state, report, 0, 0)
        return {"gff": gff, "report": report, "its": its, "its_parts": its_parts}

    log.info("S7: annotating %d ribotype sequence(s) with ITSx", n_seqs)
    subprocess.run(
        ["ITSx", "-i", str(combined), "-o", str(prefix),
         "-t", "Tracheophyta", "--cpu", str(config.threads),
         "--preserve", "T", "--save_regions", "none", "--graphical", "F"],
        check=True, capture_output=True, text=True,
    )

    # --- GFF3 ---
    pos_file = f"{prefix}.positions.txt"
    n_features = 0
    with open(gff, "w") as g:
        g.write("##gff-version 3\n")
        for name, feats in _parse_itsx_positions(pos_file):
            for feat, s, e in feats:
                if e < s:
                    continue
                g.write(f"{name}\teasy45\t{feat}\t{s}\t{e}\t.\t+\t.\tName={feat}\n")
                n_features += 1

    # --- ITS barcode outputs ---
    n_its = _write_its(combined, pos_file, its, its_parts)

    # --- report ---
    _write_report(config, state, report, n_seqs, n_features)
    log.info("S7: wrote %s (%d features), %s (%d ITS barcode(s)) and %s",
             gff.name, n_features, its.name, n_its, report.name)
    return {"gff": gff, "report": report, "its": its, "its_parts": its_parts}


def _write_its(ribotypes_fa, pos_file, its_fa, its_parts_fa):
    """Build ITS barcode outputs from ITSx positions + ribotype sequences.

    its.fasta       : one record per ribotype = ITS1-5.8S-ITS2 contiguous
                      (the standard ITS barcode used for BLAST/GenBank).
    its_parts.fasta : ITS1, 5.8S and ITS2 as separate records per ribotype
                      (for ITS2-only barcoding / separate alignment).
    Returns the number of ribotypes that yielded a full ITS region.
    """
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(ribotypes_fa, "fasta")}
    n = 0
    with open(its_fa, "w") as full, open(its_parts_fa, "w") as parts:
        for name, feats in _parse_itsx_positions(pos_file):
            seq = seqs.get(name)
            if seq is None:
                continue
            coords = {f: (s, e) for f, s, e in feats}
            i1, c58, i2 = coords.get("ITS1"), coords.get("5_8S_rRNA"), coords.get("ITS2")
            if not (i1 and i2 and i1[0] <= i2[1]):
                continue
            full.write(f">{name}_ITS ITS1-5.8S-ITS2\n{seq[i1[0] - 1:i2[1]]}\n")
            n += 1
            for label, span in (("ITS1", i1), ("5.8S", c58), ("ITS2", i2)):
                if span:
                    parts.write(f">{name}_{label}\n{seq[span[0] - 1:span[1]]}\n")
    return n


def _write_report(config, state, report_path, n_ribotypes, n_features):
    lines = ["easy45 run report", "=" * 40, f"reads:  {config.reads}",
             f"anchor: {config.anchor_ref.name}", ""]

    # ribotype table
    lines.append("Ribotypes:")
    with open(state["ribotypes"]) as fh:
        header = fh.readline()
        rows = [r.rstrip("\n").split("\t") for r in fh]
    primary = [r for r in rows if r[3] == "primary"]
    variants = [r for r in rows if r[3] == "variant"]
    foreign = [r for r in rows if r[3] == "foreign"]
    for r in primary + variants:
        extra = f", {r[5]} subs vs primary" if r[3] == "variant" else ""
        lines.append(f"  {r[3]:8s} cluster {r[0]}: {r[1]} reads "
                     f"({float(r[2])*100:.1f}%), len {r[4]}{extra}")
    lines.append(f"  ({len(rows)} clusters total; "
                 f"{len(variants)} variant(s) above threshold)")
    if foreign:
        lines.append(f"  {len(foreign)} divergent rDNA cluster(s) excluded "
                     "(<90% identity to primary; likely contaminant/endophyte "
                     "or pseudogene)")
    lines.append("")

    # IGS
    igs = state.get("igs")
    if igs and igs.exists() and igs.read_text().strip():
        n_igs = igs.read_text().count(">")
        lines.append(f"IGS: {n_igs} complete spacer(s) recovered")
    else:
        lines.append("IGS: no complete IGS recovered")
    lines.append("")

    # warnings
    warns = []
    if not primary:
        warns.append("NO ribotype passed the abundance gate (min_reads/min_freq); "
                     "units did not cluster into a dominant group - possible heavy "
                     "contamination / heterogeneity; see ribotypes.tsv")
    if config.organelle_ref is None:
        warns.append("organelle depletion was SKIPPED (no --organelle-ref)")
    if len(variants) >= 1:
        warns.append(f"{len(variants)} variant ribotype(s) detected — "
                     "possible hybrid/allopolyploid; inspect variants.fasta")
    lines.append("Warnings:")
    lines.extend(f"  - {w}" for w in warns) if warns else lines.append("  none")

    report_path.write_text("\n".join(lines) + "\n")
