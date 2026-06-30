"""S5 - consensus per cluster + ribotype calling.

For each cluster that passes the abundance gate (>= config.min_reads AND
>= config.min_freq), build a POA consensus with abpoa. The largest cluster's
consensus is the primary ribotype (consensus.fasta). Each other passing cluster
is a *candidate* that must additionally clear the homopolymer-aware divergence
check in variant.py (>= config.min_variant_sites substitution differences from
the primary) to be written to variants.fasta; candidates that differ only by
homopolymer/STR indels are recorded as noise, not variants.

Consensus is built from the margin-widened units (units_margin) so the CM trim
(Rfam SSU+LSU via cmsearch) can recover the mature 18S 5' / 26S 3' termini;
clustering (S4) used the core units.

abpoa is fed at most CONSENSUS_MAX_READS reads per cluster (HiFi ~Q30, so ~100
reads pin the consensus and POA stays tractable on the 3000+ read primary).

Input keys:  {"units_margin","membership","cluster_freq"}
Output keys: {"consensus","variants","ribotypes"}
"""

from __future__ import annotations

import logging
import subprocess

from Bio import SeqIO

from .. import variant
from ..config import Config

log = logging.getLogger("easy45")

CONSENSUS_MAX_READS = 100   # reads sampled per cluster for the POA consensus


def _abpoa_consensus(seqs, tmp_fasta):
    """Run abpoa on a list of sequence strings; return the consensus string."""
    with open(tmp_fasta, "w") as f:
        for i, s in enumerate(seqs):
            f.write(f">{i}\n{s}\n")
    proc = subprocess.run(["abpoa", str(tmp_fasta)],
                          check=True, capture_output=True, text=True)
    return "".join(ln for ln in proc.stdout.splitlines() if not ln.startswith(">"))


def _cm_trim(cons, config):
    """Trim a consensus to the mature 18S 5' / 26S 3' boundary using the Rfam
    SSU+LSU covariance models (cmsearch). Returns (trimmed_seq, ssu_start, lsu_end)
    or (cons, None, None) if a boundary could not be located.
    """
    fa = config.workdir / "s5_cmtrim_input.fasta"
    tbl = config.workdir / "s5_cmtrim.tbl"
    fa.write_text(f">c\n{cons}\n")
    subprocess.run(["cmsearch", "--noali", "--cpu", str(config.threads),
                    "--tblout", str(tbl), str(config.cm_ref), str(fa)],
                   check=True, capture_output=True, text=True)
    ssu_from = lsu_to = None
    ssu_score = lsu_score = -1.0
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 15 or f[9] != "+":
                continue
            model, sfrom, sto, score = f[2], int(f[7]), int(f[8]), float(f[14])
            if "SSU" in model and score > ssu_score:
                ssu_score, ssu_from = score, min(sfrom, sto)
            elif "LSU" in model and score > lsu_score:
                lsu_score, lsu_to = score, max(sfrom, sto)
    if ssu_from is None or lsu_to is None or lsu_to <= ssu_from:
        return cons, None, None
    return cons[ssu_from - 1:lsu_to], ssu_from, lsu_to


def _members_by_cluster(membership_path):
    by_cluster: dict[int, list[str]] = {}
    with open(membership_path) as fh:
        next(fh)  # header
        for line in fh:
            rid, cl = line.rstrip("\n").split("\t")
            by_cluster.setdefault(int(cl), []).append(rid)
    return by_cluster


def _consensus_for(cluster, members, units, tmp, config):
    ids = members.get(cluster, [])[:CONSENSUS_MAX_READS]
    cons = _abpoa_consensus([str(units[i].seq) for i in ids], tmp)
    cons, s, _ = _cm_trim(cons, config)
    if s is None:
        log.warning("S5: CM boundary not found for cluster %d; keeping untrimmed", cluster)
    return cons


def run(config: Config, state: dict) -> dict:
    units = SeqIO.index(str(state["units_margin"]), "fasta")
    members = _members_by_cluster(state["membership"])

    rows = []
    with open(state["cluster_freq"]) as fh:
        next(fh)
        for line in fh:
            cl, _centroid, size, freq = line.rstrip("\n").split("\t")
            rows.append((int(cl), int(size), float(freq)))
    rows.sort(key=lambda r: r[1], reverse=True)

    consensus_fa = config.outdir / "consensus.fasta"
    variants_fa = config.outdir / "variants.fasta"
    ribotypes_tsv = config.outdir / "ribotypes.tsv"
    tmp = config.workdir / "s5_abpoa_input.fasta"
    primary_fa = config.workdir / "s5_primary.fasta"
    cand_fa = config.workdir / "s5_candidate.fasta"
    cmp_paf = config.workdir / "s5_variant_cmp.paf"

    primary_seq = None
    confirmed = []          # (cl, size, freq, cons, verdict)
    table = []              # (cl, size, freq, status, length, subs)

    for cl, size, freq in rows:
        if not (size >= config.min_reads and freq >= config.min_freq):
            table.append((cl, size, freq, "below-threshold", "", ""))
            continue
        cons = _consensus_for(cl, members, units, tmp, config)
        if primary_seq is None:
            primary_seq = cons
            primary_fa.write_text(f">primary\n{cons}\n")
            with open(consensus_fa, "w") as f:
                f.write(f">primary_ribotype reads={size} freq={freq:.4f} "
                        f"len={len(cons)}\n{cons}\n")
            table.append((cl, size, freq, "primary", len(cons), 0))
        else:
            cand_fa.write_text(f">cand\n{cons}\n")
            v = variant.confirm_variant(primary_fa, cand_fa, cmp_paf,
                                        config.min_variant_sites)
            if v.is_real:
                confirmed.append((cl, size, freq, cons, v))
            table.append((cl, size, freq, v.category, len(cons), v.substitutions))
    units.close()

    with open(variants_fa, "w") as f:
        for cl, size, freq, cons, v in confirmed:
            f.write(f">variant_cl{cl} reads={size} freq={freq:.4f} len={len(cons)} "
                    f"subs_vs_primary={v.substitutions}\n{cons}\n")

    with open(ribotypes_tsv, "w") as t:
        t.write("cluster\tsize\tfreq\tstatus\tconsensus_len\tsubs_vs_primary\n")
        for cl, size, freq, status, length, subs in table:
            t.write(f"{cl}\t{size}\t{freq:.4f}\t{status}\t{length}\t{subs}\n")

    log.info("S5: primary ribotype + %d confirmed variant(s) "
             "(homopolymer-aware); see %s", len(confirmed), ribotypes_tsv.name)
    return {"consensus": consensus_fa, "variants": variants_fa,
            "ribotypes": ribotypes_tsv}
