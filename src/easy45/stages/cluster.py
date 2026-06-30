"""S4 — coarse clustering of per-read units.

vsearch (--cluster_size --id config.cluster_id) groups near-identical units.
This is intentionally a *coarse* grouping: it collapses obvious duplicates and
separates clearly divergent ribotypes, but the authoritative variant decision
is NOT made here. Fine-grained, homopolymer-aware variant calling happens in
``variant.py`` (driven by S5), because a single homopolymer indel can drop
percent-identity as much as real divergence.

Units are already strand-normalised (S3), so clustering runs --strand plus.

Input keys:  {"units": <FASTA>}
Output keys: {"centroids": <FASTA>, "uc": <vsearch .uc>,
              "membership": <TSV read_id->cluster>, "cluster_freq": <TSV>}
"""

from __future__ import annotations

import logging

from ..config import Config
from ..external import run as sh

log = logging.getLogger("easy45")


def _parse_uc(uc_path):
    """Return (membership {read_id: cluster}, sizes {cluster: n})."""
    membership: dict[str, int] = {}
    sizes: dict[int, int] = {}
    with open(uc_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rec = f[0]
            if rec in ("S", "H"):
                cl = int(f[1])
                read_id = f[8].split()[0]          # label, drop description
                membership[read_id] = cl
            elif rec == "C":
                sizes[int(f[1])] = int(f[2])
    return membership, sizes


def run(config: Config, state: dict) -> dict:
    units = state["units"]
    centroids = config.workdir / "s4_centroids.fasta"
    uc = config.workdir / "s4_clusters.uc"
    membership_tsv = config.workdir / "s4_membership.tsv"
    freq_tsv = config.workdir / "s4_cluster_freq.tsv"

    log.info("S4: clustering units with vsearch (--id %.3f)", config.cluster_id)
    sh(["vsearch", "--cluster_size", str(units),
         "--id", str(config.cluster_id),
         "--strand", "plus", "--sizeout",
         "--threads", str(config.threads),
         "--centroids", str(centroids), "--uc", str(uc)])

    membership, sizes = _parse_uc(uc)
    total = sum(sizes.values()) or 1

    # cluster -> representative read id (centroid = the 'S' record for that cluster)
    centroid_of: dict[int, str] = {}
    with open(uc) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                centroid_of[int(f[1])] = f[8].split()[0]

    with open(membership_tsv, "w") as m:
        m.write("read_id\tcluster\n")
        for rid, cl in membership.items():
            m.write(f"{rid}\t{cl}\n")

    ordered = sorted(sizes.items(), key=lambda kv: kv[1], reverse=True)
    with open(freq_tsv, "w") as t:
        t.write("cluster\tcentroid_read\tsize\tfreq\n")
        for cl, n in ordered:
            t.write(f"{cl}\t{centroid_of.get(cl,'?')}\t{n}\t{n/total:.4f}\n")

    top = ordered[0] if ordered else (None, 0)
    log.info("S4: %d units -> %d clusters; largest cluster = %d reads (%.1f%%)",
             total, len(sizes), top[1], 100 * top[1] / total)
    return {"centroids": centroids, "uc": uc,
            "membership": membership_tsv, "cluster_freq": freq_tsv}
