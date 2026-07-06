"""Batch mode — run the easy45 pipeline over every sample in a folder.

Auto-detects two input layouts (may be mixed):
  (a) flat       — one HiFi read file per sample directly in the folder
                   (``.fastq``/``.fq``/``.fasta``/``.fa``, optionally gzipped);
                   the sample name is the file name (minus a trailing
                   ``_hifi``/``_ccs``/``_reads`` and the extension).
  (b) subfolder  — one directory per sample holding its HiFi read file; the
                   sample name is the directory name.

Each sample -> ``outdir/<sample>/``. One sample failing does not stop the batch
(the error is logged and saved). Samples whose consensus already exists are
skipped (resume). ``batch_summary.tsv`` aggregates every sample + timing + status.
"""

from __future__ import annotations

import logging
import re
import time
import traceback
from pathlib import Path

from .config import Config, DEFAULT_ANCHOR
from .pipeline import run_pipeline

log = logging.getLogger("easy45")

_RD = re.compile(r"\.(fastq|fq|fasta|fa)(\.gz)?$", re.I)


def _sample_name(fname: str) -> str:
    base = _RD.sub("", fname)
    return re.sub(r"[._](hifi|ccs|reads)$", "", base, flags=re.I)


_FQ = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)


def _reads_in(folder: Path):
    return [p for p in sorted(folder.iterdir()) if p.is_file() and _RD.search(p.name)]


def _pick_read(rs):
    """From several read files in one sample folder, pick the best HiFi candidate:
    a HiFi/CCS-named file first, then FASTQ over FASTA (FASTQ keeps qualities)."""
    hifi = [r for r in rs if re.search(r"hifi|ccs", r.name, re.I)]
    pool = hifi or rs
    fq = [r for r in pool if _FQ.search(r.name)]
    return (fq or pool)[0]


def discover_samples(indir: Path):
    """Return [(sample_name, reads_file)] auto-detecting flat + subfolder layouts."""
    samples = [(_sample_name(f.name), f) for f in _reads_in(indir)]     # flat
    for sub in sorted(p for p in indir.iterdir() if p.is_dir()):
        rs = _reads_in(sub)
        if rs:                                                          # subfolder:
            samples.append((sub.name, _pick_read(rs)))                  # prefer HiFi/FASTQ
    seen, uniq = set(), []
    for s in samples:
        if s[0] not in seen:
            seen.add(s[0])
            uniq.append(s)
    return uniq


def _flen(fa: Path) -> str:
    """Total residues in a FASTA (pure-Python; blank if missing)."""
    if not fa.exists():
        return ""
    n = 0
    for line in fa.read_text().splitlines():
        if not line.startswith(">"):
            n += len(line.strip())
    return str(n) if n else ""


def _count_fa(fa: Path) -> str:
    if not fa.exists():
        return ""
    return str(sum(1 for line in fa.read_text().splitlines() if line.startswith(">")))


def run_batch(args) -> int:
    indir, outroot = Path(args.indir), Path(args.outdir)
    outroot.mkdir(parents=True, exist_ok=True)
    samples = discover_samples(indir)
    if not samples:
        log.error("batch: no HiFi read files found under %s", indir)
        return 1
    resume = not args.no_resume
    log.info("batch: %d sample(s) found under %s", len(samples), indir)

    rows = []
    for i, (name, reads) in enumerate(samples, 1):
        out = outroot / name
        if resume and (out / "consensus.fasta").exists():
            log.info("[%d/%d] %s — already done, skipping", i, len(samples), name)
            rows.append((name, out, "skipped(done)", 0))
            continue
        log.info("[%d/%d] %s  (%s)", i, len(samples), name, reads.name)
        t0 = time.time()
        try:
            cfg = Config(
                reads=reads, anchor_ref=args.anchor_ref or DEFAULT_ANCHOR,
                organelle_ref=args.organelle_ref, outdir=out, threads=args.threads,
                recruit_min_matchlen=args.recruit_min_matchlen, cluster_id=args.cluster_id,
                min_reads=args.min_reads, min_freq=args.min_freq,
                recover_igs=not args.no_igs, resume=resume,
            )
            run_pipeline(cfg)
            status = "OK"
        except Exception as e:  # one sample must not sink the batch
            status = f"FAIL:{type(e).__name__}"
            out.mkdir(parents=True, exist_ok=True)
            (out / "batch_error.txt").write_text(traceback.format_exc())
            log.error("[%d/%d] %s FAILED: %s (see batch_error.txt)", i, len(samples), name, e)
        rows.append((name, out, status, int(time.time() - t0)))

    _write_summary(outroot, rows)
    ok = sum(1 for r in rows if r[2] == "OK")
    log.info("batch: done — %d OK, %d skipped, %d failed of %d", ok,
             sum(1 for r in rows if r[2].startswith("skipped")),
             sum(1 for r in rows if r[2].startswith("FAIL")), len(rows))
    return 0


def _write_summary(outroot: Path, rows):
    out = outroot / "batch_summary.tsv"
    hdr = ["sample", "consensus_bp", "n_ribotypes", "ITS_bp", "IGS_bp", "seconds", "status"]
    with open(out, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for name, sdir, status, secs in rows:
            sdir = Path(sdir)
            fh.write("\t".join([
                name, _flen(sdir / "consensus.fasta"), _count_fa(sdir / "variants.fasta"),
                _flen(sdir / "its.fasta"), _flen(sdir / "igs.fasta"),
                str(secs), status]) + "\n")
    log.info("batch: wrote %s", out)
