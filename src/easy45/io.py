"""I/O helpers: sequence reading/writing and GFF output.

Sequence parsing uses Biopython (the only non-trivial pip dependency). Kept in
one module so format handling (gzip, FASTA/FASTQ autodetect) lives in one place.
"""

from __future__ import annotations

import gzip
from pathlib import Path


def open_maybe_gzip(path: Path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def detect_format(path: Path) -> str:
    """Return 'fasta' or 'fastq' by sniffing the first record character."""
    with open_maybe_gzip(path) as fh:
        first = fh.read(1)
    if first == ">":
        return "fasta"
    if first == "@":
        return "fastq"
    raise ValueError(f"Unrecognised sequence format: {path}")


def write_gff(records, path: Path) -> None:
    """Write feature records (18S/ITS1/5.8S/ITS2/26S) as GFF3. Phase 2."""
    raise NotImplementedError
