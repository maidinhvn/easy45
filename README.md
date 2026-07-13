# easy45

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/easy45)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Assembly-free recovery of the 45S nrDNA transcribed unit and ribotype variants from HiFi long reads.**

easy45 recruits HiFi reads that span a full 45S nrDNA transcribed unit
(ETS–18S–ITS1–5.8S–ITS2–26S), then reconstructs that unit *without genome
assembly*. Because each HiFi read can span an entire repeat unit, easy45 treats
every spanning read as one independent molecule — letting it report not just a
consensus, but the genuine intragenomic ribotype diversity that assembly-based
consensus collapses.

## Outputs

1. **Consensus** — the primary ribotype of the transcribed unit.
2. **Variants** — additional ribotypes that pass homopolymer-aware, read-supported
   validation (a signal of hybridisation / allopolyploidy).
3. **IGS** *(best-effort)* — the intergenic spacer, reported per length class
   when reads span it cleanly.

## Pipeline

![easy45 pipeline](docs/pipeline.png)

Recruit (minimap2) → locate genes (barrnap) → excise & orient units → cluster
(vsearch) → consensus + mature-boundary trim (abpoa + Rfam CMs) + homopolymer-aware
ribotype calling → IGS → annotate + ITS barcode (ITSx). See
[docs/pipeline.md](docs/pipeline.md) for the editable flowchart and details.

## Install

### From bioconda (recommended)

```bash
conda install -c bioconda easy45
# or, faster:
mamba install -c bioconda easy45
```

That's it — all external tools (minimap2, seqkit, vsearch, ITSx, barrnap, abpoa,
infernal) are pulled in automatically; the Python package itself stays pure-Python.

### From source (development)

```bash
conda env create -f environment.yml
conda activate easy45
easy45 check-deps
```

## Usage

The 45S anchor is bundled, so the only thing you pass is your HiFi reads. Results
are written to `easy45_out/` (created automatically):

```bash
easy45 run -i reads.fastq.gz
```

Add options as needed — a different output folder (`-o`), organelle depletion
(`-r`), more threads (`-t`):

```bash
easy45 run -i reads.fastq.gz -o results/ -r plastid_mito.fasta -t 16
```

| Flag | Required? | Meaning |
|------|-----------|---------|
| `-i, --reads` | **yes** | your HiFi reads (FASTA/FASTQ, optionally gzipped) |
| `-o, --outdir` | no | output folder — created automatically (default `easy45_out/`) |
| `-r, --organelle-ref` | no | plastid+mitochondrial genomes for organelle depletion (recommended) |
| `-a, --anchor-ref` | no | custom 45S anchor (defaults to the bundled Arabidopsis T2T unit) |
| `-t, --threads` | no | threads (default 4) |

Run `easy45 run --help` for all parameters.

### Batch mode

Process a whole folder of HiFi samples in one command — ideal for an overnight run:

```bash
easy45 batch -i hifi_folder/ -o out/ -t 16
```

`batch` auto-detects either input layout (you can even mix them in one folder):

```
# flat — sample name = file name        # subfolder — sample name = folder name
hifi_folder/                            hifi_folder/
├── sampleA.hifi.fastq.gz               ├── sampleA/
├── sampleB.fastq.gz                    │   └── reads.hifi.fastq.gz
└── sampleC.fasta                       └── sampleB/
                                            └── m64.ccs.fastq.gz
```

(trailing `_hifi`/`_ccs`/`_reads` is stripped from flat file names; when a subfolder
holds several files a HiFi/CCS-named FASTQ is preferred.)

Each sample is written to `out/<sample>/`. The run is robust for long jobs: a failed
sample is logged and skipped (the batch keeps going), and samples already finished
are skipped on re-launch (resume — pass `--no-resume` to force a redo). A top-level
**`batch_summary.tsv`** aggregates every sample's consensus length, ribotype/variant
count, ITS & IGS lengths, timing and status. All `run` options (`-r`, `--cluster-id`,
…) apply.

## Citation

If you use easy45, please cite the manuscript (in preparation) and the archived
release on Zenodo (DOI to be added on first release).

## License

MIT. See [LICENSE](LICENSE).

## Status

Functional. All pipeline stages (S0–S7) are implemented and unit-tested, and the
recovered consensus has been validated as base-identical to independent
whole-genome assemblies. A manuscript describing easy45 is in preparation; a
versioned Zenodo release and bioconda package will accompany it.
