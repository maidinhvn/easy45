# Provenance of the boundary covariance models

File: `rrna_ssu_lsu.cm`
Date obtained: 2026-06-29 (Rfam, https://rfam.org)

## What it is

A concatenation of two Rfam covariance models used by easy45 (S5) to trim the
primary/variant consensus to the precise mature 18S 5' and 26S 3' boundaries
with `cmsearch` (Infernal):

| Model | Rfam accession | Region |
|---|---|---|
| SSU_rRNA_eukarya | RF01960 | eukaryotic small-subunit (18S) rRNA |
| LSU_rRNA_eukarya | RF02543 | eukaryotic large-subunit (26S/28S) rRNA |

## Why CMs (not barrnap) define the final boundary

barrnap is a fast HMM gene *finder* (used in S2 to detect spanning reads), but
it clips a few-to-tens of bp off the 18S 5' and 26S 3' termini. Covariance
models are structure-aware and trained on broad multi-taxon Rfam seed
alignments, so they delineate the mature rRNA termini consistently and
generally — independent of any single reference sequence. Validated on a
Polyscias lallanii rDNA: cmsearch recovered the 18S 5' base and extended the
26S 3' end relative to barrnap.

## How it was obtained

    curl -s https://rfam.org/family/RF01960/cm -o RF01960.cm
    curl -s https://rfam.org/family/RF02543/cm -o RF02543.cm
    cat RF01960.cm RF02543.cm > rrna_ssu_lsu.cm

## Citation (for the paper)

> Kalvari I, et al. Rfam 14: expanded coverage of metagenomic, viral and
> microRNA families. Nucleic Acids Res. 2021;49(D1):D192-D200.

Infernal: Nawrocki EP, Eddy SR. Infernal 1.1: 100-fold faster RNA homology
searches. Bioinformatics. 2013;29(22):2933-2935.
