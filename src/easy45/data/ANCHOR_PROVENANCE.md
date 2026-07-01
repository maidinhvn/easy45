# Provenance of the default recruitment anchor

File: `default_anchor_arabidopsis_45S.fasta`
Date obtained: 2026-06-29 (NCBI E-utilities)
Reproduce with: `python scripts/make_default_anchor.py`

## What it is

One **intact, complete 45S nrDNA repeat unit** of *Arabidopsis thaliana*
(10,052 bp), spanning 18S–ITS1–5.8S–ITS2–26S–IGS, used as the default bait for
HiFi read recruitment (Stage 1). The anchor is **external to and independent of
any sample**: it does not require assembling the query genome, so it does not
compromise the assembly-free design of easy45.

## Source

| Field | Value |
|---|---|
| Organism | *Arabidopsis thaliana* (taxon:3702) |
| Isolate | CS1092 |
| Assembly region | Nucleolus organizing region **NOR2** |
| GenBank accession | **OR453402.1** |
| Extracted coordinates | `OR453402.1:36450-46501` (one repeat unit) |
| Repeat period in NOR2 | ~10,051 bp |
| Length of extracted unit | 10,052 bp |

### Primary citation (cite this in the paper)

> Fultz D, McKinlay A, Enganti R, Pikaard CS. Sequence and epigenetic
> landscapes of active and silent nucleolus organizing regions in
> *Arabidopsis*. *Sci Adv*. 2023;9(44):eadj4509. PMID: 37910609.

Direct submission: Fultz D., Howard Hughes Medical Institute / Indiana
University (submitted 14-AUG-2023).

This is the T2T (telomere-to-telomere) reference of the *Arabidopsis* NOR — the
most complete and modern rDNA reference available, chosen over the classic but
fragmented 1990 deposits (18S in X16077; 5.8S + 25S in X52320, whose 18S is only
a partial `<1..228` fragment).

## How it was obtained (methods, for the paper)

1. Queried NCBI nucleotide for complete *A. thaliana* rDNA references; identified
   the T2T NOR2 (OR453402.1, 5.42 Mb) and NOR4 (OR453401.1, 3.85 Mb) arrays.
2. Fetched the first 60 kb of OR453402.1 via NCBI E-utilities `efetch`
   (`seq_start=1, seq_stop=60000`), spanning several tandem repeat units.
3. Located repeat-unit boundaries by searching for the 5′ 40-mer of the
   *A. thaliana* 18S rRNA gene (from X16077.1, gene 88..1891). Exact matches
   occurred at NOR2 positions 4719, 36449, 46501, 56552. The first hit (4719) is
   anomalous (atypical array-edge copy); the canonical repeat period of
   ~10,051 bp is defined by the consistent later hits.
4. Excised one complete unit between two consecutive canonical 18S start
   positions (0-based slice `36449:46501`; 1-based `OR453402.1:36450-46501`),
   giving a 10,052 bp unit.
5. Verified the unit contains intact 18S, 5.8S and 25S/26S genes (seed-match
   against X16077 and X52320 gene regions).

## Caveat for downstream design

The edges of an rDNA array can contain non-canonical / variant repeat units
(cf. the anomalous 18S match at NOR2 position 4719). Unit-boundary logic
(Stages 2/6) should not assume every array position yields a textbook unit.
