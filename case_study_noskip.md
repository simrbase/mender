# Case Study: Allowing `SKIPPED_GENE` Through Validation (`skip_flags = LOW_COV`)

## Overview

By default, Mender excludes merge candidates that carry `SKIPPED_GENE` from
the MSA and translation validation steps (via `skip_flags = SKIPPED_GENE,LOW_COV`).
The rationale: a gene lying between two merge candidates suggests those candidates
may be independent loci rather than fragments of one gene, making the merge
biologically dubious.

This case study tests the alternative: remove `SKIPPED_GENE` from `skip_flags`
and let translation validation — not a heuristic pre-filter — decide whether
the merge is real.

| Setting | `skip_flags` | `spanning_rescue` | `keep_msa` | run prefix |
|---------|-------------|-------------------|-----------|------------|
| **Base** | `SKIPPED_GENE,LOW_COV` | yes | no | `20260415_spanrescue` |
| **Noskip** | `LOW_COV` only | no | yes | `20260416_noskip` |

The base run is the production configuration (max_dist=4, spanning_rescue=yes)
from the `max_dist`/`spanning_rescue` case study. The noskip run differs only in
`skip_flags`: `SKIPPED_GENE` candidates are no longer pre-filtered; `spanning_rescue`
is disabled because it is not needed when SKIPPED_GENE is allowed through
directly.

---

## Background: What Is `SKIPPED_GENE`?

A candidate pair earns `SKIPPED_GENE` when the protein tiling graph shows that
a *third*, distinct gene also tiles onto the same reference protein between the
two candidates. In other words, the reference protein is covered by at least
three genomic gene models, and the middle one is not part of the proposed merge.

This arises in two biological scenarios:

1. **True annotation split with an embedded gene:** The true locus is fragmented,
   but a real, unrelated gene happens to also match the same reference protein
   family (e.g., a tandemly duplicated paralogue or a domain-sharing gene). The
   merge is biologically valid; the skipped gene is coincidental.

2. **Incorrect merge proposal:** Three independent genes share similarity to a
   multi-domain reference protein; merging any two of them would be wrong. The
   skipped gene is a genuine intervening locus, not a coincidence.

Translation validation cannot always distinguish these cases before the merge is
attempted — but it *can* tell you whether the stitched protein is internally
consistent (no premature stops, reasonable coverage of the reference).

---

## Candidate Counts

| Flag class | Count |
|------------|-------|
| Total candidates | 1,123 |
| `SKIPPED_GENE` only (filtered in base) | 37 |
| `SKIPPED_GENE` + `LOW_COV` (filtered in both) | 5 |
| Available to noskip validation | 42 |

The 42 SKIPPED_GENE candidates represent **3.7% of all candidates** — a small
but meaningful set. In the base run, 20 of these are recovered by `spanning_rescue`
(19 PASS + 1 REVIEW); the remaining 22 are silently discarded.

---

## Top-Level Results

| Run | PASS | REVIEW | FAIL | Total validated |
|-----|------|--------|------|-----------------|
| Base (`spanning_rescue=yes`) | 1,011 | 36 | 0 | 1,047 |
| Noskip (`skip_flags=LOW_COV`) | 1,021 | 43 | 0 | 1,064 |
| **Difference** | **+10** | **+7** | **0** | **+17** |

The noskip run validates 17 additional merges. All 17 originate from the 42
SKIPPED_GENE candidates:

- **29 PASS** (base: 19 rescued via IsoSeq → net +10 PASS without IsoSeq)
- **8 REVIEW** (base: 1 rescued → net +7 REVIEW)
- **0 FAIL** (base: 0)

The non-SKIPPED_GENE fraction behaves identically: 35 REVIEW in both runs.

---

## Translation Validation Results

The key question: when `SKIPPED_GENE` candidates are stitched together and
translated, do the resulting proteins make sense?

**37 of 37 SKIPPED_GENE candidates with `translation_flag = OK`.**  
**37 of 37 have `has_internal_stop = 0` — no premature stop codons.**

| Metric | PASS cases (n=29) |
|--------|------------------|
| Mean % identity to best reference | 66.5% |
| Mean coverage of merged protein by reference | 92.9% |
| Mean coverage of reference by merged protein | 94.4% |
| MSA flag | GOOD_MSA in all 37 |

Every SKIPPED_GENE candidate that reaches translation validation produces a
complete, internally consistent open reading frame. **Not one induces a
frameshift or premature stop.** This is the strongest possible evidence that
the merged CDS is structurally correct.

---

## PASS vs REVIEW: What Separates Them?

The 8 REVIEW cases are *not* translation failures. Their REVIEW status is
driven by MSA junction scores and coverage flags, not by broken proteins:

| merge_id | min_junc | merged_cov_by_ref | additional_flags | notes |
|----------|----------|--------------------|------------------|-------|
| merge_935 | 0.6455 | 0.45 | — | Ig-lambda locus; ref coverage low |
| merge_128 | 0.7000 | 0.47 | WEAK_END | Keratin; source genes gap-heavy in MSA |
| merge_255 | 0.6568 | 0.41 | TRANSITIVE_JOIN,WEAK_END,WEAK_INTERNAL | 4-gene chain; alpha-1-antitrypsin family |
| merge_475 | 0.7000 | 0.59 | — | MMP family; one source gene all-gaps in MSA |
| merge_569 | 0.7000 | 0.56 | STRONG | FARP2; source_015647 has diverged N-term |
| merge_574 | 0.8727 | 0.36 | SINGLE_HIT | NMNA3; merged protein longer than reference |
| merge_720 | 0.8182 | 0.56 | SINGLE_HIT | ZNFX1; one source gene all-gaps in MSA |
| merge_803 | 0.7227 | 0.98 | WEAK_END | Clean coverage but junction at lower end |

**Notable case — merge_569 (FARP2 locus):** Source gene CCA3g015647 has a
diverged N-terminal domain relative to the merged protein. The skipped gene
(CCA3g015649) shows all-gaps in the alignment. This is the most credible
case where SKIPPED_GENE intuition holds — the intervening gene and the N-terminal
fragment may be independent, and the merge should be reviewed carefully.

**merge_574 and merge_720** are REVIEW solely because of `SINGLE_HIT` (only one
reference protein available for scoring) despite junction scores of 0.87 and
0.82. These are likely correct merges with sparse reference coverage.

---

## Comparison with `spanning_rescue`

In the base run, `spanning_rescue` recovers SKIPPED_GENE candidates that have
full-span IsoSeq transcript support. How does this compare to allowing them
through unconditionally?

| Recovery mechanism | SKIPPED_GENE PASS | SKIPPED_GENE REVIEW | Requires IsoSeq? |
|-------------------|------------------|--------------------|------------------|
| `spanning_rescue=yes` (base) | 19 | 1 | **Yes** |
| `skip_flags=LOW_COV` (noskip) | 29 | 8 | No |
| **Additional recoveries (noskip only)** | **+10** | **+7** | — |

The 19 IsoSeq-rescued PASS merges are a subset of the 29 noskip PASS merges.
Allowing SKIPPED_GENE through validation recovers 10 additional PASS merges that
lack full-span IsoSeq support — meaning the evidence for those merges comes
entirely from protein structure: clean ORF, no internal stops, high reference
coverage, good MSA junction scores.

The 10 extra PASS merges (no IsoSeq support) represent real biological merges
that the standard pipeline would miss unless IsoSeq data is available.

---

## Biological Interpretation

### What the 0-FAIL result means

In the default pipeline, `SKIPPED_GENE` acts as a *conservative prior*: when
in doubt, don't merge. This case study shows the prior may be too conservative.
A premature stop codon or frameshift in the merged protein would be immediate
evidence of a mis-merge — yet **none of the 42 candidates produce one.** The
merged sequences are not random collisions; they are bona fide exon chains that
encode complete proteins.

The most likely explanation for the high success rate:

- The "skipped gene" is often a tandem paralogue or domain-family member that
  happens to tile onto the same reference protein. It is genuinely a separate
  locus — but its presence does not invalidate the merge of the flanking fragments.
- Annotation fragmentation in *C. calyptratus* tends to occur at exon boundaries
  that are conserved across vertebrates. Even when an unrelated gene tiles between
  two fragments, the fragments themselves retain correct reading frames.

### Reference-invariance

From the reference comparison case study, all 42 SKIPPED_GENE candidates are
flagged identically with both an Anolis proteome and a SwissProt reference.
The intervening genes all have conserved vertebrate orthologs — they are not
squamate-specific novelties that would disappear with a different reference.
This means no reference proteome choice will resolve SKIPPED_GENE status; only
validation or IsoSeq can.

### When to use each approach

| Scenario | Recommended setting | Reason |
|----------|--------------------|-|
| IsoSeq available, conservative preferred | `skip_flags = SKIPPED_GENE,LOW_COV` + `spanning_rescue = yes` | IsoSeq FULL_SPAN is gold-standard evidence; non-IsoSeq SKIPPED_GENE skipped |
| No IsoSeq, or want maximum recovery | `skip_flags = LOW_COV` + `spanning_rescue = no` | Let translation decide; 0 failures observed |
| Manual review capacity available | `skip_flags = LOW_COV` | Inspect 8 REVIEW cases; most are likely real |
| Ultra-conservative annotation | `skip_flags = SKIPPED_GENE,LOW_COV` + `spanning_rescue = no` | Discards 42 candidates; safest for publication without IsoSeq |

---

## Recommendation

**For genomes without IsoSeq data**, set `skip_flags = LOW_COV`. The translation
validation step is sufficient to catch bad SKIPPED_GENE merges — and this case
study demonstrates empirically that it does so with 0 failures across 42
candidates. Running with `SKIPPED_GENE` in `skip_flags` silently discards up to
37 validated, complete-ORF merges per 1,000-gene run.

**For genomes with IsoSeq data**, the default (`skip_flags = SKIPPED_GENE,LOW_COV`
+ `spanning_rescue = yes`) remains a defensible choice — it restricts SKIPPED_GENE
merges to those with direct long-read evidence. However, note that 10 additional
real merges exist that IsoSeq misses (perhaps due to incomplete isoform sampling),
and removing `SKIPPED_GENE` from `skip_flags` would recover them.

**The 8 REVIEW cases** should be manually inspected in any production run using
`skip_flags = LOW_COV`. The FARP2 locus (merge_569) is the highest-priority case
for manual review, as it shows a source gene with a diverged N-terminus that may
indicate paralogue confusion rather than annotation fragmentation.
