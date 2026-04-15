# Case Study: `max_dist` and `spanning_rescue` in the CCA3 Annotation

## Overview

Four Mender runs on the *Chamaeleo calyptratus* CCA3 annotation were compared
in a 2×2 design to evaluate how `max_dist` and `spanning_rescue` interact.
All runs used the same input GFF, protein data, IsoSeq GFF, and skip_flags
(`SKIPPED_GENE,LOW_COV`). Steps 1–5 were run once per `max_dist` value;
steps 6–8 were re-run per rescue setting.

| Run | `max_dist` | `spanning_rescue` | run prefix |
|-----|-----------|-------------------|------------|
| A | 4 | no | `20260414` |
| B | 4 | yes | `20260414_spanrescue` |
| C | 2 | no | `20260415_maxdist2` |
| D | 2 | yes | `20260415_maxdist2_spanrescue` |

`max_dist` controls how far apart two genes can be in the diamond tiling graph
before they are rejected as a pair. A distance of 1 means the two genes tile
adjacently on the reference protein; a distance of 4 means up to 3 other
genes can intervene in the tiling graph and the pair is still considered.
`spanning_rescue` rescues merges that would be excluded by `skip_flags` if
they have `FULL_SPAN` IsoSeq support — a single long-read transcript spanning
the entire locus end-to-end.

---

## Top-Level Results

|  | `spanning_rescue=no` | `spanning_rescue=yes` | rescue gain |
|---|---|---|---|
| **max_dist=4** | 992 PASS / 35 REVIEW | 1011 PASS / 36 REVIEW | **+19 PASS** |
| **max_dist=2** | 994 PASS / 35 REVIEW | 1009 PASS / 36 REVIEW | **+15 PASS** |

Two findings stand out immediately:

1. **`max_dist` alone barely matters.** Without rescue, max_dist=4 and max_dist=2
   differ by only 2 PASS merges (992 vs 994). The wider chaining window generates
   more candidates, but those extra candidates are filtered out by `SKIPPED_GENE`
   at roughly the same rate, leaving net output nearly identical.

2. **`spanning_rescue` is the larger lever.** Enabling it adds 15–19 PASS merges
   depending on `max_dist`, recovering IsoSeq-confirmed merges that the
   skip_flags gate would otherwise discard. The REVIEW count is essentially
   unchanged (35 → 36) — rescue is clean.

---

## Candidate Generation

Steps 1–5 (protein tiling, IsoSeq spanning) produce the candidate table.
Flag counts per run (a candidate may carry multiple flags):

| Flag | max_dist=4 | max_dist=2 |
|------|-----------|-----------|
| Candidates total | 1,123 | 1,118 |
| `CLEAN` | 537 | 540 |
| `SINGLE_HIT` | 434 | 434 |
| `STRONG` | 46 | 44 |
| `WEAK_END` | 31 | 28 |
| `TRANSITIVE_JOIN` | 5 | 2 |
| `LOW_COV` | 59 | 60 |
| `LARGE_SPAN` | 32 | 32 |
| **`SKIPPED_GENE`** | **42** | **34** |
| **`OPPOSITE_STRAND_SKIP`** | **18** | **16** |

The max_dist=4 run generates 8 more `SKIPPED_GENE` candidates than max_dist=2.
These are chains that can only form when the wider window lets a fragment leap
over intervening same-strand genes. Since `SKIPPED_GENE` is in `skip_flags`,
those 8 extra chains are excluded before merging — which is why max_dist=4
without rescue produces *fewer* PASS merges despite more candidates.

### The OPPOSITE_STRAND_SKIP correction

The strand-aware flag logic introduced in the current code distinguishes between:

- **`SKIPPED_GENE`** — at least one skipped gene is on the *same* strand as
  the chain. Genuine concern: that gene may be an additional split fragment.
- **`OPPOSITE_STRAND_SKIP`** — *all* skipped genes are on the opposite strand.
  In vertebrate genomes, genes on complementary strands routinely interleave —
  a + strand gene sitting in the intron of a − strand gene is commonplace. These
  are unrelated to the merge and do not go into `skip_flags`.

At max_dist=4, 18 of what would previously have been 60 `SKIPPED_GENE` flags
(30%) are now correctly classified as `OPPOSITE_STRAND_SKIP` and flow through
to merging cleanly. This correction prevents 18 false alarms per run.

---

## Effect of `spanning_rescue`

`spanning_rescue = yes` rescues merges blocked by `SKIPPED_GENE` or
`OPPOSITE_STRAND_SKIP` when `isoseq_flag` is `FULL_SPAN`. It does **not**
rescue `LOW_COV` or other evidence-quality flags — a spanning read confirms
co-transcription but does not substitute for missing protein tiling support.

### IsoSeq status of SKIPPED_GENE candidates

| | max_dist=4 | max_dist=2 |
|---|---|---|
| SKIPPED_GENE candidates | 42 | 34 |
| — FULL_SPAN | 22 | 18 |
| — PARTIAL_SPAN | 3 | 2 |
| — NO_SPANNERS | 17 | 14 |

All FULL_SPAN SKIPPED_GENE candidates are rescued when `spanning_rescue = yes`.
NO_SPANNERS and PARTIAL_SPAN cases remain excluded — absence of spanning reads
is not positive evidence for the merge, and PARTIAL_SPAN support does not reach
the terminal gene.

### Rescued merge outcomes

| | max_dist=4 rescue | max_dist=2 rescue |
|---|---|---|
| Total merges rescued | 20 | 16 |
| → PASS | 19 | 15 |
| → REVIEW | 1 | 1 |
| → FAIL | 0 | 0 |

No rescued merges fail translation validation. The single REVIEW at each
distance is merge_569: a 3-gene SKIPPED_GENE,STRONG chain with a borderline
second junction score (0.70) and merged protein coverage of 55.9% — just
below `min_merged_cov = 0.60`. This is a pre-existing borderline case, not
a rescue-specific problem.

### Overlap between the two rescue runs

Of the 20 merges rescued at max_dist=4 and 16 at max_dist=2:

- **16 rescued in both** — the core SKIPPED_GENE + FULL_SPAN pairs that exist
  at both distance settings
- **4 rescued only at max_dist=4** — all PASS; gene fragments too far apart at
  the tighter distance to form a chain at all
- **0 rescued only at max_dist=2**

---

## The 4 Merges Uniquely Recoverable at max_dist=4 + spanning_rescue

These are the loci that make the max_dist=4 + rescue combination the better
choice:

| Genes | Gene / protein | ref_cov | IsoSeq reads | Flags |
|-------|----------------|---------|--------------|-------|
| CCA3g003950, 003953, 003954 | HSPA4L (Hsp70 family) | 76.2% | FULL_SPAN ×7 | SKIPPED_GENE, STRONG |
| CCA3g008405, 008408 | ELP1 (Elongator complex) | 95.9% | FULL_SPAN ×13 | SKIPPED_GENE, SINGLE_HIT |
| CCA3g014094, 014097 | SLC19A2 (thiamine transporter) | 93.8% | FULL_SPAN ×21 | SKIPPED_GENE |
| CCA3g015952, 015955, 015956 | MCF2L (RhoGEF) | 71.1% | FULL_SPAN ×5 | SKIPPED_GENE, STRONG |

All four are PASS after translation validation. Two are 3-gene chains; without
max_dist=4 the fragments are too far apart in the tiling graph to form a chain
at all. At max_dist=2, none of the four gene sets appear as candidates.

**SLC19A2** (21 spanning reads, no WEAK flags, only blocked by two
opposite-strand interleaved genes in a gene-dense interval) and **ELP1**
(13 reads, single-copy gene) are the most straightforward rescues. **HSPA4L**
and **MCF2L** involve paralogous gene families — see the caveat below.

---

## IsoSeq Alignment Reliability Caveat

All rescue decisions rely on `FULL_SPAN` status, which is determined from
alignment coordinates only — no per-read MAPQ or alignment identity is captured
in the IsoSeq GFF. For gene families with high paralogy or repetitive protein
domains, a read spanning a locus may reflect a chimeric alignment across
paralogs rather than a genuine co-transcription event.

High-risk cases in this rescue set:

| Merge | Family | Risk |
|-------|--------|------|
| HSPA4L (CCA3g003950 cluster) | Hsp70 — large, conserved, multicopy | High |
| MCF2L (CCA3g015952 cluster) | RhoGEF family — paralogous pair MCF2L/MCF2L2 | Moderate |
| ELP1 (CCA3g008405/008408) | Single-copy elongation factor | Low |
| SLC19A2 (CCA3g014094/014097) | Single-copy solute carrier | Low |

For HSPA4L and MCF2L, verify spanning reads in a genome browser or re-derive
from the original BAM filtered at MAPQ ≥ 20 before accepting the rescue.
Requiring `isoseq_min_spanning >= 3` (all four cases exceed this) reduces but
does not eliminate risk from single chimeric reads.

---

## Biological Interpretation

### What these numbers actually mean

`max_dist` controls **how fragmented a split gene can be detected as**. At
max_dist=4, a gene split into 3+ pieces with gaps in the protein tiling graph
can still form a chain. But the near-identical PASS counts without rescue
(992 vs 994) reveal something important: **most split genes in CCA3 are simple
2-fragment splits** where the two pieces tile adjacently on the reference
protein. The wider window at max_dist=4 generates more multi-fragment
candidates, but those extra candidates almost all have a same-strand gene
intervening (`SKIPPED_GENE`) and are filtered — no net gain without rescue.

`spanning_rescue` solves a distinct biological problem: **the annotator
couldn't bridge a complex locus**. The pattern is: the annotator calls genes
A, B, C as three separate loci; B sits between two fragments of what should
be a single gene. Protein evidence alone cannot resolve this — you need
independent transcriptional evidence. A `FULL_SPAN` IsoSeq read spanning the
entire locus end-to-end, without a splicing gap that would suggest B is truly
independent, provides exactly that. The +19 clean PASS merges with essentially
no REVIEW inflation confirms that spanning reads are genuinely discriminatory
in this dataset, not noise.

### The 4 uniquely recoverable loci

These are the biologically significant gains at max_dist=4 + rescue:

| Gene | Biology | IsoSeq support | Confidence |
|------|---------|----------------|------------|
| **ELP1** | Single-copy elongator complex subunit; essential for neuronal development and tRNA wobble modification | 13 spanning reads | High — no paralogy risk |
| **SLC19A2** | Single-copy thiamine transporter (SLC19 family); haploinsufficiency causes thiamine-responsive megaloblastic anaemia | 21 spanning reads | High — highest support in the rescue set |
| **MCF2L** | RhoGEF signaling; has paralog MCF2L2 in the same genome | 5 spanning reads | Moderate — verify in genome browser |
| **HSPA4L** | Hsp70 chaperone family; highly conserved, multicopy, repetitive domain structure | 7 spanning reads | Lower — highest chimeric alignment risk |

ELP1 and SLC19A2 are single-copy, functionally important genes with strong
spanning-read support — these are almost certainly genuine split-gene errors.
HSPA4L and MCF2L sit in paralogous families where a spanning read could
reflect chimeric alignment across paralogs rather than genuine co-transcription;
they warrant manual verification before finalising the annotation.

### Why the IsoSeq data earns its weight here

These 4 loci share a profile that defeats automated annotators: large genomic
spans, gene-dense intervals, and complex multi-exon internal structure. In
each case the protein tiling evidence is strong (71–96% reference coverage)
but the annotator was deterred by an intervening same-strand gene call. For a
high-quality reference annotation like CCA3, leaving these loci fragmented
would be a genuine biological error. The IsoSeq data is not merely filtering
noise — it is resolving loci that are fundamentally ambiguous from protein
evidence alone.

### Parameter guidance by scenario

| Situation | max_dist | spanning_rescue | Rationale |
|-----------|----------|-----------------|-----------|
| IsoSeq available, diverse library | 4 | yes | Full recall; spanning reads provide independent co-transcription evidence |
| IsoSeq available, repetitive or polyploid genome | 4 | yes + review paralogous families | Chimeric alignment risk is higher; inspect `STRONG` flagged rescues manually; raise `isoseq_min_spanning` |
| No IsoSeq | 2 | no | max_dist=2 gains 2 PASS by avoiding SKIPPED_GENE chains that cannot be rescued anyway; `spanning_rescue` cannot fire |
| Targeting a specific gene family | 2 | no | Conservative setting avoids chaining unrelated paralogs |

---

## Recommendation

**`max_dist = 4` with `spanning_rescue = yes` gives the best recall.**

| Setting | PASS merges | notes |
|---------|-------------|-------|
| max_dist=2, no rescue | 994 | baseline |
| max_dist=4, no rescue | 992 | −2 vs baseline |
| max_dist=2, rescue | 1009 | +15 vs baseline |
| **max_dist=4, rescue** | **1011** | **+17 vs baseline; +4 biologically important loci** |

The 4 merges uniquely recovered at max_dist=4 + rescue are real multi-fragment
gene loci that cannot be assembled at max_dist=2. The rescue mechanism provides
the IsoSeq confirmation needed to accept them despite the intervening-gene flag.

**When to use max_dist=2 instead:** if IsoSeq data are unavailable and
`spanning_rescue` cannot be used, max_dist=2 gives 2 more PASS merges by
avoiding SKIPPED_GENE chains that will be discarded anyway. The difference is
negligible in practice; the more important factor is whether IsoSeq is
available at all.

### Recommended config for production runs with IsoSeq

```ini
max_dist        = 4
skip_flags      = SKIPPED_GENE,LOW_COV
spanning_rescue = yes
```

### Recommended config for runs without IsoSeq

```ini
max_dist        = 2
skip_flags      = SKIPPED_GENE,LOW_COV
spanning_rescue = no
```
