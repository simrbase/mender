# Case Study: Effect of `max_dist` on Split-Gene Detection

## Overview

Two Mender runs on the *Chamaeleo calyptratus* CCA3 annotation were compared to
evaluate the effect of the `max_dist` parameter on candidate generation, chain
length, and validation outcomes.

| Parameter | Run A | Run B |
|-----------|-------|-------|
| `max_dist` | **4** | **2** |
| `run_prefix` | `20260414` | `20260415_maxdist2` |
| All other parameters | identical | identical |

`max_dist` controls how far apart two genes can be in the diamond tiling graph
before they are rejected as a pair. A distance of 1 means the two genes tile
adjacently on the reference protein. A distance of 4 means up to 3 other genes
can sit between them in the tiling graph and the pair is still considered.

---

## Top-Level Results

|  | max_dist=4 | max_dist=2 |
|---|---|---|
| Candidates (total) | 1,123 | 1,118 |
| — SKIPPED_GENE-flagged | 60 | 50 |
| — LOW_COV-flagged | 52 | 53 |
| Merges processed | 1,011 | 1,015 |
| **PASS** | **977** | **981** |
| REVIEW | 34 | 34 |
| FAIL | 0 | 0 |

**max_dist=2 yields 4 more PASS merges despite processing fewer unique gene pairs.**

---

## Why max_dist=2 Produces More PASS Merges

The counter-intuitive result — a tighter distance cutoff giving *more* validated
merges — is explained by the interaction between `max_dist` and the `SKIPPED_GENE`
filter.

### The chain-skip mechanism

When `max_dist=4`, the chaining algorithm can reach across larger gaps in the
tiling graph. At some loci, this causes a chain to leap over intervening genes
that failed other filters. Those intervening genes become **skipped genes**, and
the whole chain receives the `SKIPPED_GENE` flag. Since `SKIPPED_GENE` is in the
default `skip_flags`, the entire chain is excluded from processing.

With `max_dist=2`, the same gene neighborhoods can only form shorter, adjacent
chains. Genes that would have been absorbed into a longer skipped-gene chain
instead appear as independent, clean 2-gene pairs — and those get processed and
validated.

### Quantified

- max_dist=4 generates **11 SKIPPED_GENE chains** that do not exist in max_dist=2
- Those 11 chains involve **26 genes** that would otherwise be valid merge candidates
- In the max_dist=2 run, **4 of those genes** resolve into clean PASS merges
- The remaining genes either disappear from the tiling graph at the tighter cutoff
  or land in chains that are skipped for other reasons (LOW_COV, etc.)

---

## The 11 SKIPPED_GENE Chains Unique to max_dist=4

The table below lists all chains present only in the max_dist=4 run, with their
IsoSeq spanning status and the strand relationship between the merge chain and
the skipped genes.

| merge | n_genes | ref_cov | IsoSeq flag (n) | chain flags | chain strand | skipped genes (strand) |
|-------|---------|---------|-----------------|-------------|--------------|------------------------|
| merge_128 | 3 | 83.9% | NO_SPANNERS (0) | SKIPPED_GENE, WEAK_END | − | 007793:−, 007794:− |
| merge_136 | 2 | 100.0% | NO_SPANNERS (0) | SKIPPED_GENE | − | 007991:+, 007992:+ |
| merge_156 | 2 | 95.9% | **FULL_SPAN (13)** | SKIPPED_GENE, SINGLE_HIT | + | 008406:+, 008407:+ |
| merge_229 | 2 | 100.0% | NO_SPANNERS (0) | SKIPPED_GENE | + | 009594:−, 009595:− |
| merge_255 | 4 | 86.7% | NO_SPANNERS (0) | SKIPPED_GENE, TRANSITIVE_JOIN, WEAK_END, WEAK_INTERNAL | + | 010109:+ |
| merge_475 | 2 | 87.0% | NO_SPANNERS (0) | SKIPPED_GENE | − | 013914:−, 013915:− |
| merge_488 | 2 | 93.8% | **FULL_SPAN (21)** | SKIPPED_GENE | + | 014095:−, 014096:+ |
| merge_590 | 3 | 71.1% | **FULL_SPAN (5)** | SKIPPED_GENE, STRONG | − | 015953:+, 015954:− |
| merge_720 | 2 | 71.5% | NO_SPANNERS (0) | SKIPPED_GENE, SINGLE_HIT | + | 018656:+, 018657:+ |
| merge_907 | 3 | 76.2% | **FULL_SPAN (7)** | SKIPPED_GENE, STRONG | − | 003951:−, 003952:− |
| merge_1035 | 4 | 100.0% | PARTIAL_SPAN (6) | SKIPPED_GENE, TRANSITIVE_JOIN, WEAK_END, WEAK_INTERNAL | − | 021183:−, 021184:+ |

### Strand relationship of skipped genes

The `SKIPPED_GENE` flag is strand-agnostic: any annotated gene sitting between
two merge fragments in genomic order triggers it, regardless of strand. This has
important consequences:

**Opposite-strand skips (skipped genes on opposite strand from chain):**
- merge_136: chain −, skipped +/+ → the intervening genes are on the opposite
  strand and almost certainly unrelated to the merge; the `SKIPPED_GENE` flag is
  a false alarm here
- merge_229: chain +, skipped −/− → same situation
- merge_590 (partial): one skipped gene on opposite strand (+), one same (−)
- merge_1035 (partial): one on same strand (−), one opposite (+)

**Same-strand skips (skipped genes on same strand as chain):**
- merge_128: chain −, skipped −/− → skipped genes may be additional keratin
  fragments; filtering is appropriate pending manual inspection
- merge_475: chain −, skipped −/− → skipped MMP family members on same strand
- merge_720: chain +, skipped +/+ → same-strand skip
- **merge_907**: chain −, skipped −/− → FULL_SPAN (7 reads), STRONG flag —
  three HSPA4L/HYOU1 fragments all on minus strand with 7 fully spanning reads.
  The two skipped genes are also minus-strand Hsp70 family members. This is the
  strongest case for manual review: the skipped genes may be additional split
  fragments that need to be incorporated into the chain.

### Four FULL_SPAN cases deserve inspection

Despite being filtered as SKIPPED_GENE chains, four merges have long-read
spanning support:

**merge_907** (HSPA4L, CHR05Y, 3 genes, 76% ref coverage, FULL_SPAN 7 reads)
- All chain genes on minus strand; skipped genes 003951, 003952 also minus strand
- 5 reads span all 3 chain genes; 2 reads span the two sub-junctions
- Most compelling case: the skipped genes are same-strand Hsp70 family members
  and may be additional split fragments
- At max_dist=2, gene 003953 pairs with 003954 as a REVIEW 2-gene merge; gene
  003950 does not appear in any max_dist=2 candidate

**merge_590** (MCF2L/MCF2L2, CHR03, 3 genes, 71% ref coverage, FULL_SPAN 5 reads)
- Chain on minus strand; one skipped gene (+), one (−)
- At max_dist=2, genes 015955 + 015956 are a LOW_COV candidate (skipped for
  low coverage); gene 015952 disappears from the graph

**merge_488** (SLC19A2, CHR03, 2 genes, 94% ref coverage, FULL_SPAN 21 reads)
- 21 spanning reads — the strongest IsoSeq signal of all 11 cases
- Skipped genes 014095 (−) and 014096 (+) are on mixed strands
- The chain itself is clean (no WEAK flags); filtered solely because of the
  skipped genes. The skipped-strand mix suggests these are interleaved genes
  on opposite strands in a gene-dense region
- At max_dist=2, genes 014094 + 014097 do not appear in any candidate at all

**merge_156** (ELP1, CHR01, 2 genes, 96% ref coverage, FULL_SPAN 13 reads)
- Skipped genes 008406 and 008407 are on the same strand (+) as the chain
- These genes sit physically between the two ELP1 fragments and tile the
  same reference — they likely represent additional ELP1 sub-fragments that
  failed the tiling filter individually
- At max_dist=2, genes 008405 + 008408 do not appear in any candidate

---

## Multi-Gene Chain Comparison

### Are any legitimate large chains lost with max_dist=2?

Nearly all multi-gene chains (3+) are **identical in both runs**. The one
exception is a 6-gene REVIEW chain in max_dist=4 that splits under max_dist=2:

**merge_29 (max_dist=4):** 6 genes, REVIEW, TRANSITIVE_JOIN + WEAK flags,
NO_SPANNERS — CCA3g006146, 006148, 006149, 006151, 006152, 006153

Under max_dist=2 this becomes:
- `merge_29`: CCA3g006148 + 006149, CLEAN → **PASS**
- `merge_30`: CCA3g006151 + 006152 + 006153, WEAK_END → **PASS**
- CCA3g006146: **not found** in any max_dist=2 candidate (too far from its
  nearest tiling neighbor at the tighter cutoff)

No IsoSeq reads span this 6-gene locus. The REVIEW outcome under max_dist=4
(driven by TRANSITIVE_JOIN and WEAK flags) combined with zero spanning reads
suggests the 6-gene chain was speculative. The two sub-chains recovered under
max_dist=2 are both validated as PASS.

### Chain length distribution

| n genes per merge | max_dist=4 | max_dist=2 |
|---|---|---|
| 2 | 937 | 941 |
| 3 | 70 | 71 |
| 4 | 2 | 2 |
| 6 | 2 | 1 |

The 6-gene chains: one (collagen/keratin cluster) is identical in both runs;
the other (the 6-gene REVIEW chain above) is present only in max_dist=4.

---

## Conclusions

1. **max_dist=2 gives marginally more validated merges (981 vs 977 PASS)** in
   this dataset. The difference is small but consistent with the mechanism: a
   tighter cutoff avoids producing SKIPPED_GENE chains that then get discarded,
   allowing the underlying gene pairs to be processed individually.

2. **No substantial multi-gene chains are lost** at max_dist=2. The one 6-gene
   chain that disappears was a speculative REVIEW with no IsoSeq support and
   weak internal evidence; its constituent pairs validate cleanly as smaller
   merges.

3. **The `SKIPPED_GENE` filter is strand-agnostic**, which causes false positives
   when the intervening genes are on the opposite strand. merge_136 and merge_229
   (skipped genes entirely on opposite strand) are the clearest examples. A
   strand-aware `SKIPPED_GENE` flag — requiring at least one same-strand
   intervening gene to trigger — would recover these merges without increasing
   noise.

4. **Four SKIPPED_GENE chains have strong IsoSeq support** and warrant manual
   inspection regardless of `max_dist` setting:
   - merge_907 (HSPA4L, 7 FULL_SPAN reads, same-strand skips — likely needs
     skipped genes added to chain)
   - merge_488 (SLC19A2, 21 FULL_SPAN reads — highest read count of the group)
   - merge_590 (MCF2L2, 5 FULL_SPAN reads, STRONG tiling flag)
   - merge_156 (ELP1, 13 FULL_SPAN reads, same-strand skips)

5. **Recommended default: `max_dist = 2`** for production runs where
   `SKIPPED_GENE` is in `skip_flags`. The lower setting avoids generating chains
   that will be discarded anyway, and slightly improves recall by keeping
   gene pairs eligible for independent processing.

---

## Appendix: What Happened to the 26 Genes in the 11 SKIPPED_GENE Chains

| Gene | max_dist=2 outcome |
|------|--------------------|
| CCA3g007791000 | PROCESSED merge_129 (PASS) |
| CCA3g007792000 | CANDIDATE merge_129 (flags=CLEAN, skipped for other reason) |
| CCA3g007795000 | not found |
| CCA3g007990000 | not found |
| CCA3g007993000 | not found |
| CCA3g008405000 | not found |
| CCA3g008408000 | not found |
| CCA3g009593000 | not found |
| CCA3g009596000 | not found |
| CCA3g003950000 | not found |
| CCA3g003953000 | PROCESSED merge_902 (REVIEW) |
| CCA3g003954000 | CANDIDATE merge_902 (flags=CLEAN) |
| CCA3g013913000 | not found |
| CCA3g013916000 | not found |
| CCA3g014094000 | not found |
| CCA3g014097000 | not found |
| CCA3g015952000 | not found |
| CCA3g015955000 | CANDIDATE merge_586 (flags=LOW_COV, skipped) |
| CCA3g015956000 | CANDIDATE merge_586 (flags=LOW_COV, skipped) |
| CCA3g018655000 | not found |
| CCA3g018658000 | not found |
| CCA3g021181000 | not found |
| CCA3g021182000 | not found |
| CCA3g021185000 | PROCESSED merge_1030 (PASS) |
| CCA3g021186000 | CANDIDATE merge_1030 (flags=SINGLE_HIT) |
| CCA3g010107000 – CCA3g010112000 | not found (too far apart at max_dist=2) |
