# Mender — Pipeline Overview

Mender detects erroneously split gene models in a genome annotation by
cross-referencing protein homology against a reference proteome, with
optional support from PacBio IsoSeq long-read transcripts. Adjacent gene
fragments that together tile a single reference ortholog — and are
optionally confirmed by a spanning long-read transcript — are merged into
one corrected gene model.

Box style legend:

```
  ╔══════════╗   runs in every execution
  ╚══════════╝

  ┌╌╌╌╌╌╌╌╌╌╌┐   conditional — trigger shown in the step label
  └╌╌╌╌╌╌╌╌╌╌┘
```

---

```
╔══════════════════════════════════════════════════════════════╗
║  INPUTS                                                      ║
║                                                              ║
║  · GFF3 gene annotation  (sorted, parent-before-child)       ║
║  · Query proteome FASTA  (translated gene models)            ║
║  · Reference proteome FASTA  (well-annotated close relative) ║
║  · Genome FASTA  (needed for translation validation, step 8) ║
║  · IsoSeq long-read transcript GFF  [ optional ]             ║
╚══════════════════════════════════════════════════════════════╝
                              ║
╔══════════════════════════════════════════════════════════════╗
║  STEP 1 — PREPARE                                            ║
║                                                              ║
║  · Remove query proteins containing internal stop codons or  ║
║    dots — reads full multi-line records before filtering     ║
║  · Extract gene features from GFF                            ║
║  · Extract IsoSeq mRNA features  (if IsoSeq GFF provided)    ║
╚══════════════════════════════════════════════════════════════╝
                              ║
╔══════════════════════════════════════════════════════════════╗
║  STEP 2 — PROTEIN HOMOLOGY SEARCH  (diamond blastp)          ║
║                                                              ║
║  · Query proteins searched against the reference proteome    ║
║  · Identifies which annotated gene models are fragmentary    ║
║    homologs — too short to represent a full-length ortholog  ║
╚══════════════════════════════════════════════════════════════╝
                              │
┌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┐
╎  STEP 3 — ISOSEQ OVERLAP MAPPING  [ if isoseq_gff is set ]  ╎
╎                                                              ╎
╎  · bedtools intersect: finds which long-read transcript      ╎
╎    models co-localize with each annotated gene model         ╎
╎  · Produces the read-to-gene overlap table used in step 5    ╎
└╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┘
                              │
╔══════════════════════════════════════════════════════════════╗
║  STEP 4 — FIND SPLIT-GENE CANDIDATES  (find_split_genes.pl)  ║
║                                                              ║
║  Fragment filter                                             ║
║    · Flags query proteins covering ≤85% of their best        ║
║      reference hit — these look like partial gene models     ║
║                                                              ║
║  Tiling test                                                 ║
║    · Same chromosome, same strand only                       ║
║    · Nested genes excluded from rank (don't inflate dist)    ║
║    · Do two adjacent flagged fragments together span the     ║
║      same reference protein end-to-end?                      ║
║      (±15 aa tolerance in reference protein coordinates)     ║
║    · If yes: the pair is a split-gene candidate              ║
║    · Flag LARGE_SPAN  if genomic footprint > 500 kb          ║
║    · Flag LARGE_SPAN_EXTREME  if footprint > 2 Mb            ║
║                                                              ║
║  Chaining                                                    ║
║    · Link supported pairs into multi-gene loci  (A→B→C)      ║
║    · Trim / split chains at weak asymmetric junctions        ║
║    · Assign quality flags to every candidate (table below)   ║
╚══════════════════════════════════════════════════════════════╝
                              │
┌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┐
╎  STEP 5 — ISOSEQ VALIDATION  [ if isoseq_gff is set ]       ╎
╎           (uses output from steps 3 and 4)                   ╎
╎                                                              ╎
╎  · Count long-read transcripts spanning 2+ genes per locus  ╎
╎    — molecular evidence of co-transcription across           ╎
╎    adjacent gene fragments                                   ╎
╎                                                              ╎
╎  · FULL_SPAN:    ≥1 read spans all genes in the locus        ╎
╎  · PARTIAL_SPAN: reads present but none reach a terminal     ╎
╎    fragment                                                  ╎
╎  · NO_SPANNERS:  no spanning reads found                     ╎
└╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┘
                              ║
╔══════════════════════════════════════════════════════════════╗
║  STEP 6 — MERGE  (merge_split_genes.pl)                      ║
║                                                              ║
║  · Filter candidates by flags, coverage, tiling hits, and   ║
║    IsoSeq status  (see merge filters table below)            ║
║  · PARTIAL_SPAN: optionally trim unsupported terminal genes  ║
║    before merging  (fix_partial = yes)                       ║
║  · Remove source gene models from the annotation             ║
║  · Rebuild one gene per locus — join exon / CDS / UTR        ║
║    features in genomic order; recalculate CDS phase          ║
║  · Assign new IDs; tag merged genes with source=Mender       ║
╚══════════════════════════════════════════════════════════════╝
                              ║
╔══════════════════════════════════════════════════════════════╗
║  STEP 7 — GT GFF3 CHECK  (fast, non-fatal)                   ║
║                                                              ║
║  · GFF3 spec validation on Mender-created genes only         ║
║  · Isolates merge-introduced errors from pre-existing ones   ║
╚══════════════════════════════════════════════════════════════╝
                              │
┌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┐
╎  STEP 8 — TRANSLATION VALIDATION                            ╎
╎           [ run_translation_validation = yes ]               ╎
╎                                                              ╎
╎  · Translate merged CDS sequences with gffread               ╎
╎  · DIAMOND blastp merged proteins vs reference proteome      ╎
╎  · Multiple-sequence alignment of merged + source proteins;  ╎
╎    score each predicted gene-fusion junction  (mafft/kalign) ╎
╎                                                              ╎
╎  · PASS:    clean translation, good coverage, high MSA score ╎
╎  · REVIEW:  translates well but junction MSA score borderline╎
╎  · FAIL:    internal stop codons or poor junction support    ╎
╎                                                              ╎
╎  FAIL merges are replaced by their source genes in the       ╎
╎  final validated GFF                                         ╎
└╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┘
                              │
┌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┐
╎  STEP 9 — AGAT GENE-MODEL CHECK  [ run_agat = yes ]         ╎
╎           (runs on PASS GFF; or full merged GFF if step      ╎
╎            8 was skipped)                                    ╎
╎                                                              ╎
╎  · CDS/exon structure · orphan features                      ╎
╎  · Parent–child feature coherence                            ╎
└╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┘
                              ▼
╔══════════════════════════════════════════════════════════════╗
║  OUTPUTS                                                     ║
║                                                              ║
║  new_merges.gff              full merged annotation          ║
║  removed_genes.gff           source genes removed            ║
║  transl_pass.gff3            translation-validated PASS      ║
║  transl_fail.gff3            FAIL merges  (for inspection)   ║
║  transl_review.gff3          borderline merges               ║
║  new_merges_validated.gff    recommended for downstream use; ║
║                              FAIL merges restored to source  ║
╚══════════════════════════════════════════════════════════════╝
```

---

## Candidate quality flags

Flags are assigned in step 4 and carried into the merge table. Multiple
flags are comma-separated in the `flag` column.

### Protein evidence strength

| Flag | Meaning | Default action |
|---|---|---|
| `CLEAN` | Directly adjacent, ≥2 tiling hits, combined coverage ≥60% | merge |
| `STRONG` | All junctions ≥3 tiling hits — strongest protein evidence | merge |
| `SINGLE_HIT` | All junctions have exactly 1 tiling hit. For single-copy genes this is the expected maximum — inspect `hit_desc` and `pident` before skipping | review |
| `WEAK_END` | Low-evidence terminal junction survived asymmetric trimming | merge with caution |
| `WEAK_INTERNAL` | Low-evidence internal junction survived chain splitting | merge with caution |
| `LOW_COV` | Combined reference coverage <60% — more likely domain sharing than a split gene | **skip** |

### Structural caution

| Flag | Meaning | Default action |
|---|---|---|
| `SKIPPED_GENE` | A non-adjacent gene sits inside the merge locus — review before merging | **skip** |
| `TRANSITIVE_JOIN` | One or more junctions in the chain lack direct pairwise tiling; connection inferred through a shared neighbour | review |
| `MULTI_ISOFORM_JOIN` | ≥1 source gene has multiple transcripts; merged gene will contain cross-product isoforms | review |
| `LARGE_SPAN` | Merged locus exceeds `large_span_warn` (default 500 kb) | review |
| `LARGE_SPAN_EXTREME` | Merged locus exceeds `large_span_extreme` (default 2 Mb) — very few vertebrate genes span this range | **skip** recommended |

### IsoSeq support (added by step 5)

| Flag | Meaning |
|---|---|
| `FULL_SPAN` | ≥1 long-read transcript spans all genes in the locus — strong confirmation |
| `PARTIAL_SPAN` | Reads present but none reach a terminal gene — terminal fragment may not belong; use `fix_partial` |
| `NO_SPANNERS` | No spanning reads found — may reflect expression timing, not gene structure |

---

## Merge filters  (step 6 — `[merge_filters]` config section)

| Parameter | What it controls | Default |
|---|---|---|
| `skip_flags` | Exclude candidates whose flag column contains any of these | `SKIPPED_GENE,LOW_COV` |
| `flags` | Include only candidates matching this flag (`all` = no filter) | `all` |
| `min_tiling` | Minimum `max_tiling_hits` to process a candidate | `1` |
| `min_cov` | Minimum combined reference coverage % | `0` |
| `require_isoseq` | Restrict to candidates with this IsoSeq flag (`FULL_SPAN`, `PARTIAL_SPAN`, `none`) | _(all)_ |
| `isoseq_min_spanning` | Minimum number of spanning reads required | `0` |
| `fix_partial` | Auto-trim unsupported terminal genes in `PARTIAL_SPAN` candidates | `yes` |
| `asym_trim` | Trim/split chains at weak asymmetric junctions during step 4 | `yes` |
| `large_span_warn` | Genomic span (bp) above which `LARGE_SPAN` is flagged | `500000` |
| `large_span_hard` | Genomic span (bp) above which `LARGE_SPAN_EXTREME` is flagged | `2000000` |

---

## Pipeline run options

| Option | What it does |
|---|---|
| `--config myproject.cfg` | Path to the project config file (required) |
| `--steps 1,2,4,6` | Run only the specified numbered steps |
| `--dry_run` | Print all commands without executing |
| `run_translation_validation = no` | Skip step 8 |
| `run_agat = no` | Skip step 9 |

---

## Step-by-step descriptions

---

### Inputs

Mender requires three core sequence/annotation files: a GFF3 annotation of
the genome being curated, the protein sequences translated from that
annotation (the *query proteome*), and a protein FASTA from a closely
related, well-annotated species (the *reference proteome*). The GFF3 must
be sorted so that parent features (gene, mRNA) appear before their children
(exon, CDS); unsorted files will cause downstream steps to silently
mis-reconstruct merged gene models.

The genome FASTA itself is only required for step 8 (translation
validation) — if you are skipping that step it does not need to be
provided.

IsoSeq long-read transcript data is entirely optional but strongly
recommended when available. Even partial IsoSeq coverage substantially
reduces false-positive merges because co-transcription of adjacent
fragments is direct molecular evidence that they belong to a single gene,
whereas protein homology alone cannot rule out recent tandem gene
duplication.

---

### Step 1 — Prepare

**Biological rationale.** DIAMOND blastp (step 2) will fail or produce
misleading alignments if the query FASTA contains sequences with internal
stop codons (`*`) or dots (`.`). These characters appear in translated
genome annotations when the underlying CDS has a frame-shift, a sequencing
gap, or an annotation error. A truly split gene model would be expected to
translate cleanly — it is just truncated relative to its ortholog. Keeping
corrupted sequences would waste alignment time and could generate spurious
low-coverage hits that mimic real split-gene signal.

**What happens.** Each protein FASTA record is read in full (including
multi-line sequences) and the terminal stop character, if present, is
stripped — a trailing `*` or `.` is standard annotation convention and
carries no information useful to DIAMOND. The remaining sequence is then
checked: if any internal `*` or `.` is found the entire record is
discarded. Gene features are extracted from the GFF (lines containing
`\tgene\t`) into a gene-only GFF for use in bedtools overlap operations
(step 3). If an IsoSeq GFF was provided, mRNA features are similarly
extracted.

**Inputs:** `proteome_fa`, `gff`, `isoseq_gff` (optional)  
**Outputs:** cleaned protein FASTA, gene GFF, IsoSeq mRNA GFF (if used)

---

### Step 2 — Protein homology search

**Biological rationale.** The core hypothesis driving Mender is that when a
single gene has been erroneously split into two adjacent models, each
fragment will align to the same reference ortholog — but each alignment
will cover only a portion of it. A 500 aa reference protein aligned by two
250 aa query proteins, one covering residues 1–250 and the other 251–500,
is the canonical split-gene signature. Conversely, a single well-annotated
gene should align across most of its ortholog's length.

DIAMOND blastp is used because it is orders of magnitude faster than
BLAST for large proteomes while producing essentially identical
sensitivity for this type of whole-proteome vs whole-proteome search.

**What happens.** A DIAMOND database is built from the reference proteome.
The cleaned query proteome (step 1 output) is searched against it using
tabular output format 6. Beyond the standard columns, the search explicitly
requests `scovhsp` (alignment coverage as a fraction of the reference
protein length), `slen` (reference protein length), `qcovhsp` (alignment
coverage as a fraction of the query protein length), and `qlen` (query
protein length). These four columns are required by step 4: `scovhsp` and
`slen` drive the primary fragment filter, while `qcovhsp` and `qlen` were
added specifically to catch split fragments of shorter proteins that would
otherwise be missed by the reference-coverage threshold alone. The e-value
cutoff and number of threads are configurable.

**Key parameters:** `evalue` (default 1e-5), `threads`  
**Inputs:** cleaned query proteome (step 1), reference proteome FASTA  
**Outputs:** `diamond.out` — tabular blastp results

---

### Step 3 — IsoSeq overlap mapping  *(conditional)*

**Biological rationale.** A long-read IsoSeq transcript that spans two
adjacent annotated gene models is direct molecular evidence that those two
gene models are co-transcribed — i.e., they are fragments of one gene. This
step pre-computes which IsoSeq transcripts overlap which gene models so
that step 5 can count spanning reads efficiently for every candidate locus.

**What happens.** `bedtools intersect` is run with the IsoSeq mRNA GFF as
the `-a` file and the gene feature GFF as `-b`. The output is parsed to
produce a two-column read-to-gene overlap table (IsoSeq transcript ID →
gene ID). This table is the input to step 5.

**Inputs:** IsoSeq mRNA GFF (step 1), gene GFF (step 1)  
**Outputs:** `overlaps.txt` — transcript-to-gene overlap table

---

### Step 4 — Find split-gene candidates

**Biological rationale.** This is the analytical core of Mender. Two
adjacent gene models are declared a split-gene candidate when their protein
products together tile a single reference ortholog end-to-end. "Tiling"
means the DIAMOND alignments of fragment A and fragment B cover
non-overlapping, complementary regions of the same reference sequence such
that their combined coverage approaches 100%.

The tiling criterion distinguishes a split gene from other common
scenarios: a tandem duplicate (both copies align independently across the
full ortholog), a multi-domain protein that shares only one domain with its
best reference hit (low combined coverage, caught by `LOW_COV` flag), or an
unrelated protein whose alignment is coincidentally short.

**Fragment filter.** A query protein is treated as a fragment candidate if
its best DIAMOND alignment covers ≤85% of the reference protein length
(`scovhsp`). A second filter also checks query-side coverage (`qcovhsp`):
a protein that is itself short but aligns well across its own length can
still be a fragment of a larger gene — the reference-coverage threshold
alone would not flag it. Both criteria are evaluated so that fragments of
shorter proteins are not missed.

```
Reference protein (500 aa):
1                                                500
|--------------------------------------------------|

Whole gene — aligns across most of its ortholog (scovhsp = 92%):
 |===================WholeGene====================|
 → not a fragment candidate (scovhsp > 85%)

Split fragment — aligns to only the N-terminal half (scovhsp = 40%):
 |=========GeneA=========|
                          ← 60% of reference uncovered
 → fragment candidate (scovhsp ≤ 85%)
```

**Tiling test.** For each pair of genomically adjacent fragment candidates
(within `max_dist = 4` gene positions of each other on the same chromosome
and on the same strand — opposite-strand pairs are discarded immediately),
Mender checks whether their DIAMOND alignment coordinates on
the same reference protein are complementary. Before adjacency is
evaluated, nested genes — genes whose coordinates fall entirely within
another gene — are excluded from the genomic rank. Without this step, a
small intronic gene embedded inside a larger locus would inflate the
apparent distance between the two flanking split fragments and prevent the
pair from being found. Nested genes are reported in the log but are
otherwise unaffected and remain in the output annotation.

**Genomic span flags.** Even when two fragments pass the tiling test, the
physical distance they would be joined across on the genome is a useful
sanity check against false positives from gene cluster groups — families
where multiple paralogous copies sit in tandem, each aligning to the same
reference protein independently. If the total genomic footprint of the
candidate locus exceeds 500 kb the candidate is flagged `LARGE_SPAN`; if
it exceeds 2 Mb it is flagged `LARGE_SPAN_EXTREME`. Large vertebrate genes
do exist (neurexins, ROBO receptors, dystrophin can exceed 1 Mb), so
`LARGE_SPAN` alone is not grounds for rejection — but it should prompt
cross-checking with IsoSeq spanning evidence. `LARGE_SPAN_EXTREME`
candidates should be added to `skip_flags` unless independent evidence
confirms the merge. Both thresholds are configurable via
`--large_span_warn` and `--large_span_extreme`.

A related flag, `SKIPPED_GENE`, is raised when a non-nested intervening
gene sits between the two fragments (i.e., `genomic_dist > 1`). This gene
may be an additional split fragment that failed filters, an unrelated gene,
or an artifact — the `skipped_genes` column in the merge table records its
ID for manual inspection. `SKIPPED_GENE` candidates are recommended for
`skip_flags` and manual review before merging.
the same reference protein are complementary. The key measure is
`combined_cov_pct`: the span from the leftmost alignment start to the
rightmost alignment end on the reference protein, divided by the reference
protein length. This is a span measure, not the sum of individual
coverages — a 5 aa gap between the two alignments is included in the span
and slightly inflates the percentage, but for typical gaps this is
negligible. The fragments must tile the reference within a configurable
`wiggle` tolerance (default ±15 aa, measured in reference protein
coordinates; set via `wiggle` in config). A positive tiling gap means the
alignments don't quite meet; a negative gap means they overlap slightly;
both are acceptable within tolerance. The maximum number of gene positions
separating the two candidates is controlled by `max_dist` (default 4,
configurable). Values above 4 substantially increase the number of
`SKIPPED_GENE` candidates and run time. If the test passes, the pair is a
split-gene candidate.

*Example:* Reference protein = 500 aa. Gene A aligns to positions 10–200;
gene B aligns to positions 205–480. Combined span = 480 − 10 + 1 = 471 aa.
`combined_cov_pct` = 471/500 = 94.2%. Tiling gap = 205 − 200 = 5 aa
(positive, within wiggle=15). This pair passes.

```
Reference protein (500 aa):
1                                                500
|--------------------------------------------------|

GeneA aligns positions 10–200  (N-terminal fragment):
 |=========GeneA=========|

GeneB aligns positions 205–480  (C-terminal fragment):
                           |========GeneB========|
                          ^
                    tiling gap = 5 aa  (positive, < wiggle=15) — PASSES

Combined span: [10 ──────────────────────── 480] = 471/500 = 94.2%

Too large a gap (gap = 40 aa, > wiggle=15) — FAILS:
 |====GeneA====|                    |====GeneB====|
               ←────── 40 aa ──────→

Overlap too large (overlap = 20 aa, > wiggle=15) — FAILS:
 |===========GeneA===========|
               |========GeneB========|
               ←── 20 aa overlap ───→
```

**num_tiling_hits.** For a pair to record a tiling hit against a given
reference protein, both genes must independently have a fragment-level
alignment against that same protein, and those alignments must tile. The
number of distinct reference proteins that independently support the pair
in this way is `num_tiling_hits`. High counts (5+) arise when the split
boundary falls at a structurally conserved position shared across a protein
family. Low counts (1–2) are common for single-copy genes and do not mean
false positive — they simply reflect limited representation in the reference
proteome. If the two fragments have diverged enough to preferentially hit
*different* reference proteins with no overlap in their hit lists, the pair
will not appear in the output at all. This is a known blind spot: because
IsoSeq validation in step 5 only operates on candidates already identified
by the protein homology search, it cannot rescue pairs that were never
found here. Such cases would need to be identified and added to the merge
table manually.

**Chaining.** When three or more adjacent fragments all tile the same
reference protein, they are linked into a chain: A→B→C. Each junction in
the chain is evaluated independently. A junction is considered weak if it
has exactly 1 tiling hit AND the maximum tiling hits seen at any other
junction in the chain is ≥ 6 (the asymmetry threshold). Weak terminal
junctions cause the chain to be trimmed (the terminal gene is dropped).
Weak internal junctions cause the chain to be split into two separate
candidates at that point. Both trimming and splitting are applied
iteratively before the chain is written to the merge table. Use
`--no_asym_trim` to disable this behaviour — useful when the reference
proteome is phylogenetically distant and overall tiling hit counts are
low.

```
Genome (5'→3'):   [═GeneA═] ─── [══GeneB══] ─── [GeneC]

Tiling hits per junction (reference protein positions shown):
  A→B: 8 hits  |=====A=====|=====B=====|           (strong)
  B→C: 1 hit               |=====B=====|==C==|     (weak; max in chain = 8)
                                               ↑
                          B→C is weak terminal — GeneC trimmed from chain

Result:   [═GeneA═] ─── [══GeneB══]   →  merge candidate
          [GeneC]                       →  excluded; remains in annotation
```

**Quality flags.** Every candidate is annotated with one or more quality
flags summarising the evidence (see flag table above). Flags are carried
forward through all subsequent steps and appear in the final merge table.

**Key parameters:** `large_span_warn`, `large_span_extreme`, `--no_asym_trim`  
**Inputs:** DIAMOND output (step 2), GFF, reference proteome FASTA, cleaned query proteome  
**Outputs:** `merge_candidates.txt`, `split_genes_summary.txt`, `split_genes_detail.txt`

---

### Step 5 — IsoSeq validation  *(conditional)*

**Biological rationale.** Protein homology identifies candidates based on
sequence similarity; IsoSeq long reads provide independent transcriptomic
evidence. A long-read transcript that spans two adjacent gene models — i.e.,
it overlaps both models without a gap between them — directly confirms that
the two are co-transcribed as a single pre-mRNA. This evidence is
particularly valuable for candidates flagged `SINGLE_HIT` where protein
evidence is thin.

**What happens.** For each candidate locus, Mender counts IsoSeq
transcripts that overlap all genes in the chain (using the overlap table
from step 3). A transcript qualifies as *spanning* if it overlaps every
gene in the locus, not just one or two of them. Three outcomes are
possible:

- **FULL_SPAN** — at least one transcript spans all genes in the locus.
  Strong confirmation; merge is well-supported.
- **PARTIAL_SPAN** — reads are present and span some genes in the chain
  but none reach a terminal fragment. The terminal fragment may be a
  neighbouring gene that was incorrectly chained; `fix_partial = yes`
  will trim it before merging.
- **NO_SPANNERS** — no spanning reads found. Importantly, absence of spanning
  reads does not mean the merge is wrong. IsoSeq libraries are typically
  generated from specific developmental stages or tissue types; a gene that
  is not expressed in those conditions will have no reads regardless of its
  structure. The key asymmetry: if reads are long enough to span two genes
  they should also reach a third if all three are truly one transcriptional
  unit — which is why `PARTIAL_SPAN` is informative even with incomplete
  libraries, but `none` is not.

```
Genome:     ──[═══GeneA═══]──────────────────[═══GeneB═══]──

FULL_SPAN   ─────────────────────────────────────────────────→
            one read spans both models — strong co-transcription support

PARTIAL_SPAN ──────────────────────────→
            reads cover GeneA but stop before reaching GeneB terminal
            (terminal may be a mis-chained neighbour; fix_partial = yes trims it)

NO_SPANNERS ──────────→         ──────────→
            reads present but none span the full locus
            (may reflect expression timing, not gene structure)
```

**Inputs:** `merge_candidates.txt` (step 4), `overlaps.txt` (step 3)  
**Outputs:** `isoseq_validated.txt` — merge table with IsoSeq status column added

---

### Step 6 — Merge

**Biological rationale.** Once candidates have been scored by protein
evidence (step 4) and optionally validated by IsoSeq (step 5), those that
pass the merge filters (see table above) — the surviving candidates — are
used to restructure the genome annotation. "Surviving" means they were not
excluded by `skip_flags`, met any minimum `min_tiling` or `min_cov`
thresholds, and satisfied any `require_isoseq` constraint set in the
config. The source gene models (the split fragments) are removed and
replaced by a single merged gene model whose exon/CDS/UTR structure is the
ordered union of all fragments in the locus.

**What happens.** Candidates are first filtered using the parameters in the
`[merge_filters]` config section (see table above). For `PARTIAL_SPAN`
candidates, if `fix_partial = yes`, the unsupported terminal fragment is
dropped before merging. The source genes are written to
`removed_genes.gff`, then the merged gene is constructed by collecting all
child features (exon, CDS, UTR) from every source model, sorting them by
genomic coordinate, and recalculating CDS phase from scratch. New unique
IDs are assigned using the `gene_template` and `trans_template` patterns;
every merged gene is tagged `source=Mender` in the GFF attributes field so
it can be distinguished from original annotation features.

*Example:* Fragments GeneA (exons on chr1:1000–1200, 1400–1600) and GeneB
(exons on chr1:2000–2200, 2400–2600) are merged. The resulting gene has
four exons spanning chr1:1000–2600, with CDS phase recalculated from
position 1000.

**Key parameters:** all `[merge_filters]` parameters, `gene_template`, `trans_template`  
**Inputs:** validated/merge table (step 4 or 5), GFF  
**Outputs:** `new_merges.gff` (full annotation with merges), `removed_genes.gff`

---

### Step 7 — GT GFF3 check

**Why this step exists.** Merging exon/CDS features from multiple gene
models is a complex structural operation. Step 7 is a transparency check:
it tells you whether the merge introduced any GFF3 format problems before
you commit the annotation to downstream analysis. It does not change any
files or flag any genes — all decisions rest with the user.

**What happens.** First, a gene count sanity check is run: the number of
genes in the input annotation, the merged output, and `removed_genes.gff`
are counted and verified against the expected value (input − removed +
merged). A mismatch here means genes were silently dropped or duplicated
during the merge and should be investigated before proceeding.

If `run_gt = yes`, `gt gff3validator` (GenomeTools) is then run on a
temporary GFF containing only the `source=Mender` features extracted from
the merged annotation. Isolating Mender-created features means pre-existing
GFF3 errors in the input annotation will not appear in the report.

**Errors are printed to the terminal and the pipeline continues regardless
of the outcome.** No genes are flagged and no files are modified. The user
must read the output and decide whether any errors warrant manual
intervention before proceeding to steps 8 and 9.

Common merge-introduced errors to watch for, and suggested actions:

- **Exon/CDS coordinates extend beyond parent mRNA.** Locate the gene ID
  in the merged GFF and inspect the coordinates manually or in a genome
  browser. To determine whether the error was introduced by the merge or
  was present in the original annotation, run AGAT on `removed_genes.gff`
  (the original source fragments, written in step 6). If the same error
  appears there, it is pre-existing. If it only appears post-merge, remove
  the offending rows from `isoseq_validated.txt` (or `merge_candidates.txt`
  if step 5 was skipped) and re-run step 6.

- **Incorrect or missing CDS phase.** Apply the same diagnostic: run AGAT
  on `removed_genes.gff` to distinguish pre-existing phase errors from
  merge-introduced ones. merge_split_genes.pl recalculates phase from
  scratch, so newly introduced errors typically point to an unusual input
  structure. Re-running step 6 with `--flags STRONG` restricts merging to
  high-confidence candidates and reduces exposure to structurally complex
  edge cases.

- **Duplicate feature IDs.** Check that `gene_template` and
  `trans_template` in the config produce IDs that do not overlap with the
  existing annotation namespace. Adjust the templates and re-run step 6.

**Inputs:** `new_merges.gff` (step 6)  
**Outputs:** gene count summary and GT validation messages printed to terminal; no files written

---

### Step 8 — Translation validation  *(conditional)*

**Biological rationale.** A merged gene model is structurally plausible if
its exon/CDS coordinates are internally consistent, but the ultimate test
is whether the predicted protein product makes biological sense. If a merge
joins two genomic regions that are not part of the same reading frame —
for example, because a short intervening gene was mistakenly included in
the chain, or because one fragment is on the wrong strand — the merged CDS
will contain internal stop codons or fail to align well to its reference
ortholog.

Translation validation provides a per-merge quality assessment with three
outcomes: PASS, REVIEW, or FAIL. Only PASS (and optionally REVIEW) merges
are carried forward to downstream use.

**What happens.** The merged CDS sequences are translated using `gffread`
and the resulting protein sequences are searched against the reference
proteome with DIAMOND. Each merge is also subjected to a multiple-sequence
alignment (MSA) using `mafft` or `kalign`. The alignment contains the
merged protein, the individual source fragment proteins, and up to three
reference proteome hits (preferring SwissProt manually-curated sequences
when available). The position of each gene-fragment junction within the
merged protein is calculated precisely: it is the cumulative CDS length of
all source genes up to that point divided by 3 (converting bp to amino
acids). That position is then mapped to the corresponding column in the
alignment.

The junction is scored by examining a ±5 amino acid window of alignment
columns centred on that column. Three sub-scores are computed and combined
as a weighted sum (total = 0–1; default PASS threshold = 0.5). The weights
must sum to exactly 1.0 (enforced at runtime) and are configurable via
`--w_conservation`, `--w_continuity`, and `--w_gap`. The defaults were
chosen because conservation and ref continuity are considered equally
important signals of a structurally sound junction, while gap pattern —
though informative — is treated as a supporting rather than primary
criterion:

```
Merged:   ... D  A  F  K | S  M  L  R  W  T  G ...
Ref-1:    ... D  A  F  K | S  M  L  R  W  T  G ...
Ref-2:    ... D  A  F  K | T  M  L  R  W  T  G ...
Ref-3:    ... D  A  F  K | S  M  L  R  W  T  G ...
Src-A:    ... D  A  F  K | -  -  -  -  -  -  - ...   (N-terminal fragment ends)
Src-B:    ... -  -  -  - | S  M  L  R  W  T  G ...   (C-terminal fragment begins)
                          ↑
                    junction column (merged aa 215)
               |←── ±5 aa window ──→|

Conservation  (×0.4): fraction of non-gap residues matching consensus per column,
                       averaged across the window — measures residue agreement at junction
Ref continuity(×0.4): fraction of reference sequences gap-free within ±3 aa of junction
                       — measures whether junction falls in a coherent structural region
Gap pattern   (×0.2): fraction of window columns where merged has gap but references do not
                       — penalises a merged sequence that cannot join cleanly
```

- **Conservation** (weight 0.4) — at each column in the window, the
  fraction of non-gap characters that match the most common residue,
  averaged across the window. A high score means all sequences — merged,
  source fragments, and references — agree on similar residues at the
  junction region. A bad merge shows up as disagreement: the residues
  flanking the junction in the merged protein are inconsistent with what
  the source fragments and references show at that position.

- **Ref continuity** (weight 0.4) — the fraction of reference proteins
  that are completely gap-free across the junction (±3 aa). A high score
  means references flow through the junction without gaps, expected when
  the junction falls in a conserved, structurally coherent region. A low
  score means references also gap out at that position, which may indicate
  a naturally variable region rather than a merge-induced break.

- **Gap pattern** (weight 0.2) — penalises columns where the merged
  protein has a gap but the reference sequences do not. This catches cases
  where the two fragments do not join cleanly at the protein level: the
  merged sequence must introduce gaps to align, revealing a discontinuity
  at the junction that is absent from the references.

For multi-fragment merges (A→B→C) each junction is scored independently;
the minimum score across all junctions determines the overall
PASS/REVIEW/FAIL classification. FAIL merges are replaced by their
original source genes in the validated output GFF.

**Key parameters:** `run_translation_validation`, `genome_fa`, `kalign_bin` / `mafft_bin`  
**Inputs:** `new_merges.gff` (step 6), genome FASTA, reference proteome  
**Outputs:** `transl_pass.gff3`, `transl_review.gff3`, `transl_fail.gff3`, `new_merges_validated.gff`

---

### Step 9 — AGAT gene-model check  *(conditional)*

**Biological rationale.** AGAT (Another Gff Analysis Toolkit) is a
specialist GFF3 parser that enforces strict gene-model coherence rules: CDS
features must be enclosed by their parent exon, parent–child ID
relationships must be consistent, and there must be no orphan features
lacking a parent. These rules go beyond what `gt gff3validator` checks and
are particularly relevant after merging because the merge operation
assembles child features from multiple source models that may have
originally had different ID namespaces or slightly inconsistent feature
boundaries.

**What happens.** AGAT is run on the final output GFF — the translation-
validated PASS GFF if step 8 ran, or the full merged annotation if step 8
was skipped. The corrected GFF that AGAT produces is saved alongside the
input with `_agat` appended to the filename (e.g. `transl_pass_agat.gff3`
or `new_merges_agat.gff`). Stderr and stdout are combined, filtered for
lines containing "error" or "warn" (case-insensitive), and the first 20
lines are printed to the terminal. This step is purely informational — it
surfaces gene-model coherence problems for manual review but does not
modify any other pipeline files.

**Key parameters:** `run_agat`  
**Inputs:** `transl_pass.gff3` (or `new_merges.gff` if step 8 skipped)  
**Outputs:** `transl_pass_agat.gff3` (or `new_merges_agat.gff`); first 20 error/warning lines printed to terminal

---

### Outputs

The primary deliverable is `new_merges_validated.gff` — the input
annotation with split gene fragments removed and replaced by merged,
translation-validated gene models. When step 8 is skipped, `new_merges.gff`
serves the same role.

`removed_genes.gff` records every source gene model that was deleted during
merging. This file is important for provenance: if a merge is later
determined to be incorrect it can be reversed by removing the merged gene
and restoring the source features from this file.

The three translation validation files (`transl_pass.gff3`,
`transl_review.gff3`, `transl_fail.gff3`) allow independent examination of
each quality tier. REVIEW merges in particular are worth manual inspection:
they translate well but their junction MSA score is borderline, often
because the fragment boundary happens to fall in a low-complexity or
gapped region of the alignment.

