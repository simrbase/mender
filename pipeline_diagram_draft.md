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
║  · Strip internal stop codons from query protein sequences   ║
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
║    · Do two adjacent flagged fragments together span the     ║
║      same reference protein end-to-end?  (±15 aa tolerance) ║
║    · If yes: the pair is a split-gene candidate              ║
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
╎  · none:         no spanning reads found                     ╎
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
| `none` | No spanning reads found — may reflect expression timing, not gene structure |

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
