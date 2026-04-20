# Mender

**Mender** is a pipeline for detecting and correcting erroneously split gene
annotations in genome GFF3 files. Split genes are adjacent gene models that
together cover a single reference protein — a common artifact of automated
gene finding tools.

Evidence comes from two independent sources:
1. **Diamond blastp** homology against a reference proteome (protein tiling)
2. **PacBio IsoSeq** long reads that span two or more split fragments (optional but recommended)

New merged genes are tagged with `source=Mender` in the output GFF and carry
a `merged_from` attribute listing the source gene IDs.

For an annotated pipeline diagram with per-step descriptions, see
[pipeline_overview.md](pipeline_overview.md).

For a worked example comparing reference proteome choices and a recommended
two-pass strategy, see
[case_study_reference_comparison.md](case_study_reference_comparison.md).

For a worked example comparing `max_dist` settings and a detailed analysis of
SKIPPED_GENE chain behavior (including strand effects and IsoSeq evidence), see
[case_study_maxdist_comparison.md](case_study_maxdist_comparison.md).

For a consolidated table of all flags across all pipeline steps, see
[flags_reference.md](flags_reference.md).

For a full reference of every config key with defaults and valid values, see
[config_reference.md](config_reference.md).

---

## Contents

- [Quick Start](#quick-start)
- [Dependencies](#dependencies)
- [Choosing a Reference Proteome](#choosing-a-reference-proteome)
- [Pipeline Steps](#pipeline-steps)
- [How the Protein Homology Analysis Works](#how-the-protein-homology-analysis-works)
- [Scripts](#scripts)
- [Output File Column Reference](#output-file-column-reference)
- [Flag Definitions](#flag-definitions) — see also [flags_reference.md](flags_reference.md) for a consolidated single-table view
- [New Gene ID Templates](#new-gene-id-templates)
- [Without IsoSeq](#without-isoseq)
- [Config Reference](config_reference.md) — full table of all config keys with defaults

---

## Quick Start

```bash
# Copy the example config and customize it for your project
cp mender.cfg.example myproject.cfg
# edit myproject.cfg — set input GFF, proteome FA, subject FA, and optionally isoseq_gff

# Run the full pipeline
perl run_mender.pl --config myproject.cfg

# Dry run — preview what would be merged without writing files
perl run_mender.pl --config myproject.cfg --dry_run

# Run specific steps only
perl run_mender.pl --config myproject.cfg --steps 4,6
```

`mender.cfg.example` is a fully annotated template — all options are documented
there. Do not edit it directly; copy it to a project-specific name first.

---

## Dependencies

| Tool | Purpose | Required? |
|------|---------|-----------|
| `perl` | Run all scripts | yes |
| `diamond` | blastp of query proteome vs reference | yes |
| `bedtools` | IsoSeq read–gene overlap | only if `isoseq_gff` is set |
| `gffread` | Translate merged and source CDS (step 8); export protein/CDS/transcript FASTAs (step 10) | required for steps 8 and 10; skip both to avoid this dependency |
| `mafft` | Multiple sequence alignment for junction scoring (Step 8) | only if translation validation is on, `no_msa = no`, and `aligner = mafft_fast` or `mafft_auto` |
| `kalign3` | Alternative MSA aligner (faster than MAFFT; recommended) | only if translation validation is on, `no_msa = no`, and `aligner = kalign3` |
| `gt` (GenomeTools) | GFF3 spec validation of merged genes (Step 7) | optional |
| `agat` | Biological GFF consistency check of passing merges (Step 9) | optional |

All tools are expected to be in `$PATH` unless explicit paths are set in the
config file.

```bash
# Install the tools you need via conda (customize this list for your run)
conda install -c bioconda diamond bedtools gffread mafft genometools-genometools agat perl
```

---

## Choosing a Reference Proteome

The reference proteome (the `subject_fa` input) is the database against which
your query proteins are searched. Its quality and relevance directly determine
how many split genes Mender can detect and how many false positives it will
produce. The ideal reference is both closely related to your organism and
well-annotated. These goals are sometimes in tension.

### Key criteria

**Completeness.** Prefer proteomes with high BUSCO vertebrata scores (>90%).
A gene absent from the reference creates a blind spot — split fragments whose
true ortholog is missing cannot be detected.

**Annotation quality.** RefSeq and Ensembl annotations are more reliable than
automated genome-only predictions. Avoid using as a reference a proteome that
was produced by the same gene caller you are trying to fix — you may inherit
the same split-gene artifacts you are trying to correct.

**Phylogenetic distance.** Closer is better for blast sensitivity and tiling
hit counts, but too close risks inheriting the same annotation errors. The
working rule: use the closest *well-annotated* species available, not the
closest *sequenced* species.

**Isoform representation.** Reference proteomes that include multiple isoforms
per gene (Ensembl typically does) produce higher `num_tiling_hits` for
multi-copy and multi-isoform gene families, because more reference sequences
can independently confirm the same tiling pattern. Proteomes filtered to one
canonical sequence per gene give more comparable hit counts across gene
families but reduce total sensitivity. Choose based on whether you prioritize
comparability (filtered) or maximum sensitivity (all isoforms).

### Finding a reference for your organism

RefSeq and Ensembl are the primary sources for curated proteomes across most
taxa. Both provide BUSCO completeness scores for their annotated assemblies,
which is the most direct way to compare candidates: a proteome with a high
BUSCO vertebrata (or metazoa, or embryophyta, depending on your clade) score
is more complete and will produce more tiling hits than a lower-scoring one.

Search NCBI Datasets or Ensembl for annotated assemblies in the same order or
class as your target organism. Download the protein FASTA for the best-scoring
candidate. If no closely related species has a high-quality annotation, a more
distant but well-annotated relative is preferable to a close but poorly
annotated one — a fragmented or error-prone reference proteome will introduce
false split-gene candidates and suppress real ones.

For organisms with no well-annotated relative at all (early-diverging lineages,
poorly studied phyla), a broad-coverage reference such as a curated set of
orthologs from OrthoDB or a UniProt/SwissProt subset filtered to your clade
can serve as a fallback, with the understanding that `num_tiling_hits` will be
lower throughout.

---

## Pipeline Steps

These are the **numbered steps** controlled by `--steps` in `run_mender.pl`.
Use `--steps 4,6` to run only specific steps, or omit `--steps` to run all.

| Step | Name | Script / Tool | Description |
|------|------|---------------|-------------|
| 1 | prepare | run_mender.pl | Clean protein FASTA (remove internal stops), extract GFF subsets |
| 2 | diamond | `diamond blastp` | Query proteome vs reference proteome |
| 3 | bedtools | `bedtools intersect` | IsoSeq reads vs gene models — skipped if no `isoseq_gff` |
| 4 | find | `find_split_genes.pl` | Identify split gene candidates from diamond output |
| 5 | isoseq | `validate_merge_with_isoseq.pl` | Add IsoSeq spanning read evidence — skipped if no `isoseq_gff` |
| 6 | merge | `merge_split_genes.pl` | Apply merges to the GFF; write new merged genes |
| 7 | gt | `gt gff3validator` | Fast GFF3 spec check on new Mender genes only; non-fatal |
| 8 | transl | `validate_merge_translation.pl` | Translate merged proteins, score junctions, assign PASS/FAIL/REVIEW |
| 9 | agat | `agat_convert_sp_gxf2gxf.pl` | Gene-model check on PASS GFF only; slow, optional |
| 10 | fasta | `gffread` | Export protein, CDS, and transcript FASTAs from the final validated GFF |

Steps 3 and 5 are automatically skipped when `isoseq_gff` is not set.
Step 8 is skipped when `run_translation_validation = no` in config.
Step 9 is skipped when `run_agat = no` in config.

### Step 8 internal sub-steps (A–L)

Step 8 (`validate_merge_translation.pl`) runs its own internal sequence of
lettered steps. These are **not individually selectable** via `--steps` — they
always run in order as part of step 8. Use `--no_msa` or `run_translation_validation = no`
to control the scope of step 8 from the config.

| Sub-step | Description |
|----------|-------------|
| A | Parse merged GFF — collect new gene IDs, source gene lists, GFF blocks |
| B | Parse original GFF — get representative transcript and CDS length per source gene |
| C | Translate merged proteins (gffread); write clean FASTA for diamond |
| D | Translate source gene proteins (gffread) |
| E | Build or verify diamond database from ref_fa |
| F | Diamond blastp merged proteins vs reference proteome; compute coverage metrics |
| G | Load reference proteome sequences into memory |
| H | SwissProt diamond blast (if configured); load hit sequences for MSA |
| I | Load merge table for flag annotations |
| J | MSA junction scoring — align merged + source + SwissProt ref (+ best ref_fa hit) proteins; score each junction. Aligner set via `aligner` config (kalign3 recommended). |
| K | Assign PASS / FAIL / REVIEW per merge; write `transl_result`, `msa_flag`, `min_junction_score` into GFF attributes |
| L | Write outputs into `--out` directory: `report.tsv`, `pass.gff3`, `fail.gff3`, `review.gff3`, `pass_proteins.fa` |

---

## How the Protein Homology Analysis Works

Two adjacent gene models are a split-gene candidate when their protein
products together tile a single reference ortholog end-to-end — each
covering a non-overlapping portion of it. The sections below describe how
that determination is made. These all happen inside pipeline step 4
(`find_split_genes.pl`).

### What is a "pair"?

A **pair** is two genes that are adjacent (or near-adjacent) to each other in
the genome — meaning they sit next to each other in the linear order of gene
models along a chromosome, with no other non-nested gene between them
(`genomic_dist = 1`). Mender also considers pairs with up to **4 gene
positions** between them (`genomic_dist` up to 4), in case an additional gene
sits inside the split locus. Nested genes (gene models that fall entirely
within the coordinates of another gene) are excluded from the rank when
calculating distances, since they would artificially inflate the apparent
distance between flanking genes.

### Fragment detection

Each gene model in the query annotation has a translated protein sequence.
These proteins are searched against a well-annotated reference proteome using
`diamond blastp`. Mender applies a **fragment filter** to identify hits where
the annotated gene appears to be a partial homolog rather than a complete one.
Two conditions are tested — either is sufficient:

- **Classic filter**: the alignment covers ≤85% of the reference protein
  (`scovhsp ≤ 85`). The gene is too short to represent the full-length
  ortholog.
- **Soft extension**: the alignment covers nearly the entire query protein
  (`qcovhsp ≥ 90`) but less than 90% of the reference protein (`scovhsp < 90`).
  This catches split fragments of shorter proteins — for example, a 200aa
  annotated CDS that is genuinely a fragment of a 280aa gene — which fall in
  the 86–89% scovhsp range and are missed by the classic filter alone.

Proteins that pass neither condition are treated as complete homologs and
excluded from split gene analysis.

### Tiling test

For each pair of adjacent genes, Mender checks whether their diamond
alignments tile end-to-end on a shared reference protein. The key requirement
is that the two alignments cover *different, non-overlapping* parts of the
reference, and that they meet in the middle with a gap of no more than **±15
amino acids** (the *tiling gap*, or *wiggle*). This tolerance accounts for
natural fuzziness at the ends of HSPs and is configurable when running
`find_split_genes.pl` manually.

**What tiling looks like:**

```
Reference protein (500aa):
1                                                 500
|--------------------------------------------------|

Gene A protein aligns to positions 10–210  (N-terminal fragment):
 |--------Gene A--------|

Gene B protein aligns to positions 215–480  (C-terminal fragment):
                          |----------Gene B--------|

                         ^
                    tiling gap = 5aa
```

Gene A and Gene B each cover less than half of the reference on their own, but
together they span most of it. Neither is a complete homolog — they look like
two halves of the same gene.

**Tiling gap — how the two alignments meet in the middle:**

The tiling gap is the distance between the end of the left alignment and the
start of the right alignment on the reference protein. Mender allows a gap of
up to ±15aa (configurable) to account for alignment boundary fuzziness:

```
Perfect join (gap = 0):
  Ref:  |---------Gene A---------|Gene B---------|
                                 ^ meet exactly — PASSES

Small gap (gap = +5aa, within ±15aa wiggle):
  Ref:  |---------Gene A---------|·····|Gene B---|
                                  5aa gap — PASSES

Small overlap (gap = −3aa):
  Ref:  |---------Gene A-----------|
                               |---Gene B--------|
                               ^^^  3aa overlap — PASSES

Gap too large (gap = +40aa):
  Ref:  |----Gene A----|                    |--Gene B--|
                        ← 40aa gap → FAILS, not a tiling pair

Overlap too large (gap = −20aa):
  Ref:  |---------Gene A--------------------|
                    |---Gene B--------------|
                    ←  20aa overlap → FAILS, likely same domain
```

**Combined coverage — how much of the reference protein the pair spans:**

Combined coverage is the total span from the leftmost alignment start to the
rightmost alignment end, divided by the full reference protein length. It
deliberately includes the tiling gap in the span, so it may slightly overstate
true coverage — but this is intentional and consistent.

```
Reference protein (500aa):
1                                                 500
|--------------------------------------------------|

Gene A: aligns  10–210
Gene B: aligns 215–480

Combined span = 480 − 10 + 1 = 471aa
Combined coverage = 471 / 500 × 100 = 94.2%

```

A combined coverage below **60%** receives a `LOW_COV` flag. This is more
likely to indicate that the two genes share a conserved domain (e.g. a
protein family with a common motif) rather than being a genuine split gene.
The threshold is configurable in the config file (`merge_filters` section) or
via `--min_cov` when running `merge_split_genes.pl` manually.

### Confidence from multiple reference proteins

A key confidence measure is **`num_tiling_hits`**: how many *different*
reference proteins independently support the same pair as tiling fragments.

For a reference protein to count, **both** genes must have a fragment-level
blast hit against that *same* protein, and both hits must tile on it. The
candidate pool is the **intersection** of each gene's hit list — proteins that
neither gene matched are never tested. This is a strict requirement.

**What drives high `num_tiling_hits`?**
High counts arise when the reference proteome contains multiple proteins —
paralogs or multi-species orthologs — that are each similar enough to *both*
fragments at the same time. This happens when the split cuts at a structurally
conserved boundary shared across the whole protein family (e.g., between a
ligand-binding domain and a signaling domain that is the same in all family
members). Every family member then independently confirms the same tiling
pattern.

**What drives low `num_tiling_hits`?**
Low counts are expected for single-copy genes with few representatives in the
reference, or when the split falls in a region that has diverged across
paralogs so that different family members are more similar to one half than the
other. This is not a signal of a false positive — it reflects the structure of
the reference proteome, not the biology of the split.

**An important blind spot:** if the two halves have diverged enough to
preferentially hit *different* reference proteins with no overlap between their
hit lists, the pair will not appear in the output at all — the blast
intersection is empty and the tiling test never runs. IsoSeq spanning reads
are the primary way to catch these cases.

| `num_tiling_hits` | Interpretation |
|---|---|
| 1 | Thin protein evidence — verify with IsoSeq before acting |
| 2–4 | Moderate evidence |
| 5+ | Strong — split boundary is conserved across many reference proteins |
| 10+ | Very strong |

### Chaining into multi-gene groups

Adjacent pairs that share a gene are chained together. For example, if gene
A–gene B passes the tiling test and gene B–gene C also passes, they become a
single three-gene merge candidate A→B→C. Chains of 3 or more genes are then
refined before writing to the merge table:

- **Trimming**: if a terminal junction has only 1 tiling hit while the rest of
  the chain is strongly supported (≥6 hits elsewhere), the terminal gene is
  dropped. This prevents a weakly supported gene from being dragged into a
  merge by a strong neighbor. Trimming is applied iteratively until neither
  terminal qualifies for removal.
- **Splitting**: if an internal junction has only 1 tiling hit while both of
  its flanking junctions are strongly supported, the chain is split at that
  junction into two separate merge candidates.

The asymmetry threshold for both trimming and splitting (default: 6) is
configurable when running `find_split_genes.pl` manually.

If a non-adjacent gene sits inside the merge locus (i.e. between two genes in
the chain), the candidate receives a `SKIPPED_GENE` flag. That skipped gene
may be another split fragment that failed the filters, or an unrelated gene —
it is worth checking before merging.

### Quality flags

Every merge candidate receives one or more quality flags summarizing the
protein evidence. Multiple flags are comma-separated in the `flag` column.

**Evidence strength flags:**

| Flag | Meaning |
|------|---------|
| `CLEAN` | Directly adjacent, ≥2 tiling hits, coverage ≥60% — solid candidate |
| `STRONG` | All junctions ≥3 tiling hits — highest confidence protein evidence |
| `SINGLE_HIT` | All junctions have exactly 1 tiling hit. **Review flag, not a default skip.** For single-copy genes in the reference proteome, 1 tiling hit is the expected maximum. Inspect `hit_desc` and `pident` before deciding whether to skip. |
| `WEAK_END` | Low-evidence terminal junction survived trimming |
| `WEAK_INTERNAL` | Low-evidence internal junction survived splitting |
| `LOW_COV` | Combined coverage <60% — likely domain sharing, not a split gene |

**Structural concern flags:**

| Flag | Meaning |
|------|---------|
| `SKIPPED_GENE` | Non-adjacent gene inside merge locus — review manually before merging |
| `TRANSITIVE_JOIN` | One or more consecutive genes in the chain have no direct pairwise tiling evidence; the chain connection is inferred transitively. Compounds other flags: `STRONG,TRANSITIVE_JOIN` warrants the same caution as `WEAK_END`. |
| `MULTI_ISOFORM_JOIN` | At least one source gene has multiple transcripts; the merged gene will contain cross-product isoform combinations not all of which may be biologically real. |
| `LARGE_SPAN` | The merged locus exceeds the `large_span_warn` genomic span threshold (default 500 kb). Plausible for some gene families but review with IsoSeq for weak-evidence merges. |
| `LARGE_SPAN_EXTREME` | The merged locus exceeds `large_span_extreme` (default 2 Mb). Very few vertebrate genes span this range. Recommended to add to `skip_flags` unless IsoSeq confirms. |

The recommended starting point is to merge `STRONG` and `CLEAN` candidates
while skipping `SKIPPED_GENE` and `LOW_COV`. These defaults are set in
`mender.cfg.example` and can be adjusted in the `[merge_filters]` section of
the config, or via command line options to `merge_split_genes.pl`.

---

## Scripts

### `run_mender.pl`
Pipeline wrapper. Reads `mender.cfg` and runs steps in order.

```
perl run_mender.pl --config mender.cfg [--steps 1,2,3] [--dry_run]
```

### `find_split_genes.pl`
Identifies split gene candidates from diamond blastp output. For each pair of
adjacent genes that together tile a reference protein, emits a candidate row.
Adjacent pairs sharing a gene are chained into multi-gene groups.

```
perl find_split_genes.pl diamond.out annotation.gff subject.fa query.fa
```

**Output files:**
- `split_genes_summary.txt` — one row per gene pair, best supporting hit only
- `split_genes_detail.txt` — one row per tiling hit per pair (full evidence)
- `merge_candidates.txt` — one row per merge candidate (chained pairs)

### `validate_merge_with_isoseq.pl`
Adds IsoSeq spanning read evidence to the merge table. A spanning read is a
long-read transcript that overlaps 2+ genes in a merge group — independent
molecular evidence that the genes are transcribed as one unit.

```
# Prepare overlaps first:
grep -P "\tgene\t"  annotation.gff > models.gene.gff
grep -P "\tmRNA\t"  isoseq.gff     > isoseq.mrna.gff
bedtools intersect -a isoseq.mrna.gff -b models.gene.gff -wa -wb | \
    perl -p -e 's/.+\sID=([^;]+).+\sID=([^;]+).*$/$1\t$2/' > overlaps

perl validate_merge_with_isoseq.pl merge_candidates.txt overlaps > isoseq_validated.txt
```

Adds three columns to the merge table: `spanning_isoseq_count`,
`spanning_isoseq_detail`, `isoseq_flag`.

### `merge_split_genes.pl`
Merges split genes in the GFF based on the validated merge table. Source genes
are removed and replaced by a new merged gene whose transcripts are built by
joining exon/CDS/UTR features in genomic order. CDS phases are recalculated.

```
perl merge_split_genes.pl \
    --fix_partial \
    --skip_flags SKIPPED_GENE,LOW_COV \
    --gene_template  "PREFIX_g[GCOUNT:6]" \
    --trans_template "PREFIX_t[GCOUNT:6][TCOUNT:3]" \
    isoseq_validated.txt input.gff output.gff
```

**Output files (standalone use):**
- `output.gff` — complete merged annotation; all Mender genes tagged with `source=Mender`
- `removed_genes.gff` — source genes removed during merge (audit trail)
- `output.log` — merge ID → new gene ID mapping and run parameters

When run via `run_mender.pl`, the merge outputs land in `results/<prefix>/`:
- `merges.gff` — complete merged annotation
- `removed.gff` — source genes removed
- `validated.gff` — **recommended for downstream use**; FAIL merges (and optionally REVIEW) replaced with their original source genes; retained Mender genes carry `transl_result`, `msa_flag`, `min_junction_score` GFF attributes (added by step 8)
- `validated.proteins.fa` — protein sequences from the final validated GFF (step 10)
- `validated.cds.fa` — CDS sequences from the final validated GFF (step 10)
- `validated.mrna.fa` — transcript sequences from the final validated GFF (step 10)
- `run_report.txt` — full run report: all parameters, step status, result counts, output file inventory

### `split_gene_report_checks.sh`
Bash utility with prebuilt queries for the candidate and merge table outputs
(steps 4–5) and for the translation validation report (step 8). Run from the
directory containing the files you want to inspect.

```bash
bash split_gene_report_checks.sh
```

---

## Output File Column Reference

### `merge_candidates.txt` (17 columns)
```
1:merge_id              unique merge group identifier
2:num_genes             number of genes in this merge
3:genes_in_order        comma-sep gene IDs in genomic order
4:coords_in_order       comma-sep coordinates (chr.start.end)
5:ref                   chromosome
6:gene_descs            comma-sep gene descriptions
7:best_hit              reference protein best representing the chain
8:hit_desc              description of best_hit
9:best_combined_cov_pct % of reference protein covered by the chain
10:best_gap             aa gap between the two best tiling alignments
11:evalues              comma-sep e-values vs best_hit
12:junction_tiling_hits pipe-sep hit counts per junction (chains only)
13:junction_genomic_dist pipe-sep rank distances per junction (chains only)
14:skipped_genes        genes sitting inside the merge locus (worth checking)
15:max_tiling_hits      highest independent tiling hit count across any pair
16:num_pairs_in_chain   number of pairwise relationships in the chain
17:flag                 quality flag(s) — see Flag Definitions below
```

### `isoseq_validated.txt` (20 columns)
Columns 1–17 are identical to `merge_candidates.txt`. Three columns are added:
```
18:spanning_isoseq_count  number of IsoSeq reads spanning 2+ genes
19:spanning_isoseq_detail per-read detail: read_id(Nof M):gene1+gene2+...
20:isoseq_flag            FULL_SPAN | PARTIAL_SPAN | NO_SPANNERS
```

### `report.tsv` (20 columns, step 8 output)
```
1:merge_id              merge group identifier (from merge table)
2:new_gene_id           Mender gene ID in the merged GFF
3:source_genes          comma-sep source gene IDs
4:merged_protein_len    length of representative merged protein (aa)
5:has_internal_stop     1 if merged protein contains internal stop codon(s)
6:internal_stop_pos     comma-sep positions of internal stops, or "."
7:merged_cov_by_ref     fraction of merged protein covered by best diamond hit
8:ref_cov_by_merged     fraction of reference protein covered by merged protein
9:best_ref_hit          reference protein ID with highest bitscore
10:best_ref_pident      % identity to best_ref_hit
11:ref_len              length of best_ref_hit (aa)
12:n_ref_hits           number of distinct reference hits found
13:junction_scores      pipe-sep per-junction MSA scores (one per source-gene junction)
14:min_junction_score   lowest per-junction score across all junctions
15:mean_junction_score  mean per-junction score
16:msa_flag             GOOD_MSA | GOOD_MSA_LOW_REF | WEAK_JUNCTION | NO_HIT | SKIPPED
17:translation_flag     OK | FRAMESHIFT_DETECTED
18:source_flags         flag column from the merge table (step 4/5)
19:fail_reasons         pipe-sep reasons for FAIL, or "."
20:overall_result       PASS | FAIL | REVIEW
```

---

## Flag Definitions

### Protein evidence flags (step 4 output, column 17 of merge table)

| Flag | Meaning |
|------|---------|
| `CLEAN` | Directly adjacent, ≥2 tiling hits, coverage ≥60% — solid candidate |
| `STRONG` | All junctions ≥3 tiling hits — highest confidence protein evidence |
| `SINGLE_HIT` | All junctions have exactly 1 tiling hit. **Review, not a default skip**: for single-copy genes, 1 hit is the expected maximum. Check `hit_desc` and `pident`. |
| `WEAK_END` | Terminal junction has 1–2 hits but survived trimming |
| `WEAK_INTERNAL` | Internal junction has 1–2 hits but survived splitting |
| `LOW_COV` | Combined coverage <60% — likely domain sharing, not a split gene |
| `SKIPPED_GENE` | A non-adjacent gene sits inside the merge locus — check `skipped_genes` column and GFF before merging |
| `TRANSITIVE_JOIN` | One or more consecutive gene pairs in the chain have no direct pairwise tiling evidence; connection is inferred transitively through other genes. `STRONG,TRANSITIVE_JOIN` warrants the same scrutiny as `WEAK_END`. `SINGLE_HIT,TRANSITIVE_JOIN` without IsoSeq is the highest-risk category. |
| `MULTI_ISOFORM_JOIN` | At least one source gene has >1 transcript; merged gene will contain cross-product isoform combinations. Verify isoform structure before using isoform-level annotations. |
| `LARGE_SPAN` | Merged locus genomic span exceeds `large_span_warn` (default 500 kb). Review with IsoSeq for weak-evidence merges. |
| `LARGE_SPAN_EXTREME` | Merged locus span exceeds `large_span_extreme` (default 2 Mb). Recommend adding to `skip_flags` unless IsoSeq confirms. |

### IsoSeq flags (column 20 of `isoseq_validated.txt`)

| Flag | Meaning |
|------|---------|
| `FULL_SPAN` | At least one read spans all genes — strong confirmation |
| `PARTIAL_SPAN` | Reads exist but none reach a terminal gene — structural signal that a terminal gene may not belong. Use `--fix_partial` to auto-trim. |
| `NO_SPANNERS` | No spanning reads found — reads may be present on individual fragments but none bridge more than one gene. May reflect expression timing or tissue, not gene structure. Not a negative result. |

### Translation validation flags (`report.tsv`)

**`translation_flag` (column 17):**

| Value | Meaning |
|-------|---------|
| `OK` | No internal stop codons in the merged CDS translation |
| `FRAMESHIFT_DETECTED` | Internal stop codon(s) found — see `internal_stop_pos` for positions |

**`msa_flag` (column 16):**

| Value | Meaning |
|-------|---------|
| `GOOD_MSA` | Min junction MSA score ≥ threshold; no internal stop |
| `GOOD_MSA_LOW_REF` | Same as `GOOD_MSA` but fewer than `min_msa_refs` reference sequences were available — score is less reliable |
| `WEAK_JUNCTION` | Min junction score below threshold, or aligner failed |
| `NO_HIT` | No diamond hit found for the merged protein |
| `SKIPPED` | MSA was skipped (`--no_msa`) |

**`fail_reasons` (column 19):**

| Value | Meaning |
|-------|---------|
| `INTERNAL_STOP` | Merged CDS has an internal stop codon |
| `LOW_MERGED_COV` | Merged protein covers less than `min_merged_cov` of its best reference hit |
| `NO_HIT` | No reference hit found for the merged protein |

**`overall_result` (column 20):**

| Value | Meaning |
|-------|---------|
| `PASS` | Clean translation, coverage thresholds met, junction MSA score above threshold |
| `FAIL` | At least one fail reason; merged gene should not be used without review |
| `REVIEW` | No hard fail criteria but PASS criteria not fully met — inspect before use |

---

## New Gene ID Templates

New gene and transcript IDs are controlled by `--gene_template` and
`--trans_template` (or the `[ids]` section in the config file). Template
tokens:

| Token | Meaning |
|-------|---------|
| `[GCOUNT:N]` | Gene counter, zero-padded to N digits |
| `[TCOUNT:N]` | Transcript counter, zero-padded to N digits |
| `[DATE]` | Today's date in YYYYMMDD format |

**Examples:**
```ini
# Generic date-stamped (default — safe for any organism)
gene_template  = MERGE[DATE]g[GCOUNT:6]000
trans_template = MERGE[DATE]t[GCOUNT:6][TCOUNT:3]

# Organism-specific (chameleon CCA3 style)
gene_template  = CCA3g1[GCOUNT:5]000
trans_template = CCA3t1[GCOUNT:5][TCOUNT:3]

# Medicago style
gene_template  = MS.gene[GCOUNT:5]
trans_template = MS.gene[GCOUNT:5].t[TCOUNT]
```

The counter starts above the highest existing ID of that pattern found in the
input GFF, so new IDs never collide with existing ones.

---

## Without IsoSeq

IsoSeq is optional. If `isoseq_gff` is left blank in the config, steps 3 and
5 are skipped. Use stricter protein-evidence filters to compensate:

```bash
perl merge_split_genes.pl \
    --flags STRONG \
    --skip_flags SKIPPED_GENE,LOW_COV \
    merge_candidates.txt input.gff output.gff
```
