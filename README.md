# Mender

**Mender** is a pipeline for detecting and correcting erroneously split gene
annotations in genome GFF3 files. Split genes are adjacent gene models that
together cover a single reference protein — a common artifact of automated
gene finding tools.

Evidence comes from two independent sources:
1. **Diamond blastp** homology against a reference proteome (protein tiling)
2. **PacBio IsoSeq** long reads that span two or more split fragments (optional but recommended)

New merged genes are tagged with `source=Mender` and `merge_source=Mender` in
the output GFF so they are easy to identify downstream.

---

## How the Protein Homology Analysis Works

The core idea: if two adjacent genes in the genome each cover
a non-overlapping portion of the same reference protein, they are likely
two halves of a single gene that was erroneously split during automated
annotation.

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

### Step 1 — Fragment detection

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

### Step 2 — Tiling test

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

### Step 3 — Confidence from multiple reference proteins

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

### Step 4 — Chaining into multi-gene groups

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

### Step 5 — Quality flags

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

**Liberal run strategy:** set `skip_flags = LOW_COV` to pass nearly all
candidates through the merger, then use a downstream translation validation
step (translate merged CDS, BLAST against reference, check for internal stop
codons and low identity) as the primary quality gate. Flags become diagnostic
metadata explaining *why* a merged gene failed translation, rather than
pre-merge filters.

---

## Dependencies

| Tool | Purpose |
|------|---------|
| `perl` | Run all scripts |
| `diamond` | blastp of query proteome vs reference |
| `bedtools` | IsoSeq read–gene overlap (optional) |
| `gt` (GenomeTools) | GFF3 validation (optional) |
| `agat` | Biological GFF consistency check (optional) |

All tools are expected to be in `$PATH` unless explicit paths are set in the
config file.

Install only what you need. For example, if you are not using IsoSeq you can
skip `bedtools`; if you are skipping the validation step you can skip `gt` and
`agat`.

```bash
# Install the tools you need via conda (customize this list for your run)
conda install -c bioconda diamond bedtools genometools-genometools agat perl
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

### Multi-proteome strategy

Running Mender against two independent reference proteomes and taking the
intersection of supported merge candidates provides substantially stronger
evidence than a single reference run:

- Pairs supported by **both** references are the most reliable.
- Pairs supported by only one reference warrant closer inspection — the
  single-proteome hit may reflect an annotation artifact in that reference
  rather than genuine split-gene evidence.

To implement this, run Steps 1–4 twice (once per reference proteome), then
compare the two `split_genes_merge.txt` outputs by gene pair. A lightweight
post-processing script can produce a `consensus_merge_candidates.txt` with a
`proteome_support = N/M` column indicating how many references supported each
pair.

### Suggested references by taxon

| Taxon | Suggested references |
|---|---|
| Squamates (lizards, snakes) | *Anolis carolinensis* (RefSeq/Ensembl), *Pogona vitticeps* (Ensembl) |
| Mammals | Human GRCh38 (RefSeq), mouse GRCm39 (RefSeq) |
| Birds | Chicken GRCg7b (RefSeq), zebra finch |
| Fish | Zebrafish GRCz11 (RefSeq), medaka |

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

## Pipeline Steps

| Step | Name | Description |
|------|------|-------------|
| 1 | prepare | Clean protein FASTA (remove internal stops), extract GFF subsets |
| 2 | diamond | Run `diamond blastp` of query proteome vs reference proteome |
| 3 | bedtools | `bedtools intersect` to find IsoSeq reads overlapping gene models (skipped if no `isoseq_gff`) |
| 4 | find | `find_split_genes.pl` — identify split gene candidates from diamond output |
| 5 | isoseq | `validate_merge_with_isoseq.pl` — add spanning read evidence to merge table (skipped if no `isoseq_gff`) |
| 6 | merge | `merge_split_genes.pl` — apply merges to the GFF |
| 7 | check | `gt gff3validator` and `agat` validation of output GFF |

Steps 3 and 5 are automatically skipped when `isoseq_gff` is not set in the
config. In that case, step 6 uses protein evidence flags only.

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
- `split_genes_merge.txt` — one row per merge candidate (chained pairs)

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

perl validate_merge_with_isoseq.pl split_genes_merge.txt overlaps > validated_merge.txt
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
    validated_merge.txt input.gff output.gff
```

**Output files:**
- `output.gff` — merged annotation (all new features have `source=Mender`)
- `removed_genes.gff` — source genes removed (audit trail)
- `output.log` — merge ID → new gene ID mapping and run parameters

### `split_gene_report_checks.sh`
Bash utility with useful queries on all four output files. Run after step 4/5
to explore candidates before merging.

```bash
bash split_gene_report_checks.sh
```

---

## Output File Column Reference

### `split_genes_merge.txt` (17 columns)
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

### `validated_merge.txt` (20 columns)
Same as above plus:
```
18:spanning_isoseq_count  number of IsoSeq reads spanning 2+ genes
19:spanning_isoseq_detail per-read detail: read_id(Nof M):gene1+gene2+...
20:isoseq_flag            FULL_SPAN | PARTIAL_SPAN | none
```

---

## Flag Definitions

### Protein evidence flags (column 17)

**Evidence strength flags** (apply to all chain sizes):

| Flag | Meaning |
|------|---------|
| `CLEAN` | Directly adjacent, ≥2 tiling hits, coverage ≥60% — solid candidate |
| `STRONG` | All junctions ≥3 tiling hits — highest confidence protein evidence |
| `SINGLE_HIT` | All junctions have exactly 1 tiling hit. **Review, not a default skip**: for single-copy genes, 1 hit is the expected maximum. Check `hit_desc` and `pident`. |
| `WEAK_END` | Terminal junction has 1–2 hits but survived trimming |
| `WEAK_INTERNAL` | Internal junction has 1–2 hits but survived splitting |
| `LOW_COV` | Combined coverage <60% — likely domain sharing, not a split gene |

**Structural concern flags** (apply to all chain sizes):

| Flag | Meaning |
|------|---------|
| `SKIPPED_GENE` | A non-adjacent gene sits inside the merge locus — check `skipped_genes` column and GFF before merging |
| `TRANSITIVE_JOIN` | One or more consecutive gene pairs in the chain have no direct pairwise tiling evidence; connection is inferred transitively through other genes. `STRONG,TRANSITIVE_JOIN` warrants the same scrutiny as `WEAK_END`. `SINGLE_HIT,TRANSITIVE_JOIN` without IsoSeq is the highest-risk category. |
| `MULTI_ISOFORM_JOIN` | At least one source gene has >1 transcript; merged gene will contain cross-product isoform combinations. Verify isoform structure before using isoform-level annotations. |
| `LARGE_SPAN` | Merged locus genomic span exceeds `large_span_warn` (default 500 kb). Review with IsoSeq for weak-evidence merges. |
| `LARGE_SPAN_EXTREME` | Merged locus span exceeds `large_span_extreme` (default 2 Mb). Recommend adding to `skip_flags` unless IsoSeq confirms. |

### IsoSeq flags (column 20)
| Flag | Meaning |
|------|---------|
| `FULL_SPAN` | At least one read spans all genes — strong confirmation |
| `PARTIAL_SPAN` | Reads exist but none reach a terminal gene — structural signal that a terminal gene may not belong. Use `--fix_partial` to auto-trim. |
| `none` | No spanning reads found — may reflect expression stage, not gene structure. Not a negative result. |

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
    split_genes_merge.txt input.gff output.gff
```

---

## Examples

The `examples/` directory contains files from the *Chamaeleo calyptratus*
CCA3 genome annotation run (April 2026):

- `chacal.cfg` — project config used for the CCA3 run
- `pipeline_summary.txt` — full results summary and manual review notes
- `validated_merge.txt` — merge table with IsoSeq evidence (1067 candidates)
- `new_merges.log` — run log from `merge_split_genes.pl`
