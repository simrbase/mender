# Mender

**Mender** is a pipeline for detecting and correcting erroneously split gene
annotations in genome GFF3 files. Split genes are adjacent gene models that
together cover a single reference protein — a common artifact of automated
annotation tools such as Helixer.

Evidence comes from two independent sources:
1. **Diamond blastp** homology against a reference proteome (protein tiling)
2. **PacBio IsoSeq** long reads that span two or more split fragments (optional but recommended)

New merged genes are tagged with `source=Mender` and `merge_source=Mender` in
the output GFF so they are easy to identify downstream.

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
| Flag | Meaning |
|------|---------|
| `noflag` | Directly adjacent, ≥2 tiling hits, coverage ≥60% — solid candidate |
| `STRONG` | All junctions ≥3 tiling hits — highest confidence protein evidence |
| `SINGLE_HIT` | All junctions have exactly 1 tiling hit — thin evidence, rely on IsoSeq |
| `WEAK_END` | Terminal junction has 1–2 hits but survived trimming |
| `WEAK_INTERNAL` | Internal junction has 1–2 hits but survived splitting |
| `LOW_COV` | Combined coverage <60% — likely domain sharing, not a split gene |
| `SKIPPED_GENE` | A non-adjacent gene sits inside the merge locus — review manually |

### IsoSeq flags (column 20)
| Flag | Meaning |
|------|---------|
| `FULL_SPAN` | At least one read spans all genes — strong confirmation |
| `PARTIAL_SPAN` | Reads exist but none reach a terminal gene — structural signal that a terminal gene may not belong. Use `--fix_partial` to auto-trim. |
| `none` | No spanning reads found — may reflect expression stage, not gene structure |

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
