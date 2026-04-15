# Mender — Flag Reference

All flags used across the pipeline, organized by the step that generates them.

---

## Step 4 flags — protein evidence  (`flag` column, merge table)

Assigned by `find_split_genes.pl`. Carried into the merge table and into the GFF
`flag` attribute on merged genes. These are the flags controlled by `skip_flags`
and `flags` in the `[merge_filters]` config section.

| Flag | Meaning | Filterable? |
|------|---------|-------------|
| `CLEAN` | Directly adjacent pair, ≥2 tiling hits, combined reference coverage ≥60% — no structural concerns | `skip_flags`, `flags` |
| `STRONG` | All junctions have ≥3 independent tiling hits — strongest protein evidence | `skip_flags`, `flags` |
| `SINGLE_HIT` | All junctions have exactly 1 tiling hit. For single-copy genes this is the expected maximum — inspect `hit_desc` and `pident` before skipping | `skip_flags`, `flags` |
| `WEAK_END` | A terminal junction has 1–2 tiling hits but survived asymmetric trimming | `skip_flags`, `flags` |
| `WEAK_INTERNAL` | An internal junction has 1–2 tiling hits but survived chain splitting | `skip_flags`, `flags` |
| `LOW_COV` | Combined reference coverage <60% — more likely domain sharing than a split gene | `skip_flags`, `flags` |
| `SKIPPED_GENE` | A non-adjacent gene sits inside the merge locus — review the `skipped_genes` column before merging | `skip_flags`, `flags` |
| `TRANSITIVE_JOIN` | One or more consecutive gene pairs in the chain have no direct pairwise tiling evidence; the chain connection is inferred transitively. `STRONG,TRANSITIVE_JOIN` warrants the same scrutiny as `WEAK_END` | `skip_flags`, `flags` |
| `MULTI_ISOFORM_JOIN` | At least one source gene has >1 transcript; the merged gene will contain cross-product isoform combinations not all of which are biologically real | `skip_flags`, `flags` |
| `LARGE_SPAN` | Merged locus genomic span exceeds `large_span_warn` (default 500 kb). Plausible for some gene families; review with IsoSeq for weak-evidence merges | `skip_flags`, `flags` |
| `LARGE_SPAN_EXTREME` | Merged locus span exceeds `large_span_hard` (default 2 Mb). Very few vertebrate genes span this range — recommended to add to `skip_flags` unless IsoSeq confirms | `skip_flags`, `flags` |

---

## Step 5 flags — IsoSeq spanning read support  (`isoseq_flag` column)

Assigned by `validate_merge_with_isoseq.pl`. Stored in `isoseq_flag` (column 20
of `isoseq_validated.txt`). Controlled by `require_isoseq` and
`isoseq_min_spanning` in the config.

| Flag | Meaning | Filterable? |
|------|---------|-------------|
| `FULL_SPAN` | ≥1 long-read transcript spans all genes in the locus — strong co-transcription evidence | `require_isoseq` |
| `PARTIAL_SPAN` | Spanning reads exist but none reach a terminal gene — terminal fragment may not belong; use `fix_partial = yes` to auto-trim | `require_isoseq` |
| `NO_SPANNERS` | No spanning reads found — reads may be present on individual fragments but none bridge more than one gene. May reflect expression timing or tissue, not gene structure. Not a negative result. | `require_isoseq` |

---

## Step 8 flags — translation validation  (`report.tsv`)

Assigned by `validate_merge_translation.pl`. Written into `report.tsv` and into
GFF attributes (`transl_result`, `msa_flag`, `min_junction_score`) on retained
merged genes in the final `validated.gff`.

### `translation_flag` (column 17)

| Value | Meaning | Notes |
|-------|---------|-------|
| `OK` | No internal stop codons in the merged CDS translation | — |
| `FRAMESHIFT_DETECTED` | Internal stop codon(s) found — see `internal_stop_pos` | Stops within 10 aa of the N or C terminus are not flagged (partial CDS tolerance) |

### `msa_flag` (column 16)

| Value | Meaning | Notes |
|-------|---------|-------|
| `GOOD_MSA` | Min junction MSA score ≥ `min_junction_score`; clean translation | — |
| `GOOD_MSA_LOW_REF` | Same as `GOOD_MSA` but fewer than `min_msa_refs` reference sequences were available — score is less reliable | Passes the same thresholds; use as a confidence indicator |
| `WEAK_JUNCTION` | Min junction score below `min_junction_score`, or aligner failed | — |
| `NO_HIT` | No diamond hit found for the merged protein | Treated as REVIEW, not FAIL |
| `SKIPPED` | MSA scoring was skipped (`no_msa = yes` in config) | — |

### `overall_result` (column 20)

| Value | Meaning | Controlled by |
|-------|---------|---------------|
| `PASS` | Clean translation, coverage thresholds met, junction MSA score above threshold | — |
| `FAIL` | At least one hard fail criterion met; source genes are restored in `validated.gff` | — |
| `REVIEW` | No hard fail, but PASS criteria not fully met — inspect before use | `final_gff_include_review` controls whether REVIEW merges are kept in `validated.gff` |

### `fail_reasons` (column 19)

| Value | Meaning |
|-------|---------|
| `INTERNAL_STOP` | Merged CDS has an internal stop codon |
| `LOW_MERGED_COV` | Merged protein covers less than `min_merged_cov` of its best reference hit |
| `NO_HIT` | No reference hit found for the merged protein |
