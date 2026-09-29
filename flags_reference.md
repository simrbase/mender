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
| `SINGLE_HIT` | All junctions have exactly 1 tiling hit. For single-copy genes this is the expected maximum — inspect `tiling_ref_desc` and `pident` before skipping | `skip_flags`, `flags` |
| `WEAK_END` | A terminal junction has 1–2 tiling hits but survived asymmetric trimming | `skip_flags`, `flags` |
| `WEAK_INTERNAL` | An internal junction has 1–2 tiling hits but survived chain splitting | `skip_flags`, `flags` |
| `LOW_COV` | `tiling_cov_pct` <60% (coverage of the tiling reference by the best pair) — more likely domain sharing than a split gene. In `skip_flags` by default | `skip_flags`, `flags` |
| `PIECES_OVERLAP` | Two pieces of the merge align to the **same region** of some reference protein (overlap ≥ `overlap_min_frac` of the shorter alignment, default 0.5), each over ≥ `overlap_min_qcov` (default 60%) of itself. Pieces of one split gene sit side by side on every reference; pieces on top of each other are copies of one gene unit (tandem duplicates), joined only because another reference is a fused two-unit protein they tile on. Every reference both pieces hit is checked, including hits the fragment filter drops (a piece covering > 85% of a reference is a complete copy of it). The reference is named in the `overlap_ref` column. The qcov condition keeps shared domains out; a real gene with internal repeats, cut between repeats, can still be flagged. In `skip_flags` by default. Congeria: 48 of 664 candidates, 4 of them PASS merges before the flag (merge_468: two mannose receptor-like copies, each ~99% of its own single-unit Dreissena protein, joined by a fused 1225 aa one) | `skip_flags`, `flags` |
| `JUNCTION_OVERHANG` | At some junction a piece extends more than `end_tol` aa (default 30) past its alignment on that junction's tiling reference, on the junction side — sequence the reference lacks at the break (a mispredicted exon, or two different proteins) | `skip_flags`, `flags` |
| `SKIPPED_GENE` | A non-adjacent **same-strand** gene sits inside the merge locus — the skipped gene may be an additional split fragment; review the `skipped_genes` column before merging. In `skip_flags` by default. Rescuable with `spanning_rescue = yes` if `isoseq_flag` is `FULL_SPAN`. Use `spanning_rescue` with caution for paralogous gene families — a FULL_SPAN call in repetitive or multi-copy loci may reflect a chimeric alignment. | `skip_flags`, `flags` |
| `OPPOSITE_STRAND_SKIP` | A non-adjacent gene sits inside the merge locus, but **all** skipped genes are on the opposite strand — these are unrelated interleaved genes and do not affect the validity of the merge. Not in `skip_flags` by default. | `skip_flags`, `flags` |
| `TRANSITIVE_JOIN` | One or more consecutive gene pairs in the chain have no direct pairwise tiling evidence; the chain connection is inferred transitively. `STRONG,TRANSITIVE_JOIN` warrants the same scrutiny as `WEAK_END` | `skip_flags`, `flags` |
| `MULTI_ISOFORM_JOIN` | At least one source gene has >1 transcript: the placeholder carries the cross-product (N x M transcripts); a gene rebuilt by miniprot (step 6) has one transcript and records the source isoforms it did not use (`isoforms_not_used=`); without the rebuild, step 8 drops each cross-product transcript that fails on its own | `skip_flags`, `flags` |
| `LARGE_SPAN` | Merged locus genomic span exceeds `large_span_warn` (default 500 kb). Plausible for some gene families; review with IsoSeq for weak-evidence merges | `skip_flags`, `flags` |
| `LARGE_SPAN_EXTREME` | Merged locus span exceeds `large_span_hard` (default 2 Mb). Very few vertebrate genes span this range — recommended to add to `skip_flags` unless IsoSeq confirms | `skip_flags`, `flags` |

---

## Step 5 flags — IsoSeq spanning read support  (`isoseq_flag` column)

Assigned by `validate_merge_with_isoseq.pl`. Stored in `isoseq_flag` (column 22
of `isoseq_validated.txt`). Controlled by `require_isoseq` and
`isoseq_min_spanning` in the config.

| Flag | Meaning | Filterable? |
|------|---------|-------------|
| `FULL_SPAN` | ≥1 long-read transcript spans all genes in the locus — strong co-transcription evidence. **Caveat:** spanning status is based on alignment coordinates only; for gene families with high paralogy or repetitive protein domains (collagen, HSP70, zinc finger arrays, etc.) a read may appear to span the locus due to mis-alignment between paralogous copies or alignment slippage through repetitive sequence. Verify in a genome browser for such families, especially when used with `spanning_rescue`. If the original BAM is available, filter on MAPQ ≥ 20 before deriving the GFF. | `require_isoseq` |
| `PARTIAL_SPAN` | Spanning reads exist but none reach a terminal gene — terminal fragment may not belong; use `fix_partial = yes` to auto-trim | `require_isoseq` |
| `NO_SPANNERS` | No spanning reads found — reads may be present on individual fragments but none bridge more than one gene. May reflect expression timing or tissue, not gene structure. Not a negative result. | `require_isoseq` |

---

## Step 8 flags — translation validation  (`report.tsv`)

Assigned by `validate_merge_translation.pl`. Written into `report.tsv` and into the
`transl_result` GFF attribute on retained merged genes in the final `validated.gff`.

### `miniprot` (column 13) and `rebuild_ref` (column 14) — from the step 6 rebuild

| Value | Meaning |
|-------|---------|
| `SUPPORTED` | The tiling reference protein aligns across the whole merged locus as one gene (≥ `min_prot_cov` of the protein, every source fragment, no frameshift, no in-frame stop); the model was rebuilt from that alignment |
| `NOT_SUPPORTED:<reason>` | It does not: `covers_N_of_M_sources`, `frameshift`, `inframe_stop`, `protein_cov_X`, `no_alignment`, `no_reference_protein` |
| `not_run` | `run_miniprot = no` |

`rebuild_ref` is the reference protein the model was rebuilt from (the first of the merge's
`rebuild_candidates` that is SUPPORTED); for NOT_SUPPORTED it is the best-ranked candidate, whose
alignment is reported. `miniprot.tsv` has it with its description (`rebuild_ref_desc`). It is
often not the merge's `tiling_ref`: detection uses the reference the pieces tile on best, the
rebuild tries the closest reference first.

Rebuilt genes also carry `start_codon=present|from_fragment:<transcript>|incomplete`,
`stop_codon=present|from_fragment:<transcript>|absent`, and `utr5=` / `utr3=`
`from_fragment:<transcript>|none` (each from the end fragment's first transcript, longest CDS
first, that fits there). When a source gene had several transcripts: `placeholder_transcripts=`,
`source_isoforms=<gene>:<n>`, `isoforms_not_used=` (a rebuilt gene has one transcript).
A source fragment's UTR is put back only where the rebuilt CDS
starts (ends) at exactly the fragment's own start (stop) codon.

### `translation_flag` (column 15)

| Value | Meaning | Notes |
|-------|---------|-------|
| `OK` | No internal stop codons in the merged CDS translation | — |
| `FRAMESHIFT_DETECTED` | Internal stop codon(s) found — see `internal_stop_pos` | — |

### `overall_result` (column 18)

| Value | Meaning | Controlled by |
|-------|---------|---------------|
| `PASS` | No fail reason; merged protein covers ≥ `min_merged_cov` of its best reference and the reference ≥ `min_ref_cov` of it | — |
| `FAIL` | At least one fail reason; source genes are restored in `validated.gff` | — |
| `REVIEW` | No fail reason, but the merged protein covers < `min_ref_cov` of its reference | `final_gff_include_review` controls whether REVIEW merges are kept in `validated.gff` |

### `fail_reasons` (column 17)

| Value | Meaning |
|-------|---------|
| `INTERNAL_STOP` | Merged CDS has an internal stop codon |
| `FRAME_BREAK` | A fragment junction leaves the downstream fragment out of frame (`frame_break` attribute from step 6) |
| `NO_PROTEIN_SPAN` | miniprot NOT_SUPPORTED: the reference protein does not span the merge as one gene |
| `LOW_MERGED_COV` | Merged protein covers less than `min_merged_cov` of its best reference hit |
| `NO_HIT` | No reference hit found for the merged protein |
