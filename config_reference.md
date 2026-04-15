# Mender — Config Reference

All keys settable in `mender.cfg` (or any project config file), organized by
section. Copy `mender.cfg.example` to a project-specific name before editing.

---

## `[input]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `gff` | Path to genome annotation GFF3 | — | Must be sorted parent-before-child. Use `gt gff3 -sort -tidy -retainids` if needed |
| `genome_fa` | Path to genome sequence FASTA | — | Required for step 8 (translation validation) only |
| `proteome_fa` | Query protein FASTA (translated gene models) | — | Internal stop codons removed automatically in step 1 |
| `subject_fa` | Reference proteome FASTA (diamond database) | — | Prefer a well-annotated, phylogenetically close species. See README for selection guidance |
| `isoseq_gff` | IsoSeq long-read transcript GFF (mRNA features) | _(blank)_ | Leave blank to skip steps 3 and 5 |

---

## `[output]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `run_prefix` | Short identifier for this run — used to name all output files | — | Typically a date or version string; no slashes or spaces |
| `results_dir` | Directory for final per-run output files | `results` | Created automatically if absent |
| `workdir` | Directory for intermediate scratch files (steps 1–7) | `mender_workdir` | Safe to delete between runs |
| `db_dir` | Directory for shared diamond databases (subject and SwissProt) | `mender_db` | Reused across runs; not deleted between runs |

---

## `[diamond]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `threads` | Number of threads for diamond blastp | `4` | — |
| `evalue` | E-value cutoff for diamond blastp | `1e-20` | Only hits with evalue ≤ this value are kept |

---

## `[ids]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `gene_template` | Template for new merged gene IDs | `MERGE[DATE]g[GCOUNT:6]000` | Tokens: `[GCOUNT:N]` gene counter zero-padded to N digits; `[TCOUNT:N]` transcript counter; `[DATE]` today as YYYYMMDD |
| `trans_template` | Template for new merged transcript IDs | `MERGE[DATE]t[GCOUNT:6][TCOUNT:3]` | Counter starts above the highest existing ID of that pattern in the input GFF |

---

## `[merge_filters]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `skip_flags` | Exclude candidates whose `flag` column contains any of these | `SKIPPED_GENE,LOW_COV` | Comma-separated. See [flags_reference.md](flags_reference.md) for all filterable flags |
| `flags` | Include only candidates matching this flag | `all` | `all` = no include filter. Example: `flags = STRONG` to merge only high-confidence candidates |
| `min_tiling` | Minimum `max_tiling_hits` to process a candidate | `1` | 1 = no minimum. Higher values require more supporting reference proteins, but hit count reflects reference proteome depth (isoform count, copy number) as much as candidate quality — single-copy genes always score low regardless of how real the split is. Rely on flags instead |
| `min_cov` | Minimum combined reference coverage % to process a candidate | `0` | 0 = no minimum. 60 is redundant when `LOW_COV` is already in `skip_flags` |
| `require_isoseq` | Only process candidates with this exact `isoseq_flag` | _(blank = all)_ | Options: `FULL_SPAN`, `PARTIAL_SPAN`, `NO_SPANNERS` |
| `isoseq_min_spanning` | Minimum number of spanning IsoSeq reads required | `0` | 0 = no minimum. Only meaningful when `isoseq_gff` is set |
| `fix_partial` | Auto-trim unsupported terminal genes from `PARTIAL_SPAN` candidates before merging | `yes` | yes / no |
| `asym_trim` | Trim or split multi-gene chains at weak asymmetric junctions during step 4 | `yes` | yes / no. Set to `no` to see trimmed genes in output with `WEAK_END` flag instead of silently dropping them |
| `large_span_warn` | Genomic span (bp) above which `LARGE_SPAN` is flagged | `500000` | 0 = disable |
| `large_span_hard` | Genomic span (bp) above which `LARGE_SPAN_EXTREME` is flagged | `2000000` | 0 = disable |
| `max_dist` | Maximum gene rank positions between two candidate fragments to be considered a tiling pair | `4` | Nested genes excluded from the rank. Values above 4 substantially increase `SKIPPED_GENE` candidates |
| `wiggle` | Maximum allowed gap or overlap (aa) between two fragment alignments on the reference protein | `15` | Positive = gap; negative = overlap. Increase for distant reference proteomes |
| `asym_threshold` | Tiling-hit count on the strong side of a junction needed to trigger asymmetric trimming | `6` | Only applies when `asym_trim = yes` |
| `low_cov_thresh` | Combined reference coverage (%) below which `LOW_COV` is assigned | `60` | Controls the flag only; candidates are not removed unless `LOW_COV` is in `skip_flags` |

---

## `[validation]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `run_gt` | Run `gt gff3validator` on merged GFF (step 7, fast) | `yes` | yes / no. Checks strict GFF3 spec compliance; non-fatal |
| `run_translation_validation` | Run `validate_merge_translation.pl` (step 8) | `yes` | yes / no. Skip when `genome_fa` is unavailable or for a quick run |
| `min_junction_score` | Minimum per-junction MSA score for a merge to PASS | `0.5` | Range 0.0–1.0. Higher = stricter |
| `min_merged_cov` | Minimum fraction of merged protein covered by best reference hit | `0.60` | Relaxed vs the 0.85 fragment filter used in step 4 |
| `min_ref_cov` | Minimum fraction of best reference protein covered by merged protein | `0.50` | — |
| `no_msa` | Skip MSA entirely; report translation + diamond only | `no` | yes / no. Faster but loses junction scoring |
| `aligner` | MSA tool and strategy for junction scoring | `mafft_fast` | Options: `mafft_fast` (FFT-NS-2, fastest), `mafft_auto` (best strategy for input size), `kalign3` (5–10× faster than mafft_auto; recommended for large runs) |
| `w_conservation` | Weight for the conservation sub-score in junction scoring | `0.3` | Fraction of positions in ±5 aa window where merged and reference agree on majority residue. Weights must sum to 1.0 |
| `w_continuity` | Weight for the reference continuity sub-score | `0.3` | Fraction of reference sequences that are gap-free across the ±5 aa junction window |
| `w_gap` | Weight for the gap-pattern sub-score | `0.4` | Fraction of window columns where the merged protein does not introduce a gap vs references |
| `min_msa_refs` | Minimum number of reference sequences in the MSA for the junction score to be considered fully reliable | `2` | When fewer refs are available, merge is scored but flagged `GOOD_MSA_LOW_REF`. Set to 1 to suppress these warnings |
| `max_msa_refs` | Maximum number of reference sequences to include in each MSA | `3` | Diamond is run with `--max-target-seqs 5`; values above 5 have no effect unless diamond search depth is also increased |
| `keep_msa` | Write per-merge alignment files to `<results_prefix>_msa/` | `no` | yes / no. Useful for inspecting individual junction alignments; disk-intensive on large runs |
| `swissprot_fa` | Path to UniProt/SwissProt FASTA for MSA reference sequences | _(blank)_ | Falls back to `subject_fa` hits if blank or file not found |
| `swissprot_db` | Override path for pre-built SwissProt diamond database | _(blank)_ | Normally derived as `<db_dir>/swissprot.dmnd`; only set if using a DB outside `db_dir` |
| `final_gff_include_review` | Include REVIEW merges in the final `validated.gff` | `yes` | yes / no. `yes` keeps REVIEW merges (recommended when IsoSeq is available); `no` replaces them with source genes |
| `run_agat` | Run `agat_convert_sp_gxf2gxf.pl` gene-model check on PASS GFF (step 9, slow) | `yes` | yes / no. Checks CDS/exon structure, orphan features, parent–child coherence |

---

## `[paths]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `scripts_dir` | Path to Mender scripts directory | _(blank = current dir)_ | Leave blank when running from the scripts directory |
| `diamond_bin` | Path to `diamond` executable | _(blank = `$PATH`)_ | — |
| `bedtools_bin` | Path to `bedtools` executable | _(blank = `$PATH`)_ | — |
| `gt_bin` | Path to `gt` (GenomeTools) executable | _(blank = `$PATH`)_ | — |
| `agat_bin` | Path to `agat_convert_sp_gxf2gxf.pl` | _(blank = `$PATH`)_ | Include the script name in the path |
| `gffread_bin` | Path to `gffread` executable | _(blank = `$PATH`)_ | — |
| `mafft_bin` | Path to `mafft` executable | _(blank = `$PATH`)_ | — |
| `kalign_bin` | Path to `kalign` executable | _(blank = `$PATH`)_ | Only needed when `aligner = kalign3` |
| `perl_bin` | Path to `perl` executable | _(blank = `$PATH`)_ | — |
