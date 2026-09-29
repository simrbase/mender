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
| `db_dir` | Directory for shared diamond databases (the subject proteome) | `mender_db` | Reused across runs; not deleted between runs |

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
| `skip_flags` | Exclude candidates whose `flag` column contains any of these | `SKIPPED_GENE,LOW_COV,PIECES_OVERLAP` | Comma-separated. See [flags_reference.md](flags_reference.md) for all filterable flags |
| `flags` | Include only candidates matching this flag | `all` | `all` = no include filter. Example: `flags = STRONG` to merge only high-confidence candidates |
| `min_tiling` | Minimum `max_tiling_hits` to process a candidate | `1` | 1 = no minimum. Higher values require more supporting reference proteins, but hit count reflects reference proteome depth (isoform count, copy number) as much as candidate quality — single-copy genes always score low regardless of how real the split is. Rely on flags instead |
| `min_cov` | Minimum `tiling_cov_pct` to process a candidate | `0` | 0 = no minimum. 60 is redundant when `LOW_COV` is already in `skip_flags` |
| `require_isoseq` | Only process candidates with this exact `isoseq_flag` | _(blank = all)_ | Options: `FULL_SPAN`, `PARTIAL_SPAN`, `NO_SPANNERS` |
| `isoseq_min_spanning` | Minimum number of spanning IsoSeq reads required | `0` | 0 = no minimum. Only meaningful when `isoseq_gff` is set |
| `fix_partial` | Auto-trim unsupported terminal genes from `PARTIAL_SPAN` candidates before merging | `yes` | yes / no |
| `asym_trim` | Trim or split multi-gene chains at weak asymmetric junctions during step 4 | `yes` | yes / no. Set to `no` to see trimmed genes in output with `WEAK_END` flag instead of silently dropping them |
| `large_span_warn` | Genomic span (bp) above which `LARGE_SPAN` is flagged | `500000` | 0 = disable |
| `large_span_hard` | Genomic span (bp) above which `LARGE_SPAN_EXTREME` is flagged | `2000000` | 0 = disable |
| `max_dist` | Maximum gene rank positions between two candidate fragments to be considered a tiling pair | `4` | Nested genes excluded from the rank. Values above 4 substantially increase `SKIPPED_GENE` candidates |
| `wiggle` | Maximum allowed gap or overlap (aa) between two fragment alignments on the reference protein | `15` | Positive = gap; negative = overlap. Increase for distant reference proteomes |
| `asym_threshold` | Tiling-hit count on the strong side of a junction needed to trigger asymmetric trimming | `6` | Only applies when `asym_trim = yes` |
| `low_cov_thresh` | `tiling_cov_pct` (coverage of the tiling reference by the best pair) below which `LOW_COV` is assigned | `60` | Controls the flag only; candidates are not removed unless `LOW_COV` is in `skip_flags` |
| `frag_min_qcov` | Minimum % of the query protein a fragment hit must cover (qcovhsp) | `30` | Removes pure domain-sharing hits from the fragment pool. The old default 60 also dropped real splits whose piece carries sequence the reference lacks (Chamaeleo: 114 kept merges more at 30); the rebuild and step 8 judge each merge. 0 = off |
| `ref_min_qcov` | `rebuild_candidates`: references every chain gene aligns over at least this % of itself rank first | `60` | Keeps a domain-only protein (often the closest by identity) from being tried first when `frag_min_qcov` is low |
| `end_tol` | aa a piece may extend past its alignment on the junction side | `30` | Controls the `JUNCTION_OVERHANG` flag only |
| `overlap_min_frac` | `PIECES_OVERLAP`: two pieces' alignments on one reference overlap by at least this fraction of the shorter alignment | `0.5` | Controls the flag only; candidates are not removed unless `PIECES_OVERLAP` is in `skip_flags` |
| `overlap_min_qcov` | `PIECES_OVERLAP`: each piece must align over at least this % of itself on that reference | `60` | Keeps shared domains out: two genes aligning 20% of themselves to one domain are not copies |

---

## `[validation]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `run_gt` | Run `gt gff3validator` on merged GFF (step 7, fast) | `yes` | yes / no. Checks strict GFF3 spec compliance; non-fatal |
| `run_translation_validation` | Run `validate_merge_translation.pl` (step 8) | `yes` | yes / no. Skip when `genome_fa` is unavailable or for a quick run |
| `min_merged_cov` | Minimum fraction of merged protein covered by best reference hit | `0.80` | The check that separates true from false merges (Congeria negative control: 0.80 kept 89/93 real merges, 12/100 deliberately false ones) |
| `min_ref_cov` | Minimum fraction of best reference protein covered by merged protein | `0.50` | Below it the merge is REVIEW |

The MSA junction scoring keys (`min_junction_score`, `no_msa`, `aligner`, `w_conservation`, `w_continuity`, `w_gap`, `min_msa_refs`, `max_msa_refs`, `keep_msa`, `swissprot_fa`, `swissprot_db`, and `mafft_bin`/`kalign_bin` in `[paths]`) were removed on 2026-09-28: in a negative control of 100 deliberately false merges every one scored GOOD_MSA (junction scores 0.59–1.00, the same range as real merges). They are ignored with a warning.
| `final_gff_include_review` | Include REVIEW merges in the final `validated.gff` | `yes` | yes / no. `yes` keeps REVIEW merges; `no` replaces them with source genes |
| `run_agat` | Run `agat_convert_sp_gxf2gxf.pl` gene-model check on PASS GFF (step 9, slow) | `yes` | yes / no. Checks CDS/exon structure, orphan features, parent–child coherence |

---

## `[miniprot]`

| Key | What it controls | Default | Notes |
|-----|-----------------|---------|-------|
| `run_miniprot` | Rebuild each merged gene with miniprot in step 6 | `yes` | The merge's tiling reference proteins are aligned across the merged locus, best-ranked first. SUPPORTED: the placeholder model is replaced (real splice sites, missing exon recovered, start/stop codons completed from the end fragments when needed, their UTRs put back where the rebuilt start/stop is theirs). NOT_SUPPORTED: step 8 fails the merge (`NO_PROTEIN_SPAN`). Congeria negative control: 100/100 false merges NOT_SUPPORTED |
| `flank` | bp added on each side of the merged gene for the miniprot search | `5000` | — |
| `min_prot_cov` | Minimum fraction of the reference protein the alignment must cover | `0.80` | Also required: no frameshift, no in-frame stop, CDS over every source gene |
| `max_refs` | How many ranked reference candidates (`rebuild_candidates`, best first) to try | `10` (all written) | The first SUPPORTED one is used; `rebuild_ref` (miniprot.tsv, report.tsv), `miniprot_ref_rank` / `miniprot_refs_tried` record which and how many |

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
| `miniprot_bin` | Path to `miniprot` executable | _(blank = `$PATH`)_ | Needed when `run_miniprot = yes` |
| `samtools_bin` | Path to `samtools` executable | _(blank = `$PATH`)_ | Needed when `run_miniprot = yes` (region extraction) |
| `perl_bin` | Path to `perl` executable | _(blank = `$PATH`)_ | — |
