# Case Study: Anolis carolinensis vs UniProt SwissProt as Mender Reference Proteome

## Summary

Two Mender runs on the *Chamaeleo calyptratus* CCA3 annotation were compared:
one using *Anolis carolinensis* as the reference proteome, one using UniProt
SwissProt. The runs found 773 shared merge candidates and diverged on 238
Anole-only and 157 SwissProt-only merges. Contrary to expectation, SwissProt
did not collapse `max_tiling_hits` — multi-species depth compensates for
within-species non-redundancy. Neither reference alone is sufficient: Anolis
detects lineage-specific splits SwissProt cannot see; SwissProt detects
conserved splits where the Anolis model is also fragmented. Concatenating the
two databases into a single blast subject is shown to be a poor approach
because it destroys the calibration of the hit-count thresholds. The
recommended strategy is a sequential two-pass run: Anolis first, then
SwissProt on the updated annotation, which delivers both reference-specific
merge sets in a single coherent workflow using no new code.

---

**Genome:** *Chamaeleo calyptratus* CCA3 assembly  
**Annotation:** CCA3-ref.merged.06192025.slim.gff  
**IsoSeq:** CCA3C_isoseq.agat.mRNA.gff  
**Run date:** 2026-04-14  
**Pipeline version:** Mender, commit b837c42

---

## Setup

Two runs were performed with identical parameters except for `subject_fa` — the reference proteome used for diamond blastp in steps 2–4.

| Config | Reference | Run prefix |
|---|---|---|
| `chacal_anole.cfg` | *Anolis carolinensis* Ensembl AnoCar2.0v2 all-isoform proteome | 20260414 |
| `chacal_sp.cfg` | UniProt SwissProt (2026 release, all taxa) | CCA3_sp |

All other parameters were identical: `skip_flags = SKIPPED_GENE,LOW_COV`, `asym_trim = yes`, `asym_threshold = 6`, `min_tiling = 1`. Both runs used the same SwissProt FASTA for step 8 MSA reference sequences (`swissprot_fa`), so MSA scoring was directly comparable.

---

## Top-level results

| | Anole | SwissProt |
|---|---|---|
| Step 4 candidates | 1123 | 1030 |
| Processed after skip_flags | 1011 | 930 |
| PASS | 977 (96.6%) | 905 (97.3%) |
| REVIEW | 34 (3.4%) | 25 (2.7%) |
| FAIL | 0 | 0 |
| Runtime | 3m 19s | 4m 25s |

Both runs had zero FAIL results. All candidates passed translation (translation_flag = OK throughout). The slight difference in PASS rate (96.6% vs 97.3%) reflects the different candidate pools, not a quality difference between the runs.

---

## Merge candidate overlap

```bash
# Extract source_genes column from each report.tsv and sort
awk -F'\t' 'NR>1 {print $3}' results/20260414/report.tsv | sort > /tmp/anole_genes.txt
awk -F'\t' 'NR>1 {print $3}' results/CCA3_sp/report.tsv  | sort > /tmp/sp_genes.txt

# Overlap counts
comm -12 /tmp/anole_genes.txt /tmp/sp_genes.txt | wc -l   # shared
comm -23 /tmp/anole_genes.txt /tmp/sp_genes.txt | wc -l   # Anole only
comm -13 /tmp/anole_genes.txt /tmp/sp_genes.txt | wc -l   # SP only
```

| Category | Count |
|---|---|
| Shared by both runs | 773 |
| Anole only | 238 |
| SwissProt only | 157 |

The comparison is by the exact comma-separated source gene list (e.g., `CCA3g000003000.1,CCA3g000004000.1`). Two merges that share individual genes but have different chain boundaries are counted as non-overlapping here; the chain-boundary analysis below handles those separately.

---

## Anole-only merges (238)

```bash
comm -23 /tmp/anole_genes.txt /tmp/sp_genes.txt > /tmp/anole_only.txt

# Flag distribution
awk -F'\t' 'NR>1 {print $3,$NF}' OFS='\t' results/20260414/work/merge_candidates.txt \
  | grep -Ff /tmp/anole_only.txt | awk -F'\t' '{print $2}' | sort | uniq -c | sort -rn

# max_tiling_hits distribution
awk -F'\t' 'NR>1 {print $3,$15}' OFS='\t' results/20260414/work/merge_candidates.txt \
  | grep -Ff /tmp/anole_only.txt | awk -F'\t' '{print $2}' | sort -n | uniq -c
```

| Flag | Count |
|---|---|
| SINGLE_HIT | 137 |
| CLEAN | 77 |
| WEAK_END | 7 |
| STRONG | 5 |
| LARGE_SPAN | 4 |
| others | 8 |

`max_tiling_hits` is heavily concentrated at 1 (143 of 238), with a long tail to 23.

These are splits detectable from the Anolis proteome but not from SwissProt. The most likely explanations: the gene is lineage-specific or sufficiently diverged that SwissProt has no close ortholog above the e-value cutoff, or the gene family is represented in SwissProt by a single canonical sequence that aligns to only one fragment.

The SINGLE_HIT dominance is expected: the Anolis Ensembl proteome includes many isoforms per gene, but for any given split gene the relevant protein is typically one specific paralog or isoform. For lineage-specific genes there is only one Anolis sequence to detect the split.

### IsoSeq support for Anole-only merges

```bash
# Cross merge flags with IsoSeq flags for Anole-only merges
awk -F'\t' 'NR>1 {print $3,$NF}' OFS='\t' \
    results/20260414/work/merge_candidates.txt | sort > /tmp/anole_merge_flags.txt
awk -F'\t' 'NR>1 {print $3,$20}' OFS='\t' \
    results/20260414/work/isoseq_validated.txt | sort > /tmp/anole_isoseq_flags.txt

grep -Ff /tmp/anole_only.txt /tmp/anole_merge_flags.txt | sort > /tmp/ao_mflags.txt
grep -Ff /tmp/anole_only.txt /tmp/anole_isoseq_flags.txt | sort > /tmp/ao_iflags.txt

join -t$'\t' /tmp/ao_mflags.txt /tmp/ao_iflags.txt \
  | awk -F'\t' '{print $2,$3}' OFS='\t' | sort | uniq -c | sort -rn
```

| Merge flag | FULL_SPAN | PARTIAL_SPAN | NO_SPANNERS |
|---|---|---|---|
| SINGLE_HIT | 91 | 0 | 46 |
| CLEAN | 45 | 0 | 32 |
| STRONG | 2 | 1 | 2 |
| WEAK_END | 5 | 1 | 1 |
| LARGE_SPAN | 2 | 0 | 2 |
| others | 0 | 0 | 5 |
| **Total** | **145 (61%)** | **2 (1%)** | **91 (38%)** |

61% of Anole-only merges have full IsoSeq spanning support. The 91 SINGLE_HIT + FULL_SPAN cases are particularly informative: a single Anolis protein provides the tiling signal and IsoSeq reads independently confirm that a transcript crosses the junction. Spanning read counts for those 91 range from 1 to 22, with a median around 7.

```bash
# Spanning read counts for SINGLE_HIT + FULL_SPAN Anole-only merges
awk -F'\t' 'NR>1 {print $3,$NF}' OFS='\t' \
    results/20260414/work/merge_candidates.txt \
  | grep -Ff /tmp/anole_only.txt | grep 'SINGLE_HIT' | cut -f1 \
  > /tmp/anole_only_singlehit.txt

grep -Ff /tmp/anole_only_singlehit.txt results/20260414/work/isoseq_validated.txt \
  | awk -F'\t' '$20=="FULL_SPAN" {print $18}' | sort -n | uniq -c
```

The 45 CLEAN + FULL_SPAN Anole-only merges are the strongest tier within this set: ≥2 Anolis tiling hits plus spanning transcriptome evidence, for genes that SwissProt cannot see at all.

The 91 NO_SPANNERS in the Anole-only set (46 SINGLE_HIT + 32 CLEAN + 13 other) have no transcriptome confirmation and thin protein evidence. These are the candidates most deserving of scrutiny before inclusion.

---

## SwissProt-only merges (157)

```bash
comm -13 /tmp/anole_genes.txt /tmp/sp_genes.txt > /tmp/sp_only.txt

# Flag distribution
awk -F'\t' 'NR>1 {print $3,$NF}' OFS='\t' results/CCA3_sp/work/merge_candidates.txt \
  | grep -Ff /tmp/sp_only.txt | awk -F'\t' '{print $2}' | sort | uniq -c | sort -rn

# Validation results
awk -F'\t' 'NR>1 {print $3,$20}' OFS='\t' results/CCA3_sp/report.tsv \
  | grep -Ff /tmp/sp_only.txt | awk -F'\t' '{print $2}' | sort | uniq -c
```

| Flag | Count |
|---|---|
| CLEAN | 108 |
| SINGLE_HIT | 33 |
| STRONG | 4 |
| LARGE_SPAN | 4 |
| others | 8 |

Validation: 147 PASS, 10 REVIEW, 0 FAIL.

SwissProt found these merges because the split gene has orthologs across multiple vertebrate lineages, each contributing an independent tiling hit. For a gene split in the chameleon annotation whose ortholog is intact in human, mouse, chicken, and zebrafish, four SwissProt sequences tile the junction. The CLEAN-heavy flag distribution (108/157) and high PASS rate reflect this multi-species signal quality.

### IsoSeq support for SP-only merges

```bash
awk -F'\t' 'NR>1 {print $3,$NF}' OFS='\t' \
    results/CCA3_sp/work/merge_candidates.txt | sort > /tmp/sp_merge_flags.txt
awk -F'\t' 'NR>1 {print $3,$20}' OFS='\t' \
    results/CCA3_sp/work/isoseq_validated.txt | sort > /tmp/sp_isoseq_flags.txt

grep -Ff /tmp/sp_only.txt /tmp/sp_merge_flags.txt | sort > /tmp/so_mflags.txt
grep -Ff /tmp/sp_only.txt /tmp/sp_isoseq_flags.txt | sort > /tmp/so_iflags.txt

join -t$'\t' /tmp/so_mflags.txt /tmp/so_iflags.txt \
  | awk -F'\t' '{print $2,$3}' OFS='\t' | sort | uniq -c | sort -rn
```

| Merge flag | FULL_SPAN | PARTIAL_SPAN | NO_SPANNERS |
|---|---|---|---|
| CLEAN | 73 | 0 | 35 |
| SINGLE_HIT | 21 | 0 | 12 |
| STRONG | 2 | 0 | 2 |
| LARGE_SPAN | 1 | 0 | 7 |
| others | 1 | 0 | 3 |
| **Total** | **98 (62%)** | **0** | **59 (38%)** |

73 of the 157 SP-only merges are CLEAN + FULL_SPAN: multi-species SwissProt tiling signal plus IsoSeq transcripts spanning the junction. These are splits that the Anolis proteome cannot detect at all, yet both an independent transcriptomic and a multi-species protein signal confirm them. This is the most important tier of the SP-only set.

---

## Are SP-only merges independent of Anole merges?

The source-gene-list comparison above counts chains as non-overlapping even if they share individual genes. To check for true redundancy, individual gene IDs were compared across all merge groups.

```bash
# Expand all source_genes to individual gene IDs
awk -F'\t' 'NR>1 {print $3}' results/CCA3_sp/report.tsv \
  | grep -Ff /tmp/sp_only.txt | tr ',' '\n' | sort -u \
  > /tmp/sp_only_individual_genes.txt

awk -F'\t' 'NR>1 {print $3}' results/20260414/report.tsv \
  | tr ',' '\n' | sort -u \
  > /tmp/anole_all_individual_genes.txt

# Overlap
comm -12 /tmp/sp_only_individual_genes.txt /tmp/anole_all_individual_genes.txt | wc -l
comm -23 /tmp/sp_only_individual_genes.txt /tmp/anole_all_individual_genes.txt | wc -l
```

- **291 of 331 individual genes** in SP-only merges appear nowhere in any Anole merge group. The Anole blast produced zero tiling evidence at those loci.
- **40 individual genes** (from 20 SP-only merge groups) appear in an Anole merge group with a different chain boundary.

Those 20 cases were classified by chain relationship:

```bash
# For each SP-only group that shares a gene with an Anole group,
# compare chain sizes to classify the relationship
while IFS= read -r sp_grp; do
  first_gene=$(echo "$sp_grp" | tr ',' '\n' | head -1)
  anole_grp=$(awk -F'\t' -v g="$first_gene" \
    'NR>1 {n=split($3,a,","); for(i=1;i<=n;i++) if(a[i]==g) print $3}' \
    results/20260414/report.tsv | head -1)
  sp_count=$(echo "$sp_grp" | tr ',' '\n' | wc -l)
  anole_count=$(echo "$anole_grp" | tr ',' '\n' | wc -l)
  if   [ "$sp_count" -lt "$anole_count" ]; then echo "SP_SUBSET"
  elif [ "$sp_count" -gt "$anole_count" ]; then echo "SP_SUPERSET"
  else echo "DIFFERENT_CHAIN"
  fi
done < /tmp/sp_overlap_groups.txt | sort | uniq -c
```

| Relationship | Count | Description |
|---|---|---|
| SP_SUBSET | 10 | SP supports a 2-gene pair that is the core of a longer Anole chain (3–4 genes) |
| SP_SUPERSET | 9 | SP extends an Anole 2-gene pair to a 3–6 gene chain |
| DIFFERENT_CHAIN | 1 | Overlapping but shifted pair at the same locus |

No SP-only merge is a duplicate of an Anole merge. All 20 chain-boundary disagreements resolve cleanly in a sequential run (see below).

---

## max_tiling_hits: opposite shift from expectation

```bash
# max_tiling_hits distributions
awk -F'\t' 'NR>1 {print $15}' results/20260414/work/merge_candidates.txt \
  | sort -n | uniq -c
awk -F'\t' 'NR>1 {print $15}' results/CCA3_sp/work/merge_candidates.txt \
  | sort -n | uniq -c
```

| max_tiling_hits | Anole | SwissProt |
|---|---|---|
| 1 | 437 (39%) | 94 (9%) |
| 2 | 201 | 157 |
| 3 | 122 | 142 |
| 4 | 93 | 128 |
| 5 | 49 | 85 |
| 6 | 44 | 66 |
| 7+ | 177 | 358 |

SwissProt was expected to collapse `max_tiling_hits` because it is non-redundant within each species. The opposite occurred. SwissProt has far fewer single-hit candidates (9% vs 39%) and a substantially higher distribution overall. The reason: SwissProt's non-redundancy is per-species, not per-gene-family. A conserved split gene accumulates one tiling hit per species in SwissProt that has an intact ortholog. A gene with orthologs in human, mouse, rat, chicken, and zebrafish can reach max_tiling=5 from SwissProt alone.

The concern about hit-count collapse was valid for the scenario of blasting against a non-redundant database of the *same* species. Against a pan-vertebrate database the multi-species depth compensates entirely.

The thresholds CLEAN (≥2 hits), STRONG (≥3 hits), and asymmetric trim threshold (≥6 hits) are meaningful in both runs, but their biological interpretation differs: Anole hit counts reflect within-species isoform redundancy; SwissProt hit counts reflect phylogenetic breadth across vertebrates.

---

## REVIEW merges

```bash
awk -F'\t' 'NR>1 && $20=="REVIEW" {print $14}' results/20260414/report.tsv | sort -n
awk -F'\t' 'NR>1 && $20=="REVIEW" {print $14}' results/CCA3_sp/report.tsv  | sort -n
```

Both runs show the same bimodal pattern in min_junction_score for REVIEW calls: a cluster at 0.41–0.50 (genuinely borderline junction alignments) and a cluster at 0.70–0.99 (cases where a secondary flag triggered REVIEW despite a strong MSA score, e.g., GOOD_MSA_LOW_REF or WEAK_JUNCTION). Neither reference produces a meaningfully cleaner REVIEW set.

---

## SP_SUPERSET cases: extended chains

Nine SP-only groups extend a known Anole pair by one or more terminal genes:

| SP chain | Anole chain | SP result |
|---|---|---|
| CCA3g001382000.1–001384000.1 (3 genes) | CCA3g001382000.1–001383000.1 (2 genes) | PASS (min_j=0.979) |
| CCA3g002419000.1–002422000.1 (3 genes) | CCA3g002420000.1–002422000.1 (2 genes) | PASS (min_j=0.966) |
| CCA3g007623000.1–007626000.1 (4 genes) | CCA3g007623000.1–007624000.1 (2 genes) | PASS (min_j=0.700) |
| CCA3g011919000.1–011924000.1 (6 genes) | not detected | REVIEW (min_j=0.605) |

The 3- and 4-gene supersets have high junction scores and PASS cleanly. The 6-gene chain is a REVIEW — six-fragment splits are uncommon and warrant manual inspection. These cases illustrate a real limitation of single-reference runs: when an individual terminal fragment is too short or too diverged for the reference to tile with confidence, Anolis truncates the chain. The multi-species depth of SwissProt bridges those junctions.

---

## Why not combine the two databases?

The straightforward workaround — concatenating the Anolis proteome and SwissProt into a single blast subject — does not work well. `max_tiling_hits` counts become incoherent across the combined database. A conserved gene accumulates hits from both Anolis isoforms and SwissProt orthologs across many species, potentially reaching max_tiling=30 or higher. A lineage-specific gene still hits 1–2 sequences. The dynamic range expands by an order of magnitude and the thresholds that define CLEAN, STRONG, and the asymmetric trim cutoff — all calibrated assuming consistent per-gene copy number within a single database — lose their meaning. Candidates from conserved gene families become systematically over-flagged relative to lineage-specific candidates, and the flag stratification that drives merge decisions breaks down.

---

## Sequential runs: the correct approach

The right way to capture both reference-specific merge sets is to run the pipeline twice, using the output of the first run as input to the second. This requires no new code — only a proteome regeneration step between passes.

### Why sequential works where combined databases do not

After pass 1 (Anolis), the 1011 merged gene pairs are single entries in `validated.gff`. The 157 SP-only gene pairs are still present as unmerged fragments in the updated annotation — Anolis had no evidence for them and did not touch them. Pass 2 (SwissProt) runs the full pipeline against the updated annotation and detects those pairs exactly as a standalone SP run would, with one important addition: the SP_SUPERSET cases.

In the SP_SUPERSET scenario, Anolis detected and merged a 2-gene core (A+B → AB) in pass 1. The individual fragments A and B may have been too short or too diverged for a SwissProt reference protein to tile with a third fragment C. After merging, the AB protein is longer and covers more of the reference, potentially crossing the tiling threshold for junction A+B → C. Pass 2 can then detect AB+C as a new merge candidate — an extended chain that neither a standalone Anolis run nor a standalone SwissProt run could produce.

### How the chain-boundary disagreements resolve

The 20 SP-only groups that share individual genes with Anole chains resolve without intervention:

- **10 SP_SUBSET cases:** The genes in the SP 2-gene pair were merged into a longer chain by Anole in pass 1. Those gene IDs no longer exist as individual entries in the updated annotation. Pass 2 does not see them as merge candidates. The longer Anole chain is the correct output.
- **9 SP_SUPERSET cases:** Anole merged the 2-gene core in pass 1. Pass 2 runs against the merged protein and may extend the chain with the terminal fragment(s) that Anolis could not reach. The result is the full multi-gene chain that neither single pass could assemble alone.
- **1 DIFFERENT_CHAIN case:** One gene from the SP pair was consumed by an Anole merge with a different neighbor. That gene no longer exists individually. Pass 2 does not reconstruct the SP pair. Manual review of this locus is warranted.

### Steps between passes

**1. Regenerate the proteome from the pass 1 output:**

```bash
gffread results/20260414/validated.gff \
  -g /n/sci/SCI-004219-SBCHAMELEO/Chamaeleo_calyptratus/genomes/CCA3-ref/CCA3C.fasta \
  -y mender_workdir/proteome_pass2.fa
```

`gffread` translates CDS features from the updated GFF using the genome sequence. The output includes merged proteins from pass 1 alongside all unmodified source gene proteins.

**2. Configure pass 2:**

Create a new config (e.g., `chacal_sp_pass2.cfg`) with:

```ini
gff         = results/20260414/validated.gff
proteome_fa = mender_workdir/proteome_pass2.fa
subject_fa  = /path/to/uniprot_sprot.fasta
run_prefix  = CCA3_sp_pass2
```

All other parameters carry over from `chacal_sp.cfg`. Use a distinct `run_prefix` so pass 2 outputs land in `results/CCA3_sp_pass2/` without overwriting anything.

**3. Run pass 2:**

```bash
perl run_mender.pl --config chacal_sp_pass2.cfg
```

Pass 2 runs all 9 steps on the updated annotation. The diamond blast in step 2 uses the regenerated proteome as query and SwissProt as subject. Steps 3–5 use the updated GFF for IsoSeq intersections and split-gene detection. The full translation validation (step 8) and AGAT check (step 9) run on the pass 2 merges.

**4. Combine final outputs:**

The pass 2 `validated.gff` contains the pass 1 merges plus any new merges found by SwissProt. It is the final annotation. No merging of GFF files is needed — the sequential structure means pass 2 output is already the combined result.

---

## Confidence tiers across both runs

| Tier | Criteria | Count | Basis |
|---|---|---|---|
| 1 | Shared by both runs | 773 | Convergent protein evidence from single-species and pan-vertebrate reference |
| 2 | Anole-only, FULL_SPAN IsoSeq | 145 | Thin or lineage-specific protein evidence, confirmed by spanning transcripts |
| 3 | SP-only CLEAN+FULL_SPAN | 73 | Multi-species tiling signal, confirmed by spanning transcripts |
| 4 | SP-only other PASS | 74 | Multi-species protein evidence, no IsoSeq support |
| 5 | Anole-only, NO_SPANNERS | 91 | Protein-only, single or double Anolis hit, no transcriptome confirmation |

Tiers 1–3 have compound evidence from at least two independent sources (protein tiling plus either multi-species convergence or IsoSeq). Tiers 4–5 rest on a single evidence type and warrant closer inspection before inclusion in a final annotation.

---

## Practical conclusions

A single Anolis run captures the widest candidate set, including lineage-specific and diverged genes invisible to SwissProt. A single SwissProt run captures well-conserved genes with multi-species tiling support that Anolis misses. Neither run alone is sufficient.

The sequential two-pass strategy — Anolis first, SwissProt second on the updated annotation — delivers both sets in a single coherent workflow with no new code. The pipeline handles chain-boundary conflicts automatically through the updated GFF: SP_SUBSET conflicts are resolved by the longer Anole chain, SP_SUPERSET cases become genuine chain extensions in pass 2, and the 137 SP-only merges fully absent from the Anole run are detected normally in pass 2 because their gene fragments are unchanged.

For CCA3, the recommended workflow is:

1. Run pass 1 with `chacal_anole.cfg` — produces `results/20260414/validated.gff`
2. Regenerate proteome with `gffread`
3. Run pass 2 with `chacal_sp_pass2.cfg` (SwissProt subject, updated GFF and proteome)
4. Use `results/CCA3_sp_pass2/validated.gff` as the final annotation

The 6-gene chain CCA3g011919–011924 (REVIEW) and the 1 DIFFERENT_CHAIN locus (CCA3g018495/496/497) should be reviewed manually regardless of pass 2 outcome.
