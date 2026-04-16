# Case Study: Choosing a Reference Proteome for Mender

## Runs compared

| Config | Reference | Run prefix |
|--------|-----------|------------|
| `chacal_anole.cfg` | *Anolis carolinensis* Ensembl AnoCar2.0v2 all-isoform proteome | `20260414` |
| `chacal_sp.cfg` | UniProt SwissProt (2026 release, all taxa) | `CCA3_sp` |

All other parameters were identical: `skip_flags = SKIPPED_GENE,LOW_COV`,
`asym_trim = yes`, `asym_threshold = 6`, `min_tiling = 1`. Both runs used the
same SwissProt FASTA for step 8 MSA scoring (`swissprot_fa`), so translation
validation was directly comparable.

---

## Background: What the Reference Proteome Does

The reference proteome is the diamond blastp subject used in step 2 to identify
which gene models in the annotation share sequence similarity to the same protein.
Two models that each cover non-overlapping regions of the same reference protein
are the basic unit of split-gene detection. The reference therefore determines
what splits are *detectable at all*.

Two properties drive the tradeoff between reference choices:

**Coverage of the target genome's gene space** — a close relative's proteome
contains sequences for lineage-specific and diverged genes that pan-vertebrate
databases lack. A pan-vertebrate database covers conserved genes the close
relative's annotation may miss because it is itself fragmented.

**Number of independent sequences per gene** — each independent protein that
tiles a junction adds one to `max_tiling_hits`, pushing candidates from
SINGLE_HIT into CLEAN or STRONG. Anolis isoforms add depth for a single species;
SwissProt orthologs from multiple vertebrates add depth across phylogenetic
breadth. These produce different flag distributions with different biological
meanings.

---

# Part 1: Choosing Your Reference Proteome

*For users setting up a new Mender run and deciding what to use as the blast
reference.*

## Top-level results

| | Anolis | SwissProt |
|---|---|---|
| Step 4 candidates | 1,123 | 1,030 |
| Processed after skip_flags | 1,011 | 930 |
| PASS | 977 (96.6%) | 905 (97.3%) |
| REVIEW | 34 (3.4%) | 25 (2.7%) |
| FAIL | 0 | 0 |

Both runs produce 0 FAIL results. The difference in candidate count reflects
different detection coverage, not quality.

## Merge candidate overlap

| Category | Count |
|---|---|
| Shared by both runs | 773 |
| Anolis only | 238 |
| SwissProt only | 157 |

Neither reference finds all merges. 238 splits are visible only via Anolis; 157
are visible only via SwissProt. A single-reference run misses one of these sets
entirely.

## Why each reference detects a unique set

**Anolis carolinensis** detects splits in genes that have diverged from
pan-vertebrate orthologs or evolved within the squamate lineage. For a chameleon
gene with no close SwissProt entry, Anolis is the only reference that can tile
the junction. Lineage-specific gene expansions — duplicated or rapidly evolving
families unique to reptiles — are accessible only via a close relative.

**UniProt SwissProt** detects splits in ancient, conserved vertebrate genes.
The power comes from phylogenetic depth: a gene split in the chameleon whose
ortholog is intact in human, mouse, chicken, and zebrafish generates four
independent SwissProt tiling hits, each contributing to `max_tiling_hits`.
This is why SwissProt-only merges are CLEAN-dominated (108/157, 69%) while
Anolis-only merges are SINGLE_HIT-dominated (137/238, 58%): conserved genes
accumulate depth from vertebrate diversity; lineage-specific genes do not.

**The most important point:** if the Anolis annotation is itself fragmented at
a locus, Anolis cannot detect the split. The tiling algorithm requires a
reference protein that spans the junction between two chameleon fragments. If
the Anolis ortholog is split by the same annotation error, individual Anolis
fragments tile against each chameleon fragment separately but never bridge the
gap. SwissProt detects those loci because curated entries are drawn from
experimentally characterized proteins — intact, full-length sequences regardless
of what any one annotation looks like. The 157 SwissProt-only merges are not a
failure of Anolis sensitivity; they reflect a fundamental limit of using a draft
annotation as the sole reference.

## max_tiling_hits: the opposite shift from expectation

SwissProt was expected to collapse `max_tiling_hits` because it is non-redundant
within each species. The opposite occurred.

| max_tiling_hits | Anolis | SwissProt |
|---|---|---|
| 1 | 437 (39%) | 94 (9%) |
| 2 | 201 | 157 |
| 3 | 122 | 142 |
| 4 | 93 | 128 |
| 5–6 | 93 | 151 |
| 7+ | 177 | 358 |

SwissProt has far fewer single-hit candidates (9% vs 39%) and a higher
distribution overall. SwissProt's non-redundancy is per-species, not per-gene
family. A conserved split gene accumulates one tiling hit per species that has
an intact ortholog — a gene with representatives in human, mouse, rat, chicken,
and zebrafish reaches max_tiling=5 from SwissProt alone.

The thresholds CLEAN (≥2 hits), STRONG (≥3 hits), and the asymmetric trim
cutoff (≥6) are meaningful in both runs, but their biological interpretation
differs: Anolis hit counts reflect within-species isoform redundancy; SwissProt
hit counts reflect phylogenetic breadth across vertebrates.

## Why not combine the two databases?

The straightforward workaround — concatenating Anolis and SwissProt into a single
blast subject — does not work. `max_tiling_hits` counts become incoherent: a
conserved gene accumulates hits from both Anolis isoforms and SwissProt orthologs
across many species, potentially reaching max_tiling=30 or higher; a
lineage-specific gene still hits 1–2 sequences. The dynamic range expands by an
order of magnitude and the CLEAN, STRONG, and asymmetric trim thresholds — all
calibrated assuming consistent per-gene copy number within a single database —
lose their meaning. Candidates from conserved gene families are systematically
over-flagged relative to lineage-specific candidates.

## The sequential two-pass strategy

The correct approach is to run the pipeline twice in sequence: close relative
first, then SwissProt against the updated annotation. This requires no new code.

After pass 1, all Anolis-detectable merges are resolved and written to
`validated.gff`. The 157 SwissProt-only gene pairs remain as unmerged fragments
in the updated annotation — Anolis had no evidence for them and did not modify
them. Pass 2 detects them normally, treating the updated annotation as input.

The sequential approach also enables **emergent chain extension** beyond what
either single-pass run can find. In 9 cases (SP_SUPERSET), Anolis detected and
merged a 2-gene core (A+B). The individual fragments A and B were too short or
too diverged for a SwissProt reference to extend the chain to a third fragment C.
After merging, the AB protein is longer, covers more of the reference, and may
cross the tiling threshold. Pass 2 can then detect AB→C — a chain that neither
standalone run could produce.

### Setting up pass 2

**1. Regenerate the proteome from pass 1 output:**

```bash
gffread results/20260414/validated.gff \
  -g /path/to/genome.fasta \
  -y mender_workdir/proteome_pass2.fa
```

`gffread` translates CDS features from the updated GFF. The output includes
merged proteins from pass 1 alongside all unmodified gene proteins.

**2. Configure pass 2** — create `chacal_sp_pass2.cfg`:

```ini
gff         = results/20260414/validated.gff
proteome_fa = mender_workdir/proteome_pass2.fa
subject_fa  = /path/to/uniprot_sprot.fasta
run_prefix  = CCA3_sp_pass2
```

All other parameters carry over from `chacal_sp.cfg`. Use a new `run_prefix` so
pass 2 outputs land in `results/CCA3_sp_pass2/` without overwriting pass 1.

**3. Run pass 2:**

```bash
perl run_mender.pl --config chacal_sp_pass2.cfg
```

The pass 2 `validated.gff` is the final annotation — it contains pass 1 merges
plus all new pass 2 merges. No GFF merging step is needed.

## Reference choice by scenario

| Scenario | Recommended approach |
|---|---|
| Genome with a well-annotated close relative available | Close relative first, SwissProt second pass |
| Poorly annotated or deeply diverged genome | SwissProt first (or only); close-relative annotation may be too fragmented to bridge junctions |
| No IsoSeq data | Both passes; rely primarily on Tier 1 (shared) and SP-only CLEAN merges (see Part 2) |
| Gene family study with high paralogy risk | Close relative only, `min_tiling = 2`; SwissProt phylogenetic depth not needed and increases false-positive risk |

---

# Part 2: Reviewing Your Results — Confidence Tiers and REVIEW Calls

*For users who have run one or both passes and are reviewing their output.*

## Confidence tiers

Merge candidates can be stratified by the number of independent lines of evidence
supporting them. The table below uses the CCA3 dataset as a concrete example.

| Tier | Criteria | Count (CCA3) | Evidence basis |
|---|---|---|---|
| 1 | Shared by both runs | 773 | Two independent protein sources agree |
| 2 | Anole-only, FULL_SPAN IsoSeq | 145 | Lineage-specific or diverged gene, confirmed by spanning transcript |
| 3 | SP-only, CLEAN + FULL_SPAN | 73 | Multi-species protein signal + IsoSeq spanning |
| 4 | SP-only, PASS, no IsoSeq | 74 | Conserved gene, protein-only evidence |
| 5 | Anole-only, NO_SPANNERS | 91 | Single or double Anolis hit, no transcriptome support |

Tiers 1–3 (991 merges) have compound evidence from at least two independent
sources. These are the merges to accept without hesitation. Tiers 4–5 rest on a
single evidence type.

**Tier 5 warrants specific scrutiny.** These 91 Anolis-only merges have no
spanning transcripts and often a single Anolis protein tiling the junction.
At this evidence level, misassembly, transposable element insertion, or genuine
tandem duplication are plausible alternatives to a split gene. Setting
`min_tiling = 2` would remove most Tier 5 candidates; whether that is appropriate
depends on how well Anolis covers the chameleon gene space.

**IsoSeq support rates** across both unique sets are similar: 61% of Anole-only
merges (145/238) have FULL_SPAN support; 62% of SP-only merges (98/157) have
FULL_SPAN support. The NO_SPANNERS fraction is 38% in both sets — the same
biological problem appears regardless of reference.

## REVIEW calls in both runs

Both runs show the same bimodal pattern in min_junction_score for REVIEW calls:
a cluster at 0.41–0.50 (genuinely borderline junctions) and a cluster at 0.70–0.99
(cases where a secondary flag triggered REVIEW despite a clean MSA — e.g.,
GOOD_MSA_LOW_REF or WEAK_JUNCTION). Reference choice does not meaningfully
change the REVIEW rate or its cause. The triage approach from the SKIPPED_GENE
case study (open the `_aligned.fa` file; check whether source genes independently
span the full reference) applies equally to REVIEW calls from either run.

## How chain-boundary conflicts resolve in the sequential run

20 SP-only merge groups share individual gene IDs with Anole merge groups, but
have different chain boundaries. All resolve in the sequential run without
intervention:

| Conflict type | Count | How it resolves |
|---|---|---|
| SP_SUBSET: SP chain is a subset of a longer Anole chain | 10 | Anole merged those genes in pass 1; they no longer exist individually in pass 2 — the longer chain stands |
| SP_SUPERSET: SP extends an Anole pair by ≥1 gene | 9 | Pass 2 runs against the merged AB protein; may extend the chain to ABC if the merged protein crosses the tiling threshold |
| DIFFERENT_CHAIN: overlapping genes, different boundaries | 1 | One gene was consumed by a different Anole merge in pass 1; pass 2 cannot reconstruct the SP pair |

The DIFFERENT_CHAIN case (CCA3g018495/496/497) requires manual inspection to
determine which merge boundary is correct.

## Cases requiring manual review

**6-gene chain (CCA3g011919–011924):** REVIEW in the SwissProt run,
min_junction_score = 0.605. Six-fragment splits are uncommon and this junction
score is below threshold — inspect in the genome browser before including.

**DIFFERENT_CHAIN locus (CCA3g018495/496/497):** One gene was consumed by an
Anole merge with a different neighbor in pass 1. The SwissProt pair cannot be
reconstructed in pass 2. Inspect both candidate merges at this locus manually.

## SKIPPED_GENE candidates are not resolved by reference choice

All 42 Anolis SKIPPED_GENE candidates appear with the same flag in the SwissProt
run. Zero were cleared. Every intervening gene at these loci is a conserved
vertebrate gene with SwissProt entries — none are squamate-specific novelties
invisible to SwissProt. The SKIPPED_GENE flag fires on genomic position, not
blast evidence: changing the reference proteome will not remove it. The only
routes to recovering SKIPPED_GENE candidates are `spanning_rescue = yes` with
IsoSeq data, or removing `SKIPPED_GENE` from `skip_flags` (see the SKIPPED_GENE
case study).

---

## Practical recommendations

**For most vertebrate genome annotations:** run Anolis (or the closest available
relative) first, then SwissProt on the updated annotation. This is the recommended
production workflow.

**If no close-relative proteome is available:** a single SwissProt run is
sufficient for conserved genes. Lineage-specific splits will be missed, but
without a close relative there is no other option.

**After both passes:** stratify your merges by confidence tier before making
annotation decisions. Tiers 1–3 can be accepted. Tiers 4–5 warrant a closer look,
particularly any Anole-only SINGLE_HIT merges without transcriptome support.

**For CCA3, the recommended workflow is:**

1. Pass 1: `perl run_mender.pl --config chacal_anole.cfg`
2. Regenerate proteome with `gffread` from `results/20260414/validated.gff`
3. Pass 2: `perl run_mender.pl --config chacal_sp_pass2.cfg`
4. Final annotation: `results/CCA3_sp_pass2/validated.gff`

Manually review the 6-gene chain (CCA3g011919–011924) and the DIFFERENT_CHAIN
locus (CCA3g018495/496/497) before finalizing.

---

# Methods: Reproducing the Comparison Analysis

All commands assume you are in the mender working directory and that the two
runs are at `results/20260414/` (Anolis) and `results/CCA3_sp/` (SwissProt).
Adjust paths to match your run prefixes.

## 1. Merge candidate overlap

Extract the source-gene list (column 3) from each report and find the shared,
Anole-only, and SP-only sets by exact gene-list match.

```bash
awk -F'\t' 'NR>1 {print $3}' results/20260414/report.tsv | sort > /tmp/anole_genes.txt
awk -F'\t' 'NR>1 {print $3}' results/CCA3_sp/report.tsv  | sort > /tmp/sp_genes.txt

comm -12 /tmp/anole_genes.txt /tmp/sp_genes.txt | wc -l   # shared
comm -23 /tmp/anole_genes.txt /tmp/sp_genes.txt | wc -l   # Anole only
comm -13 /tmp/anole_genes.txt /tmp/sp_genes.txt | wc -l   # SP only

comm -23 /tmp/anole_genes.txt /tmp/sp_genes.txt > /tmp/anole_only.txt
comm -13 /tmp/anole_genes.txt /tmp/sp_genes.txt > /tmp/sp_only.txt
```

## 2. Flag distributions for unique sets

```bash
# Anole-only flag distribution
awk -F'\t' 'NR>1 {print $3,$NF}' OFS='\t' results/20260414/work/merge_candidates.txt \
  | grep -Ff /tmp/anole_only.txt | awk -F'\t' '{print $2}' | sort | uniq -c | sort -rn

# SP-only flag distribution
awk -F'\t' 'NR>1 {print $3,$NF}' OFS='\t' results/CCA3_sp/work/merge_candidates.txt \
  | grep -Ff /tmp/sp_only.txt | awk -F'\t' '{print $2}' | sort | uniq -c | sort -rn
```

## 3. IsoSeq support for Anole-only merges

Column 20 of `isoseq_validated.txt` holds the span flag (FULL_SPAN, PARTIAL_SPAN,
NO_SPANNERS). Column 20 of `report.tsv` holds the validation result.

```bash
awk -F'\t' 'NR>1 {print $3,$NF}' OFS='\t' \
    results/20260414/work/merge_candidates.txt | sort > /tmp/anole_merge_flags.txt
awk -F'\t' 'NR>1 {print $3,$20}' OFS='\t' \
    results/20260414/work/isoseq_validated.txt | sort > /tmp/anole_isoseq_flags.txt

grep -Ff /tmp/anole_only.txt /tmp/anole_merge_flags.txt | sort > /tmp/ao_mflags.txt
grep -Ff /tmp/anole_only.txt /tmp/anole_isoseq_flags.txt | sort > /tmp/ao_iflags.txt

join -t$'\t' /tmp/ao_mflags.txt /tmp/ao_iflags.txt \
  | awk -F'\t' '{print $2,$3}' OFS='\t' | sort | uniq -c | sort -rn
```

## 4. IsoSeq support for SP-only merges

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

## 5. max_tiling_hits distributions

Column 15 of `merge_candidates.txt` is `max_tiling_hits`.

```bash
awk -F'\t' 'NR>1 {print $15}' results/20260414/work/merge_candidates.txt \
  | sort -n | uniq -c

awk -F'\t' 'NR>1 {print $15}' results/CCA3_sp/work/merge_candidates.txt \
  | sort -n | uniq -c
```

## 6. Individual gene overlap check

This tests whether SP-only chains share any individual gene IDs with Anole chains
(not just identical chain lists).

```bash
# Expand SP-only chains to individual gene IDs
awk -F'\t' 'NR>1 {print $3}' results/CCA3_sp/report.tsv \
  | grep -Ff /tmp/sp_only.txt | tr ',' '\n' | sort -u \
  > /tmp/sp_only_individual_genes.txt

# Expand all Anole chains to individual gene IDs
awk -F'\t' 'NR>1 {print $3}' results/20260414/report.tsv \
  | tr ',' '\n' | sort -u \
  > /tmp/anole_all_individual_genes.txt

comm -12 /tmp/sp_only_individual_genes.txt /tmp/anole_all_individual_genes.txt | wc -l  # shared
comm -23 /tmp/sp_only_individual_genes.txt /tmp/anole_all_individual_genes.txt | wc -l  # SP-only exclusive
```

## 7. Chain-boundary conflict classification

For each SP-only group that shares a gene with an Anole group, compare chain
sizes to classify the relationship as SP_SUBSET, SP_SUPERSET, or DIFFERENT_CHAIN.

```bash
# Build list of SP-only groups that share genes with Anole chains
comm -12 /tmp/sp_only_individual_genes.txt /tmp/anole_all_individual_genes.txt \
  | while read gene; do
      grep -w "$gene" /tmp/sp_only.txt
    done | sort -u > /tmp/sp_overlap_groups.txt

while IFS= read -r sp_grp; do
  first_gene=$(echo "$sp_grp" | tr ',' '\n' | head -1)
  anole_grp=$(awk -F'\t' -v g="$first_gene" \
    'NR>1 {n=split($3,a,","); for(i=1;i<=n;i++) if(a[i]==g) print $3}' \
    results/20260414/report.tsv | head -1)
  sp_count=$(echo "$sp_grp"   | tr ',' '\n' | wc -l)
  anole_count=$(echo "$anole_grp" | tr ',' '\n' | wc -l)
  if   [ "$sp_count" -lt "$anole_count" ]; then echo "SP_SUBSET"
  elif [ "$sp_count" -gt "$anole_count" ]; then echo "SP_SUPERSET"
  else echo "DIFFERENT_CHAIN"
  fi
done < /tmp/sp_overlap_groups.txt | sort | uniq -c
```

## 8. REVIEW min_junction_score distributions

Column 14 of `report.tsv` is `min_junction_score`. Column 20 is the validation
result.

```bash
awk -F'\t' 'NR>1 && $20=="REVIEW" {print $14}' results/20260414/report.tsv | sort -n
awk -F'\t' 'NR>1 && $20=="REVIEW" {print $14}' results/CCA3_sp/report.tsv  | sort -n
```

## 9. SKIPPED_GENE flag cross-run check

Verify that SKIPPED_GENE candidates are identical across both runs.

```bash
awk -F'\t' 'NR>1 && $NF~/SKIPPED_GENE/ {print $3}' \
    results/20260414/work/merge_candidates.txt | sort > /tmp/anole_skipped.txt
awk -F'\t' 'NR>1 && $NF~/SKIPPED_GENE/ {print $3}' \
    results/CCA3_sp/work/merge_candidates.txt  | sort > /tmp/sp_skipped.txt

comm -12 /tmp/anole_skipped.txt /tmp/sp_skipped.txt | wc -l  # shared (expect 42)
comm -23 /tmp/anole_skipped.txt /tmp/sp_skipped.txt | wc -l  # Anole only (expect 0)
comm -13 /tmp/anole_skipped.txt /tmp/sp_skipped.txt | wc -l  # SP only (expect 0)
```
