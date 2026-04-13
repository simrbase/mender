# SPEC: validate_merge_translation.pl

## Purpose

Post-merge translation quality validator. After `merge_split_genes.pl` produces a
merged GFF, this script:

1. Extracts and translates the merged CDS using `gffread -y`
2. Detects internal stop codons → flags `FRAMESHIFT_DETECTED`
3. Runs `diamond blastp` of merged proteins against the reference proteome to
   compute bidirectional length coverage
4. Aligns the merged protein + its best reference hit(s) + the original source
   gene proteins using MAFFT; scores the alignment at each gene-gene junction
5. Writes a per-merge TSV report and a filtered "pass/fail" summary
6. Optionally writes passing merged proteins to a FASTA for downstream use

This script is designed for the **liberal run strategy**: merge aggressively (few
skip_flags), then use this script as the primary quality gate. Users keep only
merges that pass translation validation.

---

## Position in Pipeline

```
Step 1  prepare          (run_mender.pl)
Step 2  diamond          (run_mender.pl)
Step 3  bedtools         (run_mender.pl)
Step 4  find             find_split_genes.pl
Step 5  validate_isoseq  validate_merge_with_isoseq.pl   [optional]
Step 6  merge            merge_split_genes.pl
Step 7  check            gt + agat validation             [optional]
        ↓
Step 8  validate_transl  validate_merge_translation.pl   ← THIS SCRIPT
```

Step 8 is run on the OUTPUT of step 6 (the merged GFF), not on the candidates
table. It can also be run standalone whenever the user has a merged GFF.

---

## Inputs

| Argument            | Flag              | Required | Description |
|---------------------|-------------------|----------|-------------|
| Merged GFF3         | `--merged_gff`    | yes      | Output of merge_split_genes.pl |
| Original GFF3       | `--orig_gff`      | yes      | Pre-merge annotation (to extract source gene protein seqs) |
| Genome FASTA        | `--genome_fa`     | yes      | Genome sequence (for gffread -y) |
| Reference proteome  | `--ref_fa`        | yes      | Same FASTA used in diamond step (subject) |
| Diamond DB          | `--diamond_db`    | yes*     | Pre-built diamond DB from ref_fa; built automatically if absent |
| Merge table         | `--merge_table`   | recommended | TSV from find_split_genes.pl / validate_merge_with_isoseq.pl; used to annotate source genes per merge_id |
| Output prefix       | `--out`           | yes      | Prefix for all output files |
| Threads             | `--threads`       | no       | Diamond / MAFFT threads (default: 4) |
| E-value             | `--evalue`        | no       | Diamond evalue cutoff (default: 1e-10, relaxed vs. detection step) |
| Skip MAFFT          | `--no_mafft`      | no       | Disable MSA scoring; report translation + diamond only |
| Min junction score  | `--min_junction`  | no       | Minimum per-junction MSA score to PASS (default: 0.5) |
| Min merged cov      | `--min_merged_cov`| no       | Minimum merged protein coverage by best ref hit (default: 0.60) |
| Min ref cov         | `--min_ref_cov`   | no       | Minimum ref protein coverage by merged protein (default: 0.50) |
| Keep MAFFT files    | `--keep_msa`      | no       | Write one .fa + .aln per merge to `<out>_msa/` |

*Diamond DB: if `--diamond_db` is absent, the script builds one in the workdir
using the `--ref_fa`.

---

## Outputs

### `<out>_translation_report.tsv`
One row per merged gene. Columns:

```
merge_id
new_gene_id
source_genes          (comma-separated)
merged_protein_len    (aa)
has_internal_stop     (0/1)
internal_stop_pos     (comma-separated positions, or ".")
merged_cov_by_ref     (fraction of merged protein covered by best ref hit)
ref_cov_by_merged     (fraction of ref protein covered by merged protein)
best_ref_hit          (ref protein ID)
best_ref_pident       (% identity to best ref hit)
ref_len               (aa length of best ref hit)
n_ref_hits            (number of ref hits used in MSA)
junction_scores       (one score per inter-gene junction, pipe-separated)
min_junction_score    (minimum of junction_scores)
mean_junction_score   (mean of junction_scores)
msa_flag              (GOOD_MSA / WEAK_JUNCTION / NO_HIT / SKIPPED)
translation_flag      (OK / FRAMESHIFT_DETECTED)
overall_result        (PASS / FAIL / REVIEW)
fail_reasons          (pipe-separated list of reasons, or ".")
```

### `<out>_pass.gff3`
Subset of `--merged_gff` containing only PASS merges.

### `<out>_fail.gff3`
Subset containing FAIL merges (for audit).

### `<out>_review.gff3`
Subset containing REVIEW merges (borderline — inspect manually).

### `<out>_pass_proteins.fa`
Translated protein sequences for PASS merges (from gffread).

### `<out>_msa/` (if `--keep_msa`)
One subdirectory per merge_id containing:
- `<merge_id>_input.fa`   — merged + source + ref proteins, raw
- `<merge_id>_aligned.fa` — MAFFT output (--auto --reorder)

---

## Algorithm

### Step A: Extract and translate merged CDS

```
gffread --merged_gff \
    -g genome_fa \
    -y <out>_merged_proteins.fa
```

Parse `<out>_merged_proteins.fa`:
- For each protein, scan for `*` (stop codon character that gffread inserts for
  internal stops) before the terminal position.
- gffread convention: terminal stop is included as `*` at the end; internal
  stops are also `*`. Detect by finding `*` before the last character.
- Record `has_internal_stop` and positions.
- Set `translation_flag = FRAMESHIFT_DETECTED` if any internal stop found.

**Implementation note on gffread output:**
gffread -y outputs the translated CDS. If a stop codon is in-frame internally,
it appears as `*` mid-sequence. The terminal `*` is the normal stop; strip it
before counting. Internal `*` characters are the red flags.

```perl
my $seq = $prot;
$seq =~ s/\*$//;          # strip terminal stop
my @ipos = ();
while ($seq =~ /\*/g) { push @ipos, pos($seq) }
my $has_internal = @ipos ? 1 : 0;
```

### Step B: Diamond blastp — merged proteins vs. reference

Run once for all merged proteins:

```
diamond blastp \
    -q <out>_merged_proteins.fa \
    -d <diamond_db> \
    -o <out>_merged_vs_ref.tsv \
    --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore scovhsp slen qcovhsp qlen \
    --evalue <evalue> \
    --max-target-seqs 5 \
    --threads <threads>
```

For each merged protein, keep up to 5 best hits (ranked by bitscore).

Compute:
```
merged_cov_by_ref = (qend - qstart + 1) / qlen      [best hit only, or max across hits]
ref_cov_by_merged = scovhsp / 100                    [from diamond field]
```

Use the best hit (highest bitscore) as `best_ref_hit`. Use up to
`--max_msa_refs` (default 3) hits in the MSA.

Flag as `NO_HIT` if no diamond hit found.

### Step C: Extract source gene proteins (for MSA)

From `--orig_gff` + genome_fa, use gffread to extract translated proteins for
each source gene in every merge group. These are stored keyed by gene_id.

```
gffread orig_gff \
    -g genome_fa \
    -y <out>_source_proteins.fa
```

Parse into a hash: `%source_prot{gene_id} = sequence`.

### Step D: Extract reference proteins for MSA

From `--ref_fa`, extract the sequences of the top N ref hits per merge.
Store in a hash: `%ref_prot{hit_id} = sequence`.

### Step E: MAFFT junction MSA scoring

For each merge group:

1. **Build input FASTA** containing:
   - The merged protein (labeled `MERGED_<merge_id>`)
   - Each source gene protein (labeled `SOURCE_<gene_id>`)
   - Top N ref proteins from diamond hits (labeled `REF_<hit_id>`)

2. **Run MAFFT**:
   ```
   mafft --auto --reorder --quiet input.fa > aligned.fa
   ```

3. **Find junction boundaries** in the merged protein alignment:
   - The merged CDS is a concatenation of source gene CDS in order.
   - Junction position J_k (in unaligned merged aa) = sum of CDS aa lengths of
     genes 1..k (where CDS length = CDS bp / 3, rounded down).
   - Map this to the alignment column by counting non-gap characters in the
     MERGED row up to J_k.

4. **Score each junction** (one score per adjacent pair of source genes):

   **Score = weighted combination of three sub-scores:**

   a. **Conservation score at junction window** (±W aa, default W=5):
      For each column in the window, compute the fraction of non-gap sequences
      that agree with the majority character. Average across the window.
      Range: 0.0–1.0. High = well-conserved region.

   b. **Gap pattern score**:
      In the alignment window around the junction, penalize columns where the
      MERGED row has gaps but the REF rows do not (indicates inserted/missing
      residues at the join point). Score = 1 - (merged_gap_cols / window_cols).

   c. **Ref continuity score**:
      In the reference proteins, does the sequence continue smoothly across
      the position that corresponds to J_k in the merged protein?
      Check: is the REF sequence gapless across the 10 columns spanning J_k?
      Score = fraction of ref sequences that are gapless across J_k ± 3 cols.

   **Combined junction score**:
   ```
   junction_score = 0.4 * conservation + 0.4 * ref_continuity + 0.2 * gap_pattern
   ```

   These weights are adjustable via `--w_conservation`, `--w_continuity`,
   `--w_gap` (must sum to 1.0).

5. **Assign msa_flag**:
   - `GOOD_MSA`       : min_junction_score >= threshold AND no internal stop
   - `WEAK_JUNCTION`  : min_junction_score < threshold (or 0 ref hits for MSA)
   - `NO_HIT`         : no diamond hit → no MSA possible
   - `SKIPPED`        : --no_mafft flag

### Step F: Overall PASS/FAIL/REVIEW assignment

```
FAIL  if translation_flag eq FRAMESHIFT_DETECTED
FAIL  if merged_cov_by_ref < min_merged_cov AND msa_flag ne GOOD_MSA
FAIL  if msa_flag eq NO_HIT AND merged_cov_by_ref < min_merged_cov
PASS  if translation_flag eq OK
        AND merged_cov_by_ref >= min_merged_cov
        AND ref_cov_by_merged >= min_ref_cov
        AND (msa_flag eq GOOD_MSA OR --no_mafft)
REVIEW otherwise  (passes some but not all thresholds — inspect manually)
```

The intent: FAIL is reserved for clear evidence of a bad merge. REVIEW
surfaces borderline cases for manual inspection without automatically discarding
them. PASS means the merge is supported at every level checked.

---

## Junction position calculation detail

This is the trickiest part. The merged CDS in the GFF is built from exons of
N source genes concatenated in genomic order. The boundary between gene k and
gene k+1 in the merged protein corresponds to:

```
J_k (aa) = floor( sum(CDS_bp for genes 1..k) / 3 )
```

This assumes in-frame joining (the current pipeline assumption). If the merge
was actually a mid-exon split, J_k will be slightly wrong — which is exactly
the kind of error that will produce a weak junction score or an internal stop.

To compute CDS_bp per source gene in the merged context:
- Parse the merged GFF CDS features; each will have been built from the source
  gene's exons (with possible UTR trimming by merge_split_genes.pl).
- Track cumulative CDS bp across all CDS features in strand order.
- Annotate junction positions at the cumulative totals.

**Alternatively (simpler):** use the `merged_from` attribute in the merged GFF
gene record (written by merge_split_genes.pl) to get the source gene IDs. Then
use the original GFF to get each source gene's CDS length. This doesn't require
parsing the merged CDS structure — just sum up the original source CDS lengths.

Use the simpler approach unless the merged GFF's CDS features differ materially
from the concatenation of source CDS features (e.g. due to UTR trimming at
terminals). For initial implementation, use the simpler approach and add a
comment that it assumes no CDS length change from the merge process.

---

## Data structures (Perl)

```perl
# Merge groups: keyed by new_gene_id parsed from merged GFF
%merge_info = (
    new_gene_id => {
        merge_id     => "...",
        source_genes => ["geneA", "geneB", ...],  # from merged_from attribute
        chrom        => "...",
        strand       => "+/-",
        cds_len_bp   => 0,    # filled during GFF parse
    }
)

# Proteins
%merged_prot  = (new_gene_id => "MKTLL...")
%source_prot  = (gene_id     => "MKTLL...")
%ref_prot     = (hit_id      => "MKTLL...")

# Diamond hits (per merged protein)
%diamond_hits = (new_gene_id => [
    { sseqid => "...", pident => 98.0, scovhsp => 90.5, slen => 450,
      qcovhsp => 95.0, qlen => 480, evalue => 1e-150, bitscore => 900 },
    ...
])

# Results
%results = (new_gene_id => {
    has_internal_stop  => 0,
    internal_stop_pos  => [],
    merged_cov_by_ref  => 0.92,
    ref_cov_by_merged  => 0.88,
    best_ref_hit       => "ref_protein_id",
    best_ref_pident    => 97.2,
    junction_scores    => [0.87, 0.91],   # one per junction
    msa_flag           => "GOOD_MSA",
    translation_flag   => "OK",
    overall_result     => "PASS",
    fail_reasons       => [],
})
```

---

## External tool dependencies

| Tool       | Purpose                    | Required? |
|------------|----------------------------|-----------|
| gffread    | CDS translation            | yes       |
| diamond    | merged vs ref blastp       | yes       |
| mafft      | multiple sequence alignment| no (skip with --no_mafft) |
| makeblastdb| (not needed — use diamond) | no        |

Tools are expected in `$PATH` unless overridden:
- `--gffread_bin`
- `--diamond_bin`
- `--mafft_bin`

---

## run_mender.pl integration (future)

Add as step 8 to run_mender.pl. Config keys in `[validation]` section:

```ini
run_translation_validation = yes
min_junction_score         = 0.5
min_merged_cov             = 0.60
min_ref_cov                = 0.50
keep_msa                   = no
no_mafft                   = no
```

The script can also be run standalone by the user after a manual merge or after
running only a subset of pipeline steps.

---

## Edge cases and notes

1. **Single-exon source genes**: junction at position 0 (N-terminal boundary) or
   at C-terminus — these are boundary junctions and should not be scored as
   internal junctions. The merge of a single-exon gene to a multi-exon gene has
   one junction; score it.

2. **Strand**: gffread handles strand correctly when translating. No special
   handling needed for minus-strand genes.

3. **Partial CDS**: some source genes may have incomplete CDS (no start or no
   stop codon annotated). gffread still translates what is there. Internal stops
   from partial CDS at the N-terminal boundary (missing start codon means phase
   may be off) may produce false FRAMESHIFT_DETECTED. Consider: only flag
   internal stops that are more than 10aa from either end of the merged protein,
   or annotate whether the source genes had partial CDS in the original GFF.

4. **MULTI_ISOFORM_JOIN**: merged GFF may contain multiple mRNA features per
   gene. gffread will produce one protein per mRNA. Need to pick a representative
   protein per merged gene (longest non-stop CDS, or the mRNA with the most
   exons). Score each and report the best, noting how many transcripts were
   evaluated.

5. **Very long proteins** (LARGE_SPAN merges): MAFFT --auto handles well up to
   a few thousand aa. For proteins > 5000 aa consider using `--globalpair
   --maxiterate 0` for speed, or skip MSA and flag as REVIEW.

6. **No ref hit (NO_HIT)**: possible for highly diverged genes, ORFan genes, or
   annotation errors. These should be REVIEW, not automatically FAIL — the
   absence of a hit does not prove the merge is wrong, just unverifiable.

7. **Diamond version compatibility**: use the same outfmt as run_mender.pl step 2.
   The script will build its own db if needed; users should reuse the pre-built
   db from the pipeline workdir to save time.

---

## Implementation order (recommended)

1. Arg parsing + config (Getopt::Long)
2. Parse merged GFF → %merge_info (pull `merged_from` and CDS bp per source gene)
3. gffread translate → %merged_prot; internal stop scan
4. gffread translate source genes → %source_prot
5. Diamond run → %diamond_hits; compute merged_cov_by_ref / ref_cov_by_merged
6. Extract %ref_prot sequences for top N hits
7. MAFFT loop: build input FASTA, run, parse aligned FASTA, compute junction scores
8. PASS/FAIL/REVIEW assignment
9. Write TSV report + split GFF outputs + pass proteins FASTA
10. Write --keep_msa alignment files if requested
