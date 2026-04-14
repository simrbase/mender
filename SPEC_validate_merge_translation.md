# SPEC: validate_merge_translation.pl

## Purpose

Post-merge translation quality validator. After `merge_split_genes.pl` produces a
merged GFF, this script validates each merged gene model by translating its CDS,
running diamond blastp against a reference proteome, aligning merged + source +
reference proteins with MAFFT, and scoring the alignment at each gene-gene junction.

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
Step 7  gt               gt gff3validator on Mender-only GFF [optional]
        ↓
Step 8  validate_transl  validate_merge_translation.pl   ← THIS SCRIPT
        ↓
Step 9  agat             agat on pass GFF only           [optional]
```

Step 8 runs on the output of step 6 (the merged GFF). It can also be run
standalone whenever the user has a merged GFF.

---

## Internal Steps (A–L)

| Step | Name | Description |
|------|------|-------------|
| A | Parse merged GFF | Collect new gene IDs, `merged_from` source gene lists, GFF line blocks |
| B | Parse original GFF | Find representative transcript and CDS length (bp) per source gene |
| C | Translate merged proteins | gffread -y on merged GFF; pick representative per gene (prefer stop-free, then longest); write clean FASTA for diamond (`.` → `X`, `*` removed) |
| D | Translate source proteins | gffread -y on original GFF; map to source gene IDs |
| E | Diamond database | Build or verify diamond DB from ref_fa |
| F | Diamond blastp (ref) | Merged proteins vs reference proteome; compute bidirectional coverage metrics |
| G | Load ref sequences | Load reference proteome sequences into memory for MSA fallback |
| H | SwissProt blast (MSA) | If `--swissprot_fa` provided: diamond blast vs SwissProt, load hit sequences; these are used as MSA reference sequences instead of ref_fa (manually curated, no partial seqs, no `.` characters) |
| I | Load merge table | Load merge flags from find/isoseq step for TSV annotation |
| J | MAFFT junction scoring | Align merged + source + ref proteins; locate junction positions; score conservation, gap pattern, ref continuity at each junction |
| K | PASS / FAIL / REVIEW | Assign overall result per merge from translation flag + coverage + MSA flag |
| L | Write outputs | TSV report, pass/fail/review GFF3s, pass proteins FASTA |

---

## Inputs

| Argument            | Flag              | Required | Description |
|---------------------|-------------------|----------|-------------|
| Merged GFF3         | `--merged_gff`    | yes      | Output of merge_split_genes.pl |
| Original GFF3       | `--orig_gff`      | yes      | Pre-merge annotation (source gene proteins + CDS lengths) |
| Genome FASTA        | `--genome_fa`     | yes      | Genome sequence (for gffread -y) |
| Reference proteome  | `--ref_fa`        | yes      | Same FASTA used in pipeline diamond step (subject) |
| Diamond DB          | `--diamond_db`    | no       | Pre-built diamond DB from ref_fa; built automatically if absent |
| SwissProt FASTA     | `--swissprot_fa`  | no       | UniProt/SwissProt FASTA for MSA reference sequences; falls back to ref_fa if absent |
| SwissProt DB        | `--swissprot_db`  | no       | Pre-built diamond DB from swissprot_fa; built automatically if absent |
| Merge table         | `--merge_table`   | recommended | TSV from find_split_genes.pl / validate_merge_with_isoseq.pl; annotates source flags in report |
| Output prefix       | `--out`           | yes      | Prefix for all output files |
| Threads             | `--threads`       | no       | Diamond / MAFFT threads (default: 4) |
| E-value             | `--evalue`        | no       | Diamond evalue cutoff (default: 1e-10) |
| Aligner             | `--aligner`       | no       | MSA tool: `mafft_fast` (default) \| `mafft_auto` \| `kalign3`. kalign3 is ~8x faster and recommended. |
| Kalign3 path        | `--kalign_bin`    | no       | Path to kalign binary (default: `kalign` in PATH) |
| Skip MSA            | `--no_msa`        | no       | Disable MSA scoring; report translation + diamond only |
| Min junction score  | `--min_junction`  | no       | Minimum per-junction MSA score to PASS (default: 0.5) |
| Min merged cov      | `--min_merged_cov`| no       | Minimum merged protein coverage by best ref hit (default: 0.60) |
| Min ref cov         | `--min_ref_cov`   | no       | Minimum ref protein coverage by merged protein (default: 0.50) |
| Max MSA refs        | `--max_msa_refs`  | no       | Max ref hits to include per MSA (default: 3) |
| Keep MAFFT files    | `--keep_msa`      | no       | Write per-merge .fa + .aln to `<out>_msa/` |
| Score weights       | `--w_conservation`, `--w_continuity`, `--w_gap` | no | MSA sub-score weights (must sum to 1.0; defaults: 0.4, 0.4, 0.2) |

To download SwissProt:
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
diamond makedb --in uniprot_sprot.fasta --db uniprot_sprot.dmnd
```

---

## Outputs

### `<out>_report.tsv`
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
source_flags          (merge flag from find/isoseq step, e.g. STRONG, WEAK_END)
```

### `<out>_pass.gff3`
Subset of `--merged_gff` containing only PASS merges. Each gene feature has
three additional GFF9 attributes added: `transl_result=PASS`,
`msa_flag=<flag>`, `min_junction_score=<score>`.

### `<out>_fail.gff3`
Subset containing FAIL merges (audit trail). Same added attributes as pass GFF.

### `<out>_review.gff3`
Subset containing REVIEW merges (borderline — inspect manually). Same added attributes.

### `<out>_pass_proteins.fa`
Translated protein sequences for PASS merges.

### `<out>_msa/` (if `--keep_msa`)
One pair of files per merge_id:
- `<merge_id>_input.fa`   — merged + source + ref proteins (MERGED_, SOURCE_, REF_, REFPROT_), unaligned
- `<merge_id>_aligned.fa` — aligner output

---

## Algorithm

### Step A: Parse merged GFF

Read the merged GFF produced by `merge_split_genes.pl`. For each gene feature
with a `merged_from` attribute, record:
- `new_gene_id`: the `ID` attribute
- `source_genes`: the comma-separated list from `merged_from`
- `merge_id`: the `merge_id` attribute
- All GFF lines belonging to that gene (for output splitting in Step L)

### Step B: Parse original GFF for source CDS lengths

For each source gene, identify the representative transcript (longest total CDS
length in bp) and record that length. These CDS lengths are used in Step J to
locate junction positions in the merged protein.

### Step C: Translate merged proteins (gffread)

```
gffread merged_gff -g genome_fa -y <out>_merged_proteins.fa
```

For each translated protein:
- Strip terminal stop (`*`)
- Detect internal stops: any `*` remaining after stripping terminal
- Pick representative per gene: prefer stop-free, then longest
- Record `has_internal_stop` and positions

Write a clean per-gene FASTA for diamond:
- One sequence per gene (representative)
- Internal stops (`*`) removed
- Ambiguous AA (`.`, emitted by gffread for N-containing codons) replaced with `X`
  (diamond rejects `.` but accepts `X`)
- Header is the gene ID (so diamond hits come back keyed by gene ID directly)

### Step D: Translate source proteins (gffread)

```
gffread orig_gff -g genome_fa -y <out>_source_proteins.fa
```

Map each protein to its source gene ID via the original GFF's mRNA→gene
parent relationships. Keep only the representative transcript per gene.

### Step E: Diamond database

Use the pre-built `--diamond_db` if provided. Otherwise build from `--ref_fa`:
```
diamond makedb --in ref_fa --db <out>_ref.dmnd
```

### Step F: Diamond blastp (merged proteins vs reference)

```
diamond blastp \
    --db <diamond_db> \
    --query <out>_merged_proteins_clean.fa \
    --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore scovhsp slen qcovhsp qlen \
    --evalue <evalue> \
    --max-target-seqs 5 \
    --threads <threads>
```

For each merged gene, keep up to 5 best hits (ranked by bitscore, deduplicated
by subject ID). Compute coverage metrics from the best hit:

```
merged_cov_by_ref = (qend - qstart + 1) / qlen
ref_cov_by_merged = scovhsp / 100
```

These are used for PASS/FAIL thresholds in Step K. Coverage is assessed against
the reference proteome regardless of whether SwissProt is used for MSA.

### Step G: Load reference proteome sequences

Read all sequences from `--ref_fa` into memory. Used as MSA reference sequences
if SwissProt is not configured (Step H fallback).

### Step H: SwissProt diamond blast (MSA reference sequences)

If `--swissprot_fa` is provided:

1. Build or verify SwissProt diamond DB:
   ```
   diamond makedb --in swissprot_fa --db <out>_swissprot.dmnd
   ```
2. Diamond blast merged proteins vs SwissProt:
   ```
   diamond blastp --db <swissprot_db> --query <out>_merged_proteins_clean.fa ...
   ```
3. Parse hits into `%sprot_hits` (top 5 per gene by bitscore, deduplicated)
4. Load only the needed SwissProt sequences into `%sprot_prot`

**Why SwissProt for MSA?** SwissProt is manually curated: no partial sequences,
no annotation errors, no `.` characters from N-containing codons. This gives
cleaner alignments and more reliable junction scores than the species reference
proteome. Coverage metrics (Step F) still use the species reference proteome.

If `--swissprot_fa` is absent, Step H is skipped and the MAFFT loop uses
`%ref_prot` (from Step G) as MSA reference sequences.

### Step I: Load merge table

If `--merge_table` is provided, load the `merge_id → flag` mapping from the
find/isoseq TSV. Used only for annotation in the output report (Step L).

### Step J: MAFFT junction scoring

For each merged gene:

1. **Determine MSA reference hits and sequences:**
   - If SwissProt available: use `%sprot_hits{gid}` and `%sprot_prot`
   - Otherwise: use `%gene_diamond_hits{gid}` and `%ref_prot`

2. **Skip conditions:**
   - No ref diamond hit (`%gene_diamond_hits` empty): flag `NO_HIT`, skip MSA
   - `--no_msa`: flag `SKIPPED`, skip MSA
   - No merged protein sequence: flag `WEAK_JUNCTION`, skip

3. **Build MSA input FASTA** (cleaned sequences, `.` → `X`, `*` removed):
   - `MERGED_<gid>`: merged protein
   - `SOURCE_<gid>`: each source gene protein
   - `REF_<id>`: top N SwissProt hits (or ref_fa hits if SwissProt not set), up to `--max_msa_refs`
   - `REFPROT_<id>`: best ref_fa hit (the protein that originally triggered the merge detection); included when SwissProt is also providing REF_ sequences

4. **Run MSA aligner** (set via `--aligner`):
   - `kalign3` (recommended): `kalign -i input.fa -o aligned.fa`
   - `mafft_fast` (default): `mafft --retree 2 --maxiterate 0 --reorder --quiet input.fa`
   - `mafft_auto`: `mafft --auto --reorder --quiet input.fa`
     (proteins > 5000 aa: `--retree 1` added for speed)

5. **Locate junction positions in the alignment:**

   Junction J_k is the boundary between source gene k and gene k+1 in the merged
   protein. In unaligned aa coordinates:
   ```
   J_k (aa, 1-based) = floor( sum(CDS_bp for source genes 1..k) / 3 )
   ```
   Map to alignment column by counting non-gap characters in the MERGED row up
   to position J_k.

6. **Score each junction** in a window of ±5 columns:

   **a. Conservation score:**
   For each column in the window, compute the fraction of non-gap sequences
   that agree with the majority character. Average across the window.
   Range: 0.0–1.0.

   **b. Gap pattern score:**
   Penalize columns where MERGED has a gap but REF sequences do not.
   ```
   gap_pattern = 1 - (merged_gap_cols / valid_cols)
   ```

   **c. Ref continuity score:**
   Fraction of REF sequences that are gapless across J_k ± 3 columns.
   Range: 0.0–1.0. High = the reference protein continues smoothly at this
   position, suggesting a real protein boundary exists here.

   **Combined junction score:**
   ```
   junction_score = 0.4 * conservation + 0.4 * ref_continuity + 0.2 * gap_pattern
   ```
   Weights adjustable via `--w_conservation`, `--w_continuity`, `--w_gap`
   (must sum to 1.0).

7. **Assign msa_flag:**
   - `GOOD_MSA`:      min junction score >= `--min_junction` AND no internal stop
   - `WEAK_JUNCTION`: min score below threshold, or no ref sequences in alignment
   - `NO_HIT`:        no diamond hit → MSA not possible
   - `SKIPPED`:       `--no_msa` flag set

### Step K: PASS / FAIL / REVIEW assignment

```
FAIL    if has_internal_stop
FAIL    if merged_cov_by_ref < min_merged_cov AND msa_flag != GOOD_MSA
FAIL    if msa_flag == NO_HIT AND merged_cov_by_ref < min_merged_cov

PASS    if translation_flag == OK
          AND merged_cov_by_ref >= min_merged_cov
          AND ref_cov_by_merged >= min_ref_cov
          AND (msa_flag == GOOD_MSA OR --no_msa)

REVIEW  otherwise  (passes some but not all thresholds)
```

FAIL = clear evidence of a bad merge. REVIEW = borderline, inspect manually.
PASS = supported at every level checked.

### Step L: Write outputs

- `<out>_report.tsv`: one row per merged gene
- `<out>_pass.gff3` / `_fail.gff3` / `_review.gff3`: GFF3 subsets
- `<out>_pass_proteins.fa`: translated proteins for PASS merges
- `<out>_msa/` (if `--keep_msa`): per-merge alignment files

---

## Junction position calculation detail

The merged CDS is a concatenation of source gene CDS in genomic order. The
boundary between gene k and gene k+1 corresponds to:

```
J_k (aa) = floor( sum(CDS_bp for source genes 1..k) / 3 )
```

CDS lengths come from the **original GFF** (Step B), not the merged GFF. This
assumes the merge did not change CDS lengths (valid for the current
merge_split_genes.pl implementation, which concatenates exons without editing
internal boundaries). If a terminal exon was trimmed during a PARTIAL_SPAN fix,
J_k for the last junction may be slightly off.

---

## External tool dependencies

| Tool    | Purpose                     | Required? |
|---------|-----------------------------|-----------|
| gffread | CDS translation             | yes       |
| diamond | merged vs ref/sprot blastp  | yes       |
| mafft   | MSA junction scoring        | only if `--aligner mafft_fast` or `mafft_auto` and `--no_msa` not set |
| kalign3 | MSA junction scoring        | only if `--aligner kalign3` and `--no_msa` not set |

Tools expected in `$PATH` unless overridden via `--gffread_bin`, `--diamond_bin`,
`--mafft_bin`, `--kalign_bin`.

---

## run_mender.pl integration

Configured via `[validation]` section of the pipeline config file:

```ini
run_translation_validation = yes
min_junction_score         = 0.5
min_merged_cov             = 0.60
min_ref_cov                = 0.50
aligner                    = kalign3   # kalign3 | mafft_fast | mafft_auto
keep_msa                   = no
no_msa                     = no
transl_out_prefix          =           # blank = <workdir>/transl
swissprot_fa               = /path/to/uniprot_sprot.fasta
swissprot_db               =           # blank = built automatically
final_gff_include_review   = yes       # yes = keep REVIEW in validated GFF; no = replace with source genes
```

The script can also be run standalone outside of run_mender.pl.

---

## Edge cases and notes

1. **Sequence cleaning for MAFFT**: gffread emits `.` for codons containing N
   bases. Diamond rejects `.` but MAFFT also handles it poorly. All sequences
   passed to MAFFT have `.` → `X` and `*` removed. The raw sequences in
   `%merged_prot` / `%source_prot` are kept unmodified for internal stop
   detection.

2. **SwissProt hit IDs**: SwissProt FASTA headers have the format
   `>sp|P12345|GENE_HUMAN`. The first token (up to the first space) is used as
   the hit ID, which becomes the key in `%sprot_prot`. Diamond reports the same
   token in field 1 of the blastp output.

3. **Strand**: gffread handles strand correctly. No special handling needed for
   minus-strand genes.

4. **MULTI_ISOFORM_JOIN merges**: merged GFF may contain multiple mRNA per gene.
   gffread produces one protein per mRNA. The representative is chosen as the
   longest stop-free translation (or longest overall if all have stops).

5. **Very long proteins** (LARGE_SPAN merges): MAFFT --auto handles up to a few
   thousand aa well. For proteins > 5000 aa, `--retree 1` is added automatically
   for speed.

6. **NO_HIT**: absence of a diamond hit does not prove the merge is wrong —
   highly diverged or ORFan genes will have no hit. These are flagged REVIEW
   rather than FAIL.

7. **Phase errors in gt (Step 7)**: a wrong phase value in the GFF annotation
   does not affect gffread translation (gffread uses coordinates, not phase).
   A merge-introduced phase error will not produce an internal stop unless the
   exon coordinates themselves are wrong.
