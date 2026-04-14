# Mender pipeline diagram draft

This is a standalone draft so the structure can be reviewed before deciding
whether to embed it in `README.md`.

```mermaid
flowchart TD

    I["**Inputs**<br/>Genome annotation GFF3<br/>Query proteome FASTA<br/>Reference proteome FASTA<br/>Optional IsoSeq GFF<br/>Genome FASTA for translation validation"]

    P1["**1. Prepare inputs**<br/>Clean query proteins<br/>Extract gene features<br/>Extract IsoSeq mRNA features if provided"]
    P2["**2. Protein homology search**<br/>DIAMOND blastp of query proteins against a well-annotated reference proteome"]
    P3["**3. IsoSeq overlap mapping**<br/>Use bedtools to find long-read transcripts that overlap annotated genes<br/>Skipped if no IsoSeq GFF is provided"]
    P4["**4. Detect split-gene candidates**<br/>Find fragmentary homologs<br/>Test whether adjacent genes tile the same reference protein<br/>Chain supported pairs into multi-gene loci"]
    P5["**5. Add long-read transcript support**<br/>Use overlap data to count IsoSeq reads spanning 2 or more genes in a candidate locus<br/>Classify support as FULL_SPAN, PARTIAL_SPAN, or none"]
    P6["**6. Merge selected candidates**<br/>Rebuild a single gene model in genomic order<br/>Join exon, CDS, and UTR structure<br/>Recalculate CDS phase<br/>Assign new gene and transcript IDs"]
    P7["**7. GT GFF3 check**<br/>Fast syntax and spec check on Mender-created genes"]
    P8["**8. Translation validation**<br/>Translate merged CDS<br/>Compare merged proteins to reference proteins<br/>Score each inferred junction<br/>Assign PASS, FAIL, or REVIEW"]
    P9["**9. AGAT gene-model check**<br/>Check whether merged annotations remain structurally coherent as gene models"]
    O["**Outputs**<br/>Merged GFF3<br/>Removed source genes GFF3<br/>Optional translation report and PASS/FAIL/REVIEW GFFs<br/>Final validated GFF with FAIL merges restored when translation validation is run"]

    B1["**Biological signal from protein homology**<br/>A true split gene produces 2 or more adjacent fragments that hit the same orthologous protein<br/>Their alignments should cover different, non-overlapping parts of that protein<br/>Fragments that tile end-to-end suggest one biological gene was annotated as multiple models"]
    B2["**Biological signal from IsoSeq**<br/>A spanning long-read transcript supports co-transcription across adjacent fragments<br/>FULL_SPAN supports the whole candidate chain<br/>PARTIAL_SPAN suggests a terminal fragment may not belong"]
    F1["**Merge gate and review filters**<br/>Apply skip_flags, flags, min_tiling, min_cov, require_isoseq, and isoseq_min_spanning<br/>fix_partial can trim unsupported terminal genes in PARTIAL_SPAN cases<br/>asym_trim controls whether weak terminal or internal joins are automatically trimmed or split earlier"]
    F2["**Candidate annotation groups**<br/>Protein evidence strength: CLEAN, STRONG, SINGLE_HIT, LOW_COV, WEAK_END, WEAK_INTERNAL<br/>Structural caution: SKIPPED_GENE, TRANSITIVE_JOIN, MULTI_ISOFORM_JOIN, LARGE_SPAN, LARGE_SPAN_EXTREME<br/>IsoSeq support: FULL_SPAN, PARTIAL_SPAN, none"]
    OPT["**Pipeline options**<br/>Steps 3 and 5 are skipped if no IsoSeq GFF is provided<br/>Translation validation is skipped if run_translation_validation = no<br/>AGAT is skipped if run_agat = no<br/>--steps runs only selected numbered stages<br/>--dry_run prints commands without executing them"]

    I --> P1 --> P2
    P2 -->|IsoSeq provided| P3
    P2 -->|No IsoSeq| P4
    P3 --> P4
    P4 -->|IsoSeq provided| P5
    P4 -->|No IsoSeq| P6
    P5 --> P6 --> P7 --> P8 --> P9 --> O
    P7 -->|Translation validation skipped| P9
    P7 -->|Translation and AGAT skipped| O
    P8 -->|AGAT skipped| O

    B1 -. biological rationale .-> P4
    B2 -. transcript evidence .-> P5
    F1 -. selection logic .-> P6
    F2 -. summarized flags .-> P6
    OPT -. optional branches .-> P3
    OPT -. optional branches .-> P5
    OPT -. optional branches .-> P8
    OPT -. optional branches .-> P9
```

## Design notes

- The main vertical flow shows the executable pipeline in `run_mender.pl`.
- The right-hand callouts separate **biology**, **filters**, and **options** so
  the core process is still readable.
- Flags are grouped by meaning rather than shown as many separate arrows.

## If this moves into the README

- Best placement: just before `## Pipeline Steps`.
- If it feels too dense there, split it into:
  1. one main pipeline flow
  2. one smaller biological evidence panel
  3. one compact options and flags legend
