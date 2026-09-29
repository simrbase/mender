#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## find_split_genes.pl
##
## Identifies gene annotations that may have been erroneously split
## during genome annotation. The core idea: if two or more adjacent
## genes each cover a non-overlapping portion of the same reference protein,
## they are candidates to be merged into a single gene.
##
## Evidence comes from diamond blastp of the annotation's proteins vs a well-annotated
## reference proteome (e.g. Anolis). We look for adjacent gene pairs where:
##   1. Both genes hit the same reference protein as fragments (scovhsp <= 85%)
##   2. Their alignments tile on the reference protein (gap <= wiggle aa)
##   3. They are within max_dist genes of each other on the same chromosome
##
## diamond blastp --db $SUBJECT --query multi.fa \
##   --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend \
##             sstart send evalue bitscore scovhsp slen qcovhsp qlen \
##   --threads 4 --evalue 1e-20 --out diamond.out
##
## qcovhsp (query coverage per HSP) and qlen (query sequence length) are
## required as of this version. They extend the fragment filter to catch
## split fragments of shorter proteins that would otherwise be excluded by
## the scovhsp-only threshold. See FRAGMENT FILTER below.
##
## Usage: perl find_split_genes.pl [--no_asym_trim] \
##                                  [--large_span_warn N] \
##                                  [--large_span_extreme N] \
##                                  [--max_dist N] [--wiggle N] [--asym_threshold N] \
##                                  [--low_cov_thresh N] [--max_evalue E] \
##                                  [--frag_min_qcov N] [--ref_min_qcov N] [--end_tol N] \
##                                  [--overlap_min_frac F] [--overlap_min_qcov N] \
##                                  diamond.out genes.gff subject.fa query.fa
##
##   --frag_min_qcov N     a fragment hit must align over >= N% of the piece itself
##                         (default 30; 0 = off). 60 was the old default: it also drops
##                         real splits whose piece carries sequence the reference lacks
##                         (Chamaeleo: 114 merges more at 30, Congeria 24; see README)
##   --ref_min_qcov N      rebuild_candidates: references every chain gene aligns over
##                         >= N% of itself rank before the others (default 60). With a
##                         lower frag_min_qcov, references that share only a domain with
##                         a piece would otherwise rank first on identity and give the
##                         rebuild a domain protein's shape
##   --overlap_min_frac F  PIECES_OVERLAP: two pieces' alignments on one reference
##                         overlap by >= F of the shorter alignment (default 0.5)
##   --overlap_min_qcov N  PIECES_OVERLAP: each piece aligns over >= N% of itself
##                         on that reference (default 60)
##
## Output files:
##
## split_genes_summary.txt
##   One row per gene pair that passes all filters. Best supporting hit only.
##   Columns:
##     g1, g2             : the two gene IDs (genomic order)
##     g1_coord, g2_coord : chromosome.start.end for each gene
##     genomic_dist       : number of rank positions between g1 and g2
##                          (1 = directly adjacent, 2 = one gene between them, etc)
##     g1_desc, g2_desc   : descriptions from query fasta headers
##     tiling_ref         : reference protein with highest combined coverage
##                          pct for this specific gene pair. Selected purely
##                          by combined span / subject length — no description
##                          preference applied here.
##     tiling_ref_desc    : description of the tiling_ref protein
##     tiling_cov_pct     : % of tiling_ref covered by g1+g2 alignments together
##                          (span from min sstart to max send / subject length)
##     tiling_gap         : aa between the two alignments on tiling_ref
##                          positive = gap between alignments (they don't quite meet)
##                          zero     = alignments meet perfectly end to end
##                          negative = alignments overlap slightly
##     g1_subj, g2_subj   : alignment coordinates on the subject (sstart-send)
##     g1_cov, g2_cov     : individual scovhsp% for each gene vs tiling_ref
##     g1_evalue,g2_evalue: diamond evalues for each gene vs tiling_ref
##     num_tiling_hits    : how many different reference proteins support this
##                          pair as tiling fragments — higher = more confident
##
## split_genes_detail.txt
##   One row per supporting tiling hit per gene pair, sorted by coverage desc.
##   Same columns as summary but for every tiling hit, not just the best.
##   Use this to see the full evidence for any pair from the summary.
##   Columns: g1, g2, genomic_dist, tiling_ref, tiling_ref_desc, tiling_cov_pct,
##            tiling_gap, g1_subj, g2_subj, g1_cov, g2_cov,
##            g1_evalue, g2_evalue
##
## split_genes_merge.txt
##   One row per merge candidate. Adjacent pairs sharing a gene are chained
##   into multi-gene groups (e.g. A-B + B-C becomes A,B,C in one row).
##   Chains are trimmed and split based on junction tiling hit asymmetry
##   before writing — see CHAIN TRIMMING AND SPLITTING below.
##   Two kinds of reference protein are reported, one per step:
##     tiling_ref         : DETECTION. The reference the pieces tile on with the
##                          highest combined coverage (for a chain: that of its best
##                          pair). tiling_cov_pct, tiling_gap and evalues are all
##                          measured on it, and LOW_COV is judged on it.
##     rebuild_candidates : REBUILD. The reference proteins hitting two or more chain
##                          genes, ranked: most chain genes first, then identity-
##                          weighted coverage (the fraction of the reference's
##                          residues the fragments align identically: sum(pident x
##                          aligned span) / length, so a close reference beats a
##                          distant one at the same coverage), then a real
##                          description over a bare ID, then the id. The miniprot
##                          rebuild tries them in order; the one it uses is reported
##                          as rebuild_ref in miniprot.tsv and report.tsv.
##   tiling_ref is usually one of the candidates (unless 10 others rank above it) but
##   often not the first: the closest
##   reference (often the species' own transcript) rebuilds best, the best-tiling one
##   is the evidence that the pieces are consecutive parts of one protein.
##   Columns:
##     merge_id              : unique identifier
##     num_genes             : how many genes to merge
##     genes_in_order        : comma-sep gene IDs in genomic order
##     coords_in_order       : comma-sep coordinates in genomic order
##     ref                   : chromosome
##     gene_descs            : comma-sep descriptions in genomic order
##     tiling_ref            : reference protein the best pair tiles on (see above)
##     tiling_ref_desc       : description of tiling_ref
##     tiling_cov_pct        : combined coverage of tiling_ref by the best pair
##     tiling_gap            : tiling gap of the best pair on tiling_ref
##     evalues               : comma-sep evalues in genomic order vs tiling_ref
##                             '?' means that gene has no fragment hit vs tiling_ref
##     junction_tiling_hits  : pipe-sep tiling hit counts per junction in order
##                             e.g. "22|8" for a 3-gene chain with two junctions
##                             only present for chains with 2+ pairs, else NA
##     junction_genomic_dist : pipe-sep genomic rank distances per junction
##                             e.g. "1|2" means first junction is direct neighbors,
##                             second junction skips one gene
##                             only present for chains with 2+ pairs, else NA
##     skipped_genes         : comma-sep gene IDs of any genes that sit between
##                             a junction pair (genomic_dist > 1). These genes
##                             were not flagged as merge candidates but sit
##                             inside your merge locus — worth checking if they
##                             are also split gene fragments that failed filters.
##                             'none' if all junctions are directly adjacent.
##     max_tiling_hits       : highest num_tiling_hits seen across any pair
##     num_pairs_in_chain    : number of pairwise tiling relationships in chain
##     flag                  : quality flag(s), comma-sep if multiple apply
##                             see FLAG DEFINITIONS below
##     rebuild_candidates    : comma-sep reference proteins, ranked (see above, at
##                             most 10); the miniprot rebuild tries them in order
##     overlap_ref           : the reference on which two pieces align to the same
##                             region (PIECES_OVERLAP); '.' if none
##
## =============================================================================
## KEY CONCEPTS
## =============================================================================
##
## TILING_COV_PCT
##   Brief: What fraction of the reference protein is spanned by the two genes
##          together.
##
##   Detail: This is NOT the sum of the two individual coverages. It is the
##           span from the leftmost alignment start to the rightmost alignment
##           end on the reference protein, divided by the reference protein
##           length. Any gap between the two alignments is included in the span,
##           so tiling_cov_pct slightly overstates true coverage when there
##           is a gap. For a 15aa gap on a 500aa protein this is ~3% inflation.
##
##   Diagram:
##
##     reference protein (500aa):
##     |--------------------------------------------------|
##
##     gene1 aligns to positions 10-200:
##     |------gene1------|
##
##     gene2 aligns to positions 205-480:
##                        |----------gene2-----------|
##
##     combined span = 480 - 10 + 1 = 471aa
##     tiling_cov_pct = 471 / 500 * 100 = 94.2%
##     tiling_gap = 205 - 200 = 5aa (positive = gap between alignments)
##
##   Contrast with summing individual coverages:
##     gene1 scovhsp = 191/500 = 38.2%
##     gene2 scovhsp = 276/500 = 55.2%
##     naive sum     = 93.4%  (close but not same, would double-count overlaps)
##
## -----------------------------------------------------------------------------
##
## TILING_GAP
##   Brief: How well the two alignments meet in the middle on the reference
##          protein.
##
##   Detail: tiling_gap = sstart_of_rightmost - send_of_leftmost
##           Positive value  = gap between alignments (they don't quite meet)
##           Zero            = alignments meet perfectly end to end
##           Negative value  = alignments overlap slightly
##           We allow abs(tiling_gap) <= wiggle (default 15aa) to account for
##           alignment fuzziness at the ends of HSPs.
##
##   Diagram:
##
##     Perfect join (gap = 0):
##     |------gene1------|gene2------------|
##                       ^ meet exactly
##
##     Small gap (gap = 5, within wiggle=15, PASSES):
##     |------gene1------| - - |gene2------|
##                        5aa gap
##
##     Small overlap (gap = -3, within wiggle=15, PASSES):
##     |------gene1--------|
##                      |gene2------------|
##                      3aa overlap
##
##     Too large a gap (gap = 40, exceeds wiggle=15, REJECTED):
##     |------gene1------|                        |gene2---|
##                        <-------- 40aa -------->
##
##     Too large an overlap (gap = -20, exceeds wiggle=15, REJECTED):
##     |------gene1-----------------------------|
##                   |gene2--------------------|
##                   <------ 20aa overlap ----->
##
## -----------------------------------------------------------------------------
##
## NUM_TILING_HITS
##   Brief: How many different reference proteins independently support this
##          pair as tiling fragments. Higher = more confident.
##
##   Detail: For a reference protein to contribute a tiling hit, BOTH gene1
##           AND gene2 must have a fragment-level blast hit against that same
##           protein, AND those two hits must tile end-to-end on it. The set
##           of candidate proteins is the intersection of each gene's hit list,
##           so only proteins that both genes independently matched are tested.
##
##           High num_tiling_hits arises when the reference proteome contains
##           multiple proteins (paralogs or multi-species orthologs) that are
##           each similar enough to BOTH fragments simultaneously — meaning
##           the split cuts at a structurally conserved boundary that is
##           shared across all those proteins (e.g., between a ligand-binding
##           and a signaling domain that is conserved across a whole family).
##
##           Low num_tiling_hits does NOT mean false positive. It commonly
##           occurs for single-copy genes that have few representatives in
##           the reference proteome, or when the split site falls in a region
##           that has diverged across paralogs so that different family members
##           look more like one half than the other. Such cases are still real
##           splits — they simply lack the breadth of reference coverage to
##           drive counts higher.
##
##           An important blind spot: if the two halves have diverged enough
##           to preferentially hit DIFFERENT reference proteins (no overlap in
##           their hit lists), the pair will never appear in output at all —
##           the blast intersection is empty. IsoSeq spanning reads are the
##           primary way to catch these.
##
##           1 hit  = thin protein evidence; cross-check with IsoSeq
##           2-4    = moderate evidence
##           5+     = strong; the split boundary is conserved across many
##                    reference proteins
##           10+    = very strong
##
##   Diagram:
##
##     reference protein A:  |--gene1--|--gene2--|   gap=2   <- tiling hit
##     reference protein B:  |--gene1--|--gene2--|   gap=-1  <- tiling hit
##     reference protein C:  |--gene1--|--gene2--|   gap=5   <- tiling hit
##     reference protein D:  |--gene1-----------| (gene2 gap=80, rejected)
##
##     num_tiling_hits = 3 (A, B, C passed; D rejected)
##     tiling_ref = whichever of A/B/C gives highest tiling_cov_pct
##
## -----------------------------------------------------------------------------
##
## GENOMIC_DIST
##   Brief: Number of rank positions between the two genes on the chromosome,
##          after excluding nested genes from the rank.
##
##   Detail: genomic_dist=1 means the genes are directly adjacent in the
##           annotation with no other non-nested gene between them.
##           genomic_dist=2 means one gene sits between them, and so on.
##           We search up to max_dist=4. Note this counts annotated gene
##           positions, not base pairs — two genes could be 100kb apart
##           in sequence but have genomic_dist=1 if no other gene lies between.
##
##   Example from data:
##     CCA3g000130000.1 end:   7462203
##     CCA3g000132000.1 start: 7462205   <- 2bp gap in sequence
##     genomic_dist = 1                  <- directly adjacent in rank
##     (gene 131 was nested inside 130 and excluded from rank)
##
## -----------------------------------------------------------------------------
##
## EVALUES (in merge table)
##   Brief: Comma-separated diamond evalues, one per gene in genomic order,
##          all vs the chain's tiling_ref reference protein.
##
##   Detail: These come from your diamond blastp search of the annotation's proteins
##           vs the reference proteome — NOT from the evalue strings embedded
##           in your fasta headers (those are from a different search).
##           A '?' means that gene had no diamond hit vs the chain tiling_ref
##           that passed the evalue/coverage filters. In a chain of 3+ genes
##           tiling_ref comes from the best pair, and a gene outside that pair
##           may tile against a different reference protein. Check
##           split_genes_detail.txt for that gene pair to see what it does
##           tile against.
##
## -----------------------------------------------------------------------------
##
## SKIPPED_GENES
##   Brief: Gene IDs of any annotated genes that sit between a junction pair
##          when genomic_dist > 1.
##
##   Detail: When two genes in a merge candidate are not directly adjacent
##           (genomic_dist > 1), there are annotated genes sitting inside
##           the merged locus that were not identified as merge candidates
##           themselves. These skipped genes may be:
##             1. Genuinely separate genes unrelated to the merge
##             2. Split gene fragments that failed your filters (e.g. scovhsp
##                > 85%, no shared tiling hit, or too diverged to blast well)
##             3. Small annotation artifacts (TE fragments, tiny ORFs)
##           Check skipped genes in your GFF and blast output to determine
##           which case applies. If a skipped gene turns out to be another
##           fragment, it should be added to the merge manually.
##
## -----------------------------------------------------------------------------
##
## CHAIN TRIMMING AND SPLITTING
##   Brief: Multi-gene chains are refined by removing weak terminal genes and
##          splitting at weak internal junctions before writing to the merge
##          table.
##
##   Asymmetry threshold: $asym_threshold (default 6)
##   A junction is "weak" if it has exactly 1 tiling hit AND the maximum
##   tiling hits across any junction in the chain is >= asym_threshold.
##
##   TRIMMING (iterative, applied first):
##     If the first junction is weak, drop gene 1 from the chain.
##     If the last junction is weak, drop the last gene from the chain.
##     Repeat until neither terminal junction qualifies for dropping.
##
##     Example: chain hits = 1|8|6  (3 junctions, 4 genes)
##       First junction (1) is weak, max=8, 8>=6 -> drop gene 1
##       Remaining: 8|6  -> neither end is 1 -> done
##       Result: 3-gene chain with junctions 8|6
##
##   SPLITTING (applied after trimming):
##     Scan internal junctions. If a junction has 1 tiling hit AND both its
##     left and right neighbor junctions have >= asym_threshold hits,
##     split the chain at that junction into two separate merge rows.
##
##     Example: chain hits = 10|1|17  (3 junctions, 4 genes A,B,C,D)
##       Internal junction B-C has 1 hit, left=10, right=17, both >=6
##       -> split into A,B (junction=10) and C,D (junction=17)
##
##   Note: trimming and splitting only apply to chains with 3+ genes (2+ pairs).
##         2-gene pairs are never trimmed or split.
##
##   --no_asym_trim flag:
##     Disables all trimming and splitting. Chains are returned whole.
##     Weak terminal junctions are still flagged as WEAK_END. Use this when:
##       - Your reference proteome is distant (fewer tiling hits overall) and
##         the asym_threshold of 6 is rarely reached, making trimming too
##         aggressive.
##       - You want to inspect which genes would be trimmed before deciding
##         whether to discard them.
##     Compare flag distributions (WEAK_END count, chain lengths) between
##     runs with and without this flag to gauge its impact on your dataset.
##
## -----------------------------------------------------------------------------
##
## FLAG DEFINITIONS
##   One or more flags per merge row, comma-separated if multiple apply.
##   Flags document properties of the merge candidate — they are not all
##   equivalent in risk. Read interpretations below before deciding which
##   to add to skip_flags in the merge step.
##
##   Flags that apply to all chain sizes (2-gene pairs and longer):
##
##     LOW_COV            : tiling_cov_pct < 60%. The genes together cover
##                          less than 60% of the reference protein, suggesting
##                          domain sharing rather than a genuine split gene.
##                          Recommended to skip_flags in production runs.
##
##     PIECES_OVERLAP     : two pieces of the merge align to the SAME region of
##                          some reference protein: the overlap is >= overlap_min_frac
##                          (default 0.5) of the shorter of the two alignments, and
##                          each piece aligns over >= overlap_min_qcov (default 60%)
##                          of ITSELF there (E <= max_evalue). Pieces of one split
##                          gene are consecutive parts of a protein and sit side by
##                          side on every reference; two pieces on top of each
##                          other are copies of one gene unit (tandem duplicates),
##                          joined only because some other reference is a fused
##                          two-unit protein they tile on.
##                          Every reference hit by two pieces is checked, including
##                          the hits the fragment filter drops (a piece covering
##                          > 85% of a reference is a complete copy of it) -- the
##                          tiling test rejects an overlap on one reference but
##                          never looks at the others; this flag does.
##                          The qcov condition keeps shared domains out (two genes
##                          aligning 20% of themselves to one kinase domain are
##                          not copies). Risk: a real gene with internal repeats,
##                          cut so that each piece carries some repeats, can also
##                          overlap on a short single-repeat reference.
##                          Example (Congeria merge_468): pieces of 617 and 609 aa
##                          tile on a fused 1225 aa Dreissena protein (gap 15) but
##                          each covers ~99% of two single-unit Dreissena proteins,
##                          on the same residues.
##                          In skip_flags by default.
##
##     JUNCTION_OVERHANG  : at some junction a piece extends more than end_tol aa
##                          (default 30) past its alignment on that junction's
##                          tiling reference, on the junction side: sequence the
##                          reference lacks at the break (a mispredicted exon at
##                          the break, or two different proteins).
##
##     SKIPPED_GENE       : a non-adjacent gene sits inside the merge locus
##                          (genomic_dist > 1 for at least one junction) AND at
##                          least one skipped gene is on the SAME strand as the
##                          chain. The skipped gene may be an additional split
##                          fragment that failed filters, a different gene, or
##                          an artifact. Check skipped_genes column and GFF
##                          before merging. Recommended to skip_flags and review
##                          manually.
##
##     OPPOSITE_STRAND_SKIP : a non-adjacent gene sits inside the merge locus,
##                          but ALL skipped genes are on the OPPOSITE strand.
##                          These are unrelated interleaved genes and do not
##                          affect the validity of the merge. Not included in
##                          skip_flags by default.
##
##     MULTI_ISOFORM_JOIN : at least one source gene has more than one annotated
##                          transcript. The merged gene will be built by cross-
##                          product (N*M transcript combinations). Not all
##                          cross-product transcripts are necessarily real —
##                          they require transcript-level evidence to validate.
##                          Review isoform structure in a genome browser before
##                          relying on isoform-level annotations from this merge.
##                          The miniprot rebuild (step 6) keeps one transcript and
##                          records the unused source isoforms (isoforms_not_used=).
##
##     LARGE_SPAN         : the total genomic footprint of the merged locus
##                          exceeds large_span_warn (default 500 kb). Biologically
##                          plausible for some large gene families (e.g., neurexins,
##                          ROBO receptors, dystrophin) but cross-check with IsoSeq
##                          spanning evidence before acting on weak-evidence merges.
##
##     LARGE_SPAN_EXTREME : the total genomic footprint exceeds large_span_extreme
##                          (default 2 Mb). Only a handful of vertebrate genes span
##                          this range. Strong FP risk unless FULL_SPAN IsoSeq
##                          confirms. Recommended to add to skip_flags unless you
##                          have independent evidence.
##
##   For 2-gene pairs (single junction, no trimming/splitting):
##
##     SINGLE_HIT         : the junction has exactly 1 tiling hit — only one
##                          reference protein supports this merge.
##
##                          IMPORTANT — this is a REVIEW flag, not a default
##                          skip flag. Interpretation depends on the reference
##                          proteome and the gene family:
##
##                          When the reference proteome contains only one copy of
##                          a gene (single-copy gene), that gene has at most one
##                          protein to tile against. For such genes, 1 tiling hit
##                          is the expected maximum — it reflects the proteome's
##                          gene copy number, not evidence weakness. A SINGLE_HIT
##                          merge of a clearly described single-copy gene with a
##                          strong evalue and good combined coverage is often as
##                          reliable as a multi-hit merge from a gene family.
##
##                          Genuine risk arises when the single hit is against a
##                          protein from a large multi-domain family (EGF repeats,
##                          Ig domains, WD40, kinase domains, etc.) where domain
##                          sharing among unrelated proteins could produce a
##                          coincidental tiling alignment. In these cases, check
##                          pident: low identity (<30%) on a single hit is a red
##                          flag for domain-sharing false positives.
##
##                          Recommended action: inspect tiling_ref_desc and pident. If
##                          the hit is a well-described single-copy gene with
##                          strong evalue, proceed. If the hit is a generic domain
##                          protein with low pident, require IsoSeq confirmation.
##
##     (no flag)          : 2+ tiling hits, coverage >= 60%, directly adjacent —
##                          solid candidate
##
##   For 3+ gene chains (after trimming and splitting):
##
##     STRONG             : all junctions have >= 3 tiling hits — high confidence.
##
##     TRANSITIVE_JOIN    : one or more consecutive gene pairs in this chain have
##                          NO direct pairwise tiling evidence between them. Their
##                          junction is "?" in junction_tiling_hits. These genes
##                          were placed in the same chain because each is
##                          independently connected to another gene in the chain
##                          through a different reference protein — not because
##                          they themselves share a common reference where their
##                          alignments tile end-to-end.
##
##                          In graph terms: the chain is connected, but the edge
##                          between at least one consecutive gene pair is absent
##                          and the connection is inferred transitively through
##                          shared membership in overlapping pairs.
##
##                          Two genes with no shared tiling reference may be:
##                            (a) Genuine split fragments whose split site falls
##                                in a region too diverged for shared blast hits.
##                                IsoSeq spanning the ? junction confirms the merge.
##                            (b) Two different genes (paralogs, different family
##                                members, syntenic cluster neighbors) that each
##                                tile with a shared third gene but not with each
##                                other. This is the mechanism behind the ROBO
##                                case (same paralog type, different strands joined
##                                via shared N-terminal alignment) and the GABRB2
##                                case (gamma and beta GABA receptor subunits joined
##                                via shared Cys-loop fold on a third reference).
##
##                          Key diagnostic: if the gene flanking the ? junction has
##                          "?" in the evalues column vs the chain tiling_ref, that
##                          gene has no direct evidence against the chain's own
##                          reference protein — a strong signal it does not belong.
##                          Also check that descriptions on both sides of the ?
##                          match. FULL_SPAN IsoSeq crossing the ? junction is the
##                          most reliable confirmation.
##
##                          TRANSITIVE_JOIN compounds other flags.
##                          STRONG,TRANSITIVE_JOIN warrants the same caution as a
##                          WEAK_END case. SINGLE_HIT,TRANSITIVE_JOIN with no
##                          IsoSeq is the highest-risk category in the merge table.
##
##     WEAK_END           : a terminal junction has 1-2 hits but survived trimming
##                          (asymmetry threshold not met). Low evidence at one end.
##
##     WEAK_INTERNAL      : an internal junction has 1-2 hits but both neighbors
##                          did not both meet the threshold for splitting.
##
##     SINGLE_HIT         : all junctions have exactly 1 tiling hit.
##                          Same interpretation caveats as for 2-gene pairs above.
##
##     LOW_COV            : tiling_cov_pct < 60%
##
##     PIECES_OVERLAP     : two pieces (any two, not only neighbours) align to the
##                          same region of some reference. See the definition above.
##
##     SKIPPED_GENE       : at least one junction spans a non-adjacent same-strand
##                          gene. See SKIPPED_GENE description in the 2-gene section.
##
##     OPPOSITE_STRAND_SKIP : all skipped genes at non-adjacent junctions are on
##                          the opposite strand; does not affect merge validity.
##
##   Interpreting combinations:
##     STRONG                         -> act on this merge
##     STRONG,SKIPPED_GENE            -> strong evidence but check skipped gene;
##                                       it may need to be added to the merge
##     STRONG,OPPOSITE_STRAND_SKIP    -> skipped gene is on opposite strand;
##                                       safe to merge (spanning_rescue can help)
##     STRONG,TRANSITIVE_JOIN         -> inspect the ? junction carefully;
##                                       require IsoSeq if possible
##     SINGLE_HIT                     -> inspect tiling_ref_desc and pident; see above
##     SINGLE_HIT,TRANSITIVE_JOIN     -> high risk; require IsoSeq FULL_SPAN
##     WEAK_END,LOW_COV               -> low confidence; review carefully
##     STRONG,LOW_COV                 -> protein evidence consistent but incomplete;
##                                       may be domain sharing, not a split gene
##     LARGE_SPAN,SINGLE_HIT          -> high FP risk; require IsoSeq
##     LARGE_SPAN_EXTREME             -> add to skip_flags unless IsoSeq confirms
##     MULTI_ISOFORM_JOIN             -> merge proceeds; verify transcript structure
##     PIECES_OVERLAP (any others)    -> likely tandem copies, not a split gene;
##                                       skipped by default, review by eye
##
## =============================================================================

my $no_asym_trim       = 0;
my $large_span_warn    = 500000;   # bp: flag LARGE_SPAN above this
my $large_span_extreme = 2000000;  # bp: flag LARGE_SPAN_EXTREME above this
my $wiggle         = 15;  # aa wiggle room for tiling junction
my $max_dist       =  4;  # max genomic rank distance between gene pair
my $asym_threshold =  6;  # min hits on strong side to trigger trim/split
my $low_cov_thresh = 60;  # tiling_cov_pct below this gets LOW_COV flag
my $max_evalue     = 1e-20;  # a fragment hit must have evalue <= this
my $frag_min_qcov  = 30;  # a fragment must align over >= this % of ITSELF (qcovhsp)
my $ref_min_qcov   = 60;  # rebuild_candidates: references all chain genes align over >= this % of themselves first
my $end_tol        = 30;  # aa a fragment may extend past its alignment at the junction side
my $overlap_min_frac = 0.5;  # PIECES_OVERLAP: overlap >= this fraction of the shorter alignment
my $overlap_min_qcov = 60;   # PIECES_OVERLAP: each piece aligns over >= this % of ITSELF there

GetOptions(
    "no_asym_trim"         => \$no_asym_trim,
    "large_span_warn=i"    => \$large_span_warn,
    "large_span_extreme=i" => \$large_span_extreme,
    "max_dist=i"           => \$max_dist,
    "wiggle=i"             => \$wiggle,
    "asym_threshold=i"     => \$asym_threshold,
    "low_cov_thresh=f"     => \$low_cov_thresh,
    "max_evalue=f"         => \$max_evalue,
    "frag_min_qcov=f"      => \$frag_min_qcov,
    "ref_min_qcov=f"       => \$ref_min_qcov,
    "end_tol=i"            => \$end_tol,
    "overlap_min_frac=f"   => \$overlap_min_frac,
    "overlap_min_qcov=f"   => \$overlap_min_qcov,
) or die usage_text();

sub usage_text {
    return "usage: $0 [--no_asym_trim] [--large_span_warn N] [--large_span_extreme N] " .
           "[--max_dist N] [--wiggle N] [--asym_threshold N] [--low_cov_thresh N] " .
           "[--max_evalue E] [--frag_min_qcov N] [--ref_min_qcov N] [--end_tol N] " .
           "[--overlap_min_frac F] [--overlap_min_qcov N] blast gff subject_fa query_fa\n";
}
my $usage_str = usage_text();
my $blast      = shift or die $usage_str;
my $gff        = shift or die $usage_str;
my $subject_fa = shift or die $usage_str;
my $query_fa   = shift or die $usage_str;
my $strong_thresh  =  3;  # min junction hits for STRONG flag
my $frag_scov_max  = 85;  # scovhsp <= this is always a fragment candidate
my $frag_qcov_min  = 90;  # qcovhsp >= this enables the soft scovhsp extension
my $frag_scov_soft = 90;  # with qcovhsp >= frag_qcov_min, keep if scovhsp < this
my $max_ref_candidates = 10;  # ranked reference proteins written per merge (rebuild_candidates)

# ---------------------------------------------------------------------------
# Load reference protein descriptions from subject fasta headers
# Parses gene_symbol and description fields from Ensembl-style headers
# ---------------------------------------------------------------------------
my %desc;
open FA, $subject_fa or die "cant open subject fasta: $subject_fa $!\n";
while (my $line = <FA>){
    chomp $line;
    next unless $line =~ /^>/;
    my ($id)     = $line =~ /^>(\S+)/;
    my ($symbol) = $line =~ /gene_symbol:(\S+)/;
    my ($d)      = $line =~ /description:(.+?)(?:\s+\[|$)/;
    my ($gene)   = $line =~ /\bgene:(ENS\S+)/;
    my @parts;
    for my $part ($symbol, $d) { push @parts, $part if defined $part }
    my $joined   = join(" | ", @parts);
    $desc{$id}   = $joined ne "" ? $joined : ($gene // "");
}
close FA;

# ---------------------------------------------------------------------------
# Build transcript->gene mapping directly from GFF mRNA features.
# Also count transcripts per gene — used for MULTI_ISOFORM_JOIN flag.
# ---------------------------------------------------------------------------
my %isoforms;           # transcript_id -> gene_id
my %gene_n_transcripts; # gene_id -> number of annotated transcripts
open GFF0, $gff or die "cant open gff: $gff $!\n";
while (my $line = <GFF0>) {
    chomp $line;
    next if $line =~ /^#/;
    next unless $line =~ /\tmRNA\t/;
    my @f = split "\t", $line;
    my ($tid)     = $f[8] =~ /ID=([^;]+)/;
    my ($gene_id) = $f[8] =~ /Parent=([^;]+)/;
    next unless defined $tid && defined $gene_id;
    $isoforms{$tid} = $gene_id;
    $gene_n_transcripts{$gene_id}++;
}
close GFF0;

# Load query fasta descriptions into %desc, keyed by both transcript and gene ID
# Must run after %isoforms is populated for transcript->gene resolution
open FA, $query_fa or die "cant open query fasta: $query_fa $!\n";
while (my $line = <FA>){
    chomp $line;
    next unless $line =~ /^>/;
    my ($id) = $line =~ /^>(\S+)/;
    my ($d)  = $line =~ /^>\S+\s+(.+)/;
    $desc{$id} = $d // "";
    my $gene_id = $isoforms{$id} // $id;
    $desc{$gene_id} //= $d // "";
}
close FA;

# ---------------------------------------------------------------------------
# Load gene coordinates from GFF, gene features only
# ---------------------------------------------------------------------------
my %locs;  # gene_id -> {ref, start, end, strand, coord}
open GFF, $gff or die "cant open gff: $gff $!\n";
while (my $line = <GFF>){
    chomp $line;
    next if $line =~ /^#/;
    next if $line !~ /\tgene\t/;
    my @f = split "\t", $line;
    my ($gene_id) = $f[8] =~ /ID=([^;]+)/;
    next unless defined $gene_id;
    $locs{$gene_id}{ref}    = $f[0];
    $locs{$gene_id}{start}  = $f[3];
    $locs{$gene_id}{end}    = $f[4];
    $locs{$gene_id}{strand} = $f[6];
    $locs{$gene_id}{coord}  = "$f[0].$f[3].$f[4]";
}
close GFF;

# ---------------------------------------------------------------------------
# Comparators for sort, named so the code never reads sort's package variables: with the
# ($$) prototype, sort passes the two elements in @_
# ---------------------------------------------------------------------------
# genes by chromosome, start, the longer gene first (it is the parent), id
sub by_genome_position($$) {
    my ($left, $right) = @_;
    return $locs{$left}{ref}   cmp $locs{$right}{ref}
        || $locs{$left}{start} <=> $locs{$right}{start}
        || $locs{$right}{end}  <=> $locs{$left}{end}
        || $left cmp $right;
}
# genes of one chromosome by start, then id
sub by_start_then_id($$) {
    my ($left, $right) = @_;
    return $locs{$left}{start} <=> $locs{$right}{start} || $left cmp $right;
}
sub by_number_ascending($$)  { my ($left, $right) = @_; return $left <=> $right }
sub by_number_descending($$) { my ($left, $right) = @_; return $right <=> $left }
sub by_pct_descending($$)    { my ($left, $right) = @_; return $right->{pct} <=> $left->{pct} }
sub by_tiling_pct_descending($$) { my ($left, $right) = @_; return $right->{tiling_pct} <=> $left->{tiling_pct} }
sub by_num_tiling_descending($$) { my ($left, $right) = @_; return $right->{num_tiling} <=> $left->{num_tiling} }

# the first element of a sorted copy (sort in list context, then take [0])
sub first_sorted_desc { my @sorted = sort by_number_descending @_; return $sorted[0] }
sub first_sorted_asc  { my @sorted = sort by_number_ascending @_;  return $sorted[0] }

# Sort all genes by chromosome then start position
my @by_start = sort by_genome_position keys %locs;

# ---------------------------------------------------------------------------
# Identify nested genes: genes whose coordinates fall entirely within another
# gene. These are excluded from the rank array because they would inflate the
# apparent genomic distance between flanking genes.
# Genes are visited in start order keeping, per chromosome, the gene reaching
# furthest so far: a gene ending at or before that one's end lies inside it. (A
# fixed look-back window missed genes inside long genes with more than 10 genes
# in their introns.)
# ---------------------------------------------------------------------------
my %nested;  # nested_gene -> parent_gene
my ($reach_ref, $reach_gene) = ("", undef);
for my $gene (@by_start) {
    if ($locs{$gene}{ref} ne $reach_ref) {
        ($reach_ref, $reach_gene) = ($locs{$gene}{ref}, $gene);
        next;
    }
    if ($locs{$gene}{end} <= $locs{$reach_gene}{end}) {
        $nested{$gene} = $reach_gene;
    } else {
        $reach_gene = $gene;
    }
}

# Build rank array and rank lookup excluding nested genes
my @all_sorted;
for my $gene (@by_start) { push @all_sorted, $gene unless exists $nested{$gene} }
my %rank;
for my $i (0 .. $#all_sorted) {
    $rank{ $all_sorted[$i] } = $i;
}

print STDERR "Total genes: ",           scalar keys %locs,   "\n";
print STDERR "Nested genes excluded: ", scalar keys %nested, "\n";
print STDERR "Genes in rank: ",         scalar @all_sorted,  "\n";

# ---------------------------------------------------------------------------
# Parse diamond output
# Only store hits where evalue <= 1e-20 AND scovhsp <= 85%
# scovhsp <= 85% is the fragment filter — genes with higher coverage are
# complete homologs, not split gene fragments
#
# Two data structures:
#   subjects{hit}{gene} = {cov, sstart, send, evalue}
#     stores the best-coverage hit per gene per subject protein
#   gene_hits{gene}{hit} = 1
#     inverted index for fast intersection when finding shared hits
#   copy_hsp{gene}{hit} = {bits, sstart, send}
#     every hit (fragment filter or not) where the gene aligns over
#     >= overlap_min_qcov of itself, best bitscore per gene and subject;
#     used only by the PIECES_OVERLAP check
# ---------------------------------------------------------------------------
my %subjects;
my %gene_hits;
my %copy_hsp;
my %subject_len;  # subject_id -> length in aa (from slen field)

open BLAST, $blast or die "cant open blast output: $blast $!\n";
while (my $line = <BLAST>){
    chomp $line;
    next if $line =~ /^#/;
    my ($query, $subject, $per_identity, $alignment_length, $mismatches,
        $gap_opens, $qstart, $qend, $sstart, $send, $evalue, $bitscore,
        $scovhsp, $slen, $qcovhsp, $qlen) = split "\t", $line;

    # resolve transcript to gene level
    my $gene_id = $isoforms{$query} // $query;
    $subject_len{$subject} = $slen;

    # PIECES_OVERLAP evidence, taken before the fragment filter: a piece covering
    # > 85% of a reference is a complete copy of it, which is what this check needs
    if ($evalue <= $max_evalue && defined $qcovhsp && $qcovhsp >= $overlap_min_qcov
        && (!exists $copy_hsp{$gene_id}{$subject} || $bitscore > $copy_hsp{$gene_id}{$subject}{bits})) {
        $copy_hsp{$gene_id}{$subject} = { bits => $bitscore, sstart => $sstart, send => $send };
    }

    # ---------------------------------------------------------------------------
    # FRAGMENT FILTER
    # Keep hits where the aligned gene appears to be a partial homolog, not a
    # complete ortholog. Two conditions define a fragment:
    #
    #   Classic:  scovhsp <= frag_scov_max (default 85%)
    #     The alignment covers at most 85% of the reference protein. The gene
    #     is too short to represent the full-length ortholog.
    #
    #   Soft extension (requires qcovhsp in diamond output):
    #     qcovhsp >= frag_qcov_min (default 90%)
    #     AND scovhsp < frag_scov_soft (default 90%)
    #     The alignment covers nearly the entire query protein (the annotated
    #     CDS is internally complete) but accounts for less than 90% of the
    #     reference protein. This catches split fragments of shorter proteins
    #     — e.g., a 200aa annotated CDS that is genuinely a fragment of a 280aa
    #     gene — which fall in the 86-89% scovhsp range and would be missed by
    #     the classic filter alone.
    #
    # The classic and soft conditions are OR'd: either is sufficient.
    # If qcovhsp is absent from the blast output (undef), only the classic
    # condition applies and behavior is identical to pre-qcovhsp versions.
    # ---------------------------------------------------------------------------
    my $is_fragment = ($scovhsp <= $frag_scov_max)
        || (defined $qcovhsp && $qcovhsp >= $frag_qcov_min
            && $scovhsp < $frag_scov_soft);
    # A split fragment is a piece of the reference protein: it aligns over most of
    # ITSELF. A long multidomain protein aligning over one domain of the reference
    # also has scovhsp <= 85 but is domain sharing, not a fragment.
    $is_fragment = 0 if defined $qcovhsp && $qcovhsp < $frag_min_qcov;

    if ($evalue <= $max_evalue && $is_fragment) {
        # keep only the best-coverage HSP per gene-subject pair (ties: lower evalue)
        if (!exists $subjects{$subject}{$gene_id} ||
             $subjects{$subject}{$gene_id}{cov} < $scovhsp ||
            ($subjects{$subject}{$gene_id}{cov} == $scovhsp && $evalue < $subjects{$subject}{$gene_id}{evalue})){
            $subjects{$subject}{$gene_id}{cov}    = $scovhsp;
            $subjects{$subject}{$gene_id}{sstart} = $sstart;
            $subjects{$subject}{$gene_id}{send}   = $send;
            $subjects{$subject}{$gene_id}{evalue} = $evalue;
            $subjects{$subject}{$gene_id}{qstart} = $qstart;
            $subjects{$subject}{$gene_id}{qend}   = $qend;
            $subjects{$subject}{$gene_id}{qlen}   = $qlen;
            $subjects{$subject}{$gene_id}{pident} = $per_identity;
            $subjects{$subject}{$gene_id}{qcov}   = $qcovhsp;
        }
        $gene_hits{$gene_id}{$subject} = 1;
    }
}
close BLAST;

print STDERR "Subject proteins with fragment hits: ", scalar keys %subjects,  "\n";
print STDERR "Genes with fragment hits: ",            scalar keys %gene_hits, "\n";

# ---------------------------------------------------------------------------
# Output files
# ---------------------------------------------------------------------------
open SUMMARY, ">split_genes_summary.txt" or die "cant open summary output $!\n";
open DETAIL,  ">split_genes_detail.txt"  or die "cant open detail output $!\n";

print SUMMARY join("\t",
    "g1", "g2",
    "g1_coord", "g2_coord",
    "genomic_dist",
    "g1_desc", "g2_desc",
    "tiling_ref", "tiling_ref_desc",
    "tiling_cov_pct",
    "tiling_gap",
    "g1_subj", "g2_subj",
    "g1_cov", "g2_cov",
    "g1_evalue", "g2_evalue",
    "num_tiling_hits"
), "\n";

print DETAIL join("\t",
    "g1", "g2",
    "genomic_dist",
    "tiling_ref", "tiling_ref_desc",
    "tiling_cov_pct",
    "tiling_gap",
    "g1_subj", "g2_subj",
    "g1_cov", "g2_cov",
    "g1_evalue", "g2_evalue"
), "\n";

# Collect all reported pairs for downstream chaining into merge candidates
# Also build a lookup: "g1\tg2" -> {num_tiling_hits, genomic_dist}
my @all_pairs;
my $n_not_collinear = 0;   # tiling tests skipped: the two pieces in the wrong order for the strand
my %pair_overhang;  # "g1\tg2" -> aa past the alignment at the junction (best hit)
my %pair_tiling;  # "g1\tg2" -> num_tiling_hits
my %pair_dist;    # "g1\tg2" -> genomic_dist

# ---------------------------------------------------------------------------
# Main loop: for each ranked gene, check the next max_dist neighbors
# For each neighbor pair, find shared reference protein hits and test tiling
# ---------------------------------------------------------------------------
for my $i (0 .. $#all_sorted - 1) {
    my $g1 = $all_sorted[$i];
    next unless exists $gene_hits{$g1};  # skip genes with no fragment hits

    my $limit = ($i + $max_dist) > $#all_sorted ? $#all_sorted : $i + $max_dist;

    for my $j ($i+1 .. $limit) {
        my $g2 = $all_sorted[$j];
        next unless $locs{$g1}{ref} eq $locs{$g2}{ref};
        next unless $locs{$g1}{strand} eq $locs{$g2}{strand};  # same-strand only
        next unless exists $gene_hits{$g2};

        # Find shared reference hits by intersecting inverted indexes.
        # Iterating the smaller set is faster than scanning all subjects.
        my ($small, $big) =
            scalar keys %{$gene_hits{$g1}} < scalar keys %{$gene_hits{$g2}}
            ? ($gene_hits{$g1}, $gene_hits{$g2})
            : ($gene_hits{$g2}, $gene_hits{$g1});

        my @shared;
        for my $hit (keys %$small) { push @shared, $hit if exists $big->{$hit} }
        @shared = sort @shared;   # sorted: ties never depend on hash order
        next unless @shared;

        my $tiling_ref    = undef;
        my $tiling_pct    = 0;
        my $tiling_gap    = undef;
        my $best_g1_range = undef;
        my $best_g2_range = undef;
        my $best_g1_cov   = undef;
        my $best_g2_cov   = undef;
        my $best_g1_eval  = undef;
        my $best_g2_eval  = undef;
        my $num_tiling    = 0;
        my $best_ev_sum   = undef;
        my $best_overhang = 0;
        my @tiling_hits;  # all hits that pass tiling filter, for detail output

        for my $hit (@shared) {
            my $s1 = $subjects{$hit}{$g1}{sstart};
            my $e1 = $subjects{$hit}{$g1}{send};
            my $s2 = $subjects{$hit}{$g2}{sstart};
            my $e2 = $subjects{$hit}{$g2}{send};

            # Tiling gap: how far apart are the two alignments on the subject?
            # Positive = gap between them, negative = they overlap slightly.
            # We allow abs(gap) <= wiggle to account for alignment fuzziness.
            # Collinearity: g1 is upstream of g2 in the genome. On + the upstream
            # gene encodes the N-terminal part, so g1 must align before g2 on the
            # reference; on - the reverse. Pieces in the wrong order are not one
            # gene split in two (inverted tandem copies, domain shuffles).
            my $g1_first = $s1 < $s2 ? 1 : 0;
            if ($g1_first != ($locs{$g1}{strand} eq "-" ? 0 : 1)) {
                $n_not_collinear++;
                next;
            }
            my $gap = ($s1 < $s2) ? $s2 - $e1 : $s1 - $e2;
            next if abs($gap) > $wiggle;

            $num_tiling++;

            # Junction overhang: the N-terminal piece should align up to its own
            # C-terminus, the C-terminal piece from its own N-terminus; residues past
            # the alignment at the junction are sequence the reference lacks there
            my ($n_gene, $c_gene) = $g1_first ? ($g1, $g2) : ($g2, $g1);
            my $n_hit = $subjects{$hit}{$n_gene};
            my $c_hit = $subjects{$hit}{$c_gene};
            my $overhang = 0;
            $overhang = $n_hit->{qlen} - $n_hit->{qend} if defined $n_hit->{qlen} && $n_hit->{qlen} ne "";
            my $c_over = defined $c_hit->{qstart} ? $c_hit->{qstart} - 1 : 0;
            $overhang = $c_over if $c_over > $overhang;

            # Combined span: from the leftmost sstart to the rightmost send
            my $min_s = $s1 < $s2 ? $s1 : $s2;
            my $max_e = $e1 > $e2 ? $e1 : $e2;
            my $span  = $max_e - $min_s + 1;
            my $pct   = $subject_len{$hit}
                      ? 100 * $span / $subject_len{$hit}
                      : 0;

            push @tiling_hits, {
                hit    => $hit,
                pct    => $pct,
                overhang => $overhang,
                gap    => $gap,
                g1_r   => "$s1-$e1",
                g2_r   => "$s2-$e2",
                g1_cov => $subjects{$hit}{$g1}{cov},
                g2_cov => $subjects{$hit}{$g2}{cov},
                g1_ev  => $subjects{$hit}{$g1}{evalue},
                g2_ev  => $subjects{$hit}{$g2}{evalue},
            };

            # best = highest combined coverage; ties by summed evalue, then subject id
            my $ev_sum = $subjects{$hit}{$g1}{evalue} + $subjects{$hit}{$g2}{evalue};
            if (!defined $tiling_ref || $pct > $tiling_pct
                || ($pct == $tiling_pct && $ev_sum < $best_ev_sum)) {
                $best_ev_sum   = $ev_sum;
                $best_overhang = $overhang;
                $tiling_pct    = $pct;
                $tiling_ref    = $hit;
                $tiling_gap    = $gap;
                $best_g1_range = "$s1-$e1";
                $best_g2_range = "$s2-$e2";
                $best_g1_cov   = $subjects{$hit}{$g1}{cov};
                $best_g2_cov   = $subjects{$hit}{$g2}{cov};
                $best_g1_eval  = $subjects{$hit}{$g1}{evalue};
                $best_g2_eval  = $subjects{$hit}{$g2}{evalue};
            }
        }

        next unless defined $tiling_ref;

        my $gdist    = $j - $i;
        my $g1_desc  = $desc{$g1} // "?";
        my $g2_desc  = $desc{$g2} // "?";
        my $tiling_ref_desc = $desc{$tiling_ref} // "?";

        print SUMMARY join("\t",
            $g1, $g2,
            $locs{$g1}{coord}, $locs{$g2}{coord},
            $gdist,
            $g1_desc, $g2_desc,
            $tiling_ref, $tiling_ref_desc,
            sprintf("%.1f", $tiling_pct),
            $tiling_gap,
            $best_g1_range, $best_g2_range,
            $best_g1_cov, $best_g2_cov,
            $best_g1_eval, $best_g2_eval,
            $num_tiling
        ), "\n";

        # Detail: all tiling hits sorted by combined coverage descending
        for my $th (sort by_pct_descending @tiling_hits) {
            print DETAIL join("\t",
                $g1, $g2,
                $gdist,
                $th->{hit}, $desc{$th->{hit}} // "?",
                sprintf("%.1f", $th->{pct}),
                $th->{gap},
                $th->{g1_r}, $th->{g2_r},
                $th->{g1_cov}, $th->{g2_cov},
                $th->{g1_ev}, $th->{g2_ev}
            ), "\n";
        }

        # store pair for chaining and junction lookups
        push @all_pairs, {
            g1         => $g1,
            g2         => $g2,
            gdist      => $gdist,
            tiling_ref => $tiling_ref,
            tiling_pct => $tiling_pct,
            tiling_gap => $tiling_gap,
            g1_cov     => $best_g1_cov,
            g2_cov     => $best_g2_cov,
            g1_eval    => $best_g1_eval,
            g2_eval    => $best_g2_eval,
            num_tiling => $num_tiling,
            overhang   => $best_overhang,
        };
        $pair_tiling{"$g1\t$g2"} = $num_tiling;
        $pair_overhang{"$g1\t$g2"} = $best_overhang;
        $pair_dist{"$g1\t$g2"}   = $gdist;
    }
}

close SUMMARY;
close DETAIL;

# ---------------------------------------------------------------------------
# Chain pairs into multi-gene merge candidates using union-find.
# Any two pairs that share a gene are merged into one group.
# Example: pair A-B and pair B-C both contain B, so they become group A,B,C.
# ---------------------------------------------------------------------------
my @chains;
for my $pair (@all_pairs) {
    my ($g1, $g2) = ($pair->{g1}, $pair->{g2});

    my @matching_idx;
    for my $idx (0 .. $#chains) {
        push @matching_idx, $idx if $chains[$idx]{genes}{$g1} || $chains[$idx]{genes}{$g2};
    }

    if (!@matching_idx) {
        push @chains, {
            genes => { $g1 => 1, $g2 => 1 },
            pairs => [ $pair ]
        };
    } else {
        my %merged_genes = ($g1 => 1, $g2 => 1);
        my @merged_pairs = ($pair);
        for my $idx (@matching_idx) {
            %merged_genes = (%merged_genes, %{ $chains[$idx]{genes} });
            push @merged_pairs, @{ $chains[$idx]{pairs} };
        }
        for my $idx (reverse @matching_idx) {
            splice @chains, $idx, 1;
        }
        push @chains, {
            genes => \%merged_genes,
            pairs => \@merged_pairs
        };
    }
}

# ---------------------------------------------------------------------------
# Helper: given an ordered gene list, get junction tiling hits array
# Returns array of hit counts (undef if pair not found)
# ---------------------------------------------------------------------------
sub get_junction_hits {
    my @genes = @_;
    my @hits;
    for my $i (0 .. $#genes - 1) {
        push @hits, $pair_tiling{"$genes[$i]\t$genes[$i+1]"};
    }
    return @hits;
}

# ---------------------------------------------------------------------------
# Helper: given an ordered gene list, get junction genomic distances
# Returns array of distances (undef if pair not found)
# ---------------------------------------------------------------------------
# the values with undef shown as "?"
sub unknown_as_question {
    my @shown;
    for my $value (@_) { push @shown, $value // "?" }
    return @shown;
}

sub get_junction_dists {
    my @genes = @_;
    my @dists;
    for my $i (0 .. $#genes - 1) {
        push @dists, $pair_dist{"$genes[$i]\t$genes[$i+1]"};
    }
    return @dists;
}

# ---------------------------------------------------------------------------
# Helper: given an ordered gene list, get any skipped gene IDs
# A skipped gene is one that sits between a junction pair (dist > 1)
# Returns list of skipped gene IDs in genomic order
# ---------------------------------------------------------------------------
sub get_skipped_genes {
    my @genes = @_;
    my @skipped;
    for my $i (0 .. $#genes - 1) {
        my $g1   = $genes[$i];
        my $g2   = $genes[$i+1];
        my $dist = $pair_dist{"$g1\t$g2"} // 1;
        next unless $dist > 1;
        # genes between g1 and g2 in the rank array
        my $r1 = $rank{$g1};
        my $r2 = $rank{$g2};
        for my $r ($r1+1 .. $r2-1) {
            push @skipped, $all_sorted[$r];
        }
    }
    return @skipped;
}

# ---------------------------------------------------------------------------
# Helper: apply trimming and splitting to a gene list, return list of
# refined gene lists (may be multiple if chain was split)
#
# Trimming: iteratively drop terminal genes whose junction has 1 tiling hit
#   AND the max junction hits in the chain >= asym_threshold
# Splitting: after trimming, split at internal junctions with 1 tiling hit
#   where both neighboring junctions have >= asym_threshold hits
# ---------------------------------------------------------------------------
sub refine_chain {
    my @genes = @_;

    # 2-gene pairs are never trimmed or split
    return ([@genes]) if @genes <= 2;

    # When --no_asym_trim is active, skip all trimming and splitting.
    # Chains are returned whole; WEAK_END junctions will be flagged but not removed.
    return ([@genes]) if $no_asym_trim;

    # --- TRIMMING PASS (iterative) ---
    my $trimmed = 1;
    while ($trimmed && @genes > 2) {
        $trimmed = 0;
        my @hits = get_junction_hits(@genes);
        my @known_hits;
        for my $hit_count (@hits) { push @known_hits, $hit_count if defined $hit_count }
        my $max_hits = first_sorted_desc(@known_hits) // 0;

        # check first junction
        if (defined $hits[0] && $hits[0] == 1 && $max_hits >= $asym_threshold) {
            shift @genes;
            $trimmed = 1;
            next;  # restart loop after trim
        }

        # check last junction
        if (defined $hits[-1] && $hits[-1] == 1 && $max_hits >= $asym_threshold) {
            pop @genes;
            $trimmed = 1;
        }
    }

    # if trimmed down to 1 gene, discard
    return () if @genes < 2;

    # 2-gene result after trimming — return as-is, no splitting needed
    return ([@genes]) if @genes <= 2;

    # --- SPLITTING PASS ---
    # Scan internal junctions (not first, not last) for weak links
    # A junction splits if: hits==1 AND both neighbors >= asym_threshold
    my @hits = get_junction_hits(@genes);
    my @split_at;  # indices in @hits where we should split

    for my $i (1 .. $#hits - 1) {  # internal junctions only
        next unless defined $hits[$i] && $hits[$i] == 1;
        my $left  = $hits[$i-1] // 0;
        my $right = $hits[$i+1] // 0;
        if ($left >= $asym_threshold && $right >= $asym_threshold) {
            push @split_at, $i;
        }
    }

    # no splits needed
    return ([@genes]) unless @split_at;

    # split the gene list at the identified junctions
    my @pieces;
    my $start = 0;
    for my $split_i (@split_at) {
        push @pieces, [@genes[$start .. $split_i]];
        $start = $split_i + 1;
    }
    push @pieces, [@genes[$start .. $#genes]];

    # discard any pieces with fewer than 2 genes
    my @kept;
    for my $piece (@pieces) { push @kept, $piece if @$piece >= 2 }
    return @kept;
}

# ---------------------------------------------------------------------------
# Helper: rank the reference proteins of a chain. The ranked list
# (rebuild_candidates) is what rebuild_with_miniprot.pl tries, in order, to
# rebuild the merged gene.
#   1. most chain genes hit (the evidence that joins them)
#   1b. every one of those genes aligns over >= ref_min_qcov of ITSELF (a
#      reference that shares only a domain with a piece ranks after the ones the
#      pieces align to as a whole; matters when frag_min_qcov is lower)
#   2. identity-weighted coverage: the fraction of the reference's residues
#      aligned identically by the chain's fragments, sum(pident * aligned span)
#      / reference length. A close reference scores higher than a distant one
#      at the same coverage; a distant reference makes spurious frameshifts and
#      misplaced splice sites when aligned to the genome.
#   3. a described hit (readability), then the reference id (deterministic)
# Only references hitting two or more chain genes are candidates.
# ---------------------------------------------------------------------------
sub identity_weighted_cov {
    my ($hit, @genes) = @_;
    my $slen = $subject_len{$hit} or return 0;
    my $identical = 0;
    for my $gene (@genes) {
        my $hsp = $subjects{$hit}{$gene} or next;
        my $pident = $hsp->{pident};
        $pident = 0 unless defined $pident && $pident ne "";
        $identical += $pident / 100 * ($hsp->{send} - $hsp->{sstart} + 1);
    }
    my $frac = $identical / $slen;
    return $frac > 1 ? 1 : $frac;
}

# a named comparator: sort passes the two elements in @_ (the ($$) prototype)
sub by_reference_rank($$) {
    my ($left, $right) = @_;
    return $right->{n} <=> $left->{n}
        || $right->{whole} <=> $left->{whole}
        || $right->{id_cov} <=> $left->{id_cov}
        || $right->{has_desc} <=> $left->{has_desc}
        || $left->{hit} cmp $right->{hit};
}

sub rank_chain_references {
    my ($ordered_ref) = @_;
    my @ordered = @$ordered_ref;

    my %candidate_hits;
    for my $gene (@ordered) {
        for my $hit (keys %{ $gene_hits{$gene} }) {
            $candidate_hits{$hit} = 1;
        }
    }

    my @ranked;
    for my $hit (keys %candidate_hits) {
        my @with_hit;
        for my $gene (@ordered) {
            push @with_hit, $gene if exists $subjects{$hit}{$gene};
        }
        next unless @with_hit > 1;
        my $has_desc = ($desc{$hit} // "") =~ /[\s\|]/ ? 1 : 0;
        my $whole = 1;
        for my $gene (@with_hit) {
            my $qcov = $subjects{$hit}{$gene}{qcov};
            $whole = 0 if defined $qcov && $qcov ne "" && $qcov < $ref_min_qcov;
        }
        push @ranked, { hit => $hit, n => scalar @with_hit, has_desc => $has_desc, whole => $whole,
                        id_cov => identity_weighted_cov($hit, @with_hit) };
    }
    @ranked = sort by_reference_rank @ranked;
    return @ranked;
}

# ---------------------------------------------------------------------------
# Helper: PIECES_OVERLAP -- do two pieces of a merge align to the same region of
# some reference? Pieces of one split gene sit side by side on every reference;
# copies of one gene unit sit on top of each other. Uses copy_hsp (hits where a
# piece aligns over >= overlap_min_qcov of itself, fragment filter or not).
# Returns the first such reference (sorted ids), or undef.
# ---------------------------------------------------------------------------
sub pieces_overlap_ref {
    my @genes = @_;
    for my $i (0 .. $#genes - 1) {
        my $left_hits = $copy_hsp{ $genes[$i] } or next;
        for my $j ($i + 1 .. $#genes) {
            my $right_hits = $copy_hsp{ $genes[$j] } or next;
            for my $hit (sort keys %$left_hits) {
                my $right = $right_hits->{$hit} or next;
                my $left  = $left_hits->{$hit};
                my $from    = $left->{sstart} > $right->{sstart} ? $left->{sstart} : $right->{sstart};
                my $to      = $left->{send}   < $right->{send}   ? $left->{send}   : $right->{send};
                my $overlap = $to - $from + 1;
                next if $overlap <= 0;
                my $left_len  = $left->{send}  - $left->{sstart}  + 1;
                my $right_len = $right->{send} - $right->{sstart} + 1;
                my $shorter   = $left_len < $right_len ? $left_len : $right_len;
                return $hit if $overlap >= $overlap_min_frac * $shorter;
            }
        }
    }
    return undef;
}

# ---------------------------------------------------------------------------
# Helper: compute flags for a merge row
# ---------------------------------------------------------------------------
sub compute_flags {
    my ($ordered_ref, $tiling_pct, $overlap_ref) = @_;
    my @ordered = @$ordered_ref;
    my @flags;

    # LARGE_SPAN / LARGE_SPAN_EXTREME: total genomic footprint of the merged locus.
    # Calculated from the outermost coordinates of all genes in the chain.
    # A large span is biologically plausible for some gene families but should be
    # reviewed — very large spans combined with weak evidence are high-risk merges.
    my (@starts, @ends);
    for my $gene (@ordered) { push @starts, $locs{$gene}{start}; push @ends, $locs{$gene}{end} }
    my $chain_start = first_sorted_asc(@starts);
    my $chain_end   = first_sorted_desc(@ends);
    my $chain_span  = $chain_end - $chain_start + 1;
    if    ($chain_span > $large_span_extreme) { push @flags, "LARGE_SPAN_EXTREME" }
    elsif ($chain_span > $large_span_warn)    { push @flags, "LARGE_SPAN"         }

    # LOW_COV applies to all chain sizes
    push @flags, "LOW_COV" if $tiling_pct < $low_cov_thresh;

    # PIECES_OVERLAP: two pieces are copies of one gene unit (see pieces_overlap_ref)
    push @flags, "PIECES_OVERLAP" if defined $overlap_ref;

    # SKIPPED_GENE / OPPOSITE_STRAND_SKIP: any junction with genomic dist > 1.
    # If ALL skipped genes are on the opposite strand from the chain, flag as
    # OPPOSITE_STRAND_SKIP (interleaved unrelated genes; not a concern for the
    # merge itself).  If ANY skipped gene is on the same strand, flag as
    # SKIPPED_GENE (may be an additional split fragment — review before merging).
    my @skipped = get_skipped_genes(@ordered);
    if (@skipped) {
        my $chain_strand = $locs{$ordered[0]}{strand};
        my $any_same     = 0;
        for my $gene (@skipped) { $any_same++ if ($locs{$gene}{strand} // "") eq $chain_strand }
        push @flags, ($any_same ? "SKIPPED_GENE" : "OPPOSITE_STRAND_SKIP");
    }

    # MULTI_ISOFORM_JOIN: any source gene has more than one annotated transcript.
    # The merge script builds transcripts by cross-product (N*M combinations).
    # Not all cross-product transcripts are necessarily biologically real —
    # they require transcript-level evidence (IsoSeq isoform data) to validate.
    my $multi_isoform = 0;
    for my $gene (@ordered) { $multi_isoform++ if ($gene_n_transcripts{$gene} // 1) > 1 }
    push @flags, "MULTI_ISOFORM_JOIN" if $multi_isoform;

    # JUNCTION_OVERHANG: at some junction a piece extends more than end_tol aa past
    # its best alignment on the junction side -- sequence the reference lacks there
    # (a mispredicted exon at the break, or two different proteins)
    my $overhanging = 0;
    for my $i (0 .. $#ordered - 1) { $overhanging++ if ($pair_overhang{"$ordered[$i]\t$ordered[$i+1]"} // 0) > $end_tol }
    push @flags, "JUNCTION_OVERHANG" if $overhanging;

    # 2-gene pairs: only SINGLE_HIT possible beyond the above
    if (@ordered == 2) {
        my $hits = $pair_tiling{"$ordered[0]\t$ordered[1]"} // 0;
        push @flags, "SINGLE_HIT" if $hits == 1;
        return @flags ? join(",", @flags) : "";
    }

    # 3+ gene chains
    my @hits    = get_junction_hits(@ordered);
    my @numeric;
    my $n_unknown = 0;
    for my $hit_count (@hits) {
        if (defined $hit_count) { push @numeric, $hit_count } else { $n_unknown++ }
    }

    # TRANSITIVE_JOIN: one or more consecutive genes in this chain have no direct
    # pairwise tiling evidence between them (junction = "?" in junction_tiling_hits).
    # They were placed in the same chain only because each is independently connected
    # to another gene via a different reference protein — not because they themselves
    # share a common reference where their alignments tile end-to-end.
    # See FLAG DEFINITIONS for full description and biological implications.
    push @flags, "TRANSITIVE_JOIN" if $n_unknown;

    unless (@numeric) {
        return @flags ? join(",", @flags) : "";
    }

    my $min_hits = first_sorted_asc(@numeric);
    my $max_hits = first_sorted_desc(@numeric);

    if ($min_hits >= $strong_thresh) {
        push @flags, "STRONG";
    } elsif ($min_hits == 0 && $max_hits == 0) {
        # all unknown — no numeric flag to add
    } else {
        my $all_one  = 1;
        my $weak_end = 0;
        my $weak_int = 0;

        for my $i (0 .. $#hits) {
            my $h = $hits[$i] // 0;
            $all_one = 0 if $h != 1;
            next unless $h <= 2;

            my $is_terminal = ($i == 0 || $i == $#hits);
            if ($is_terminal) {
                # Flag as WEAK_END when the chain survived with a weak terminal:
                # - max_hits below threshold (trimming wouldn't fire anyway), OR
                # - h == 2 (2-hit terminal, never trimmed), OR
                # - no_asym_trim active (trimming disabled; flag all weak terminals
                #   including those that would have been trimmed in default mode)
                $weak_end = 1 if $max_hits < $asym_threshold || $h > 1 || $no_asym_trim;
            } else {
                $weak_int = 1;
            }
        }

        push @flags, "SINGLE_HIT"    if $all_one;
        push @flags, "WEAK_END"      if $weak_end && !$all_one;
        push @flags, "WEAK_INTERNAL" if $weak_int && !$all_one;
    }
    return @flags ? join(",", @flags) : "";
}

# ---------------------------------------------------------------------------
# Write merge table — apply trimming/splitting before writing each chain
# ---------------------------------------------------------------------------
open MERGE, ">split_genes_merge.txt" or die "cant open merge output $!\n";  # renamed to merge_candidates.txt by run_mender.pl
print MERGE join("\t",
    "merge_id",
    "num_genes",
    "genes_in_order",
    "coords_in_order",
    "ref",
    "gene_descs",
    "tiling_ref", "tiling_ref_desc",
    "tiling_cov_pct",
    "tiling_gap",
    "evalues",
    "junction_tiling_hits",
    "junction_genomic_dist",
    "skipped_genes",
    "max_tiling_hits",
    "num_pairs_in_chain",
    "flag",
    "rebuild_candidates",
    "overlap_ref"
), "\n";

my $merge_id   = 1;
my $n_trimmed  = 0;
my $n_split    = 0;
my $n_overlap  = 0;

# chains in genome order of their first gene
for my $chain (@chains) {
    my @genes = sort by_start_then_id keys %{ $chain->{genes} };
    $chain->{first_gene} = $genes[0];
}
sub by_chain_position($$) {
    my ($left, $right) = @_;
    my ($left_gene, $right_gene) = ($left->{first_gene}, $right->{first_gene});
    return $locs{$left_gene}{ref}   cmp $locs{$right_gene}{ref}
        || $locs{$left_gene}{start} <=> $locs{$right_gene}{start}
        || $left_gene cmp $right_gene;
}
for my $chain (sort by_chain_position @chains) {

    # genes in genomic order before refinement
    my @raw_ordered = sort by_start_then_id keys %{ $chain->{genes} };

    # apply trimming and splitting — may return multiple pieces
    my @pieces = refine_chain(@raw_ordered);
    $n_trimmed++ if scalar @pieces == 1 && scalar @{$pieces[0]} < scalar @raw_ordered;
    $n_split++   if scalar @pieces > 1;

    for my $piece (@pieces) {
        my @ordered = @$piece;
        my $ref     = $locs{ $ordered[0] }{ref};

        # rank the reference proteins for this piece: the rebuild candidates
        my @ranked_refs = rank_chain_references(\@ordered);
        my $last_ref = $#ranked_refs < $max_ref_candidates - 1 ? $#ranked_refs : $max_ref_candidates - 1;
        my @ranked_ids;
        for my $ranked (@ranked_refs[0 .. $last_ref]) { push @ranked_ids, $ranked->{hit} }

        # the best pair within this piece: its tiling reference, coverage and gap
        my %in_piece;
        for my $gene (@ordered) { $in_piece{$gene} = 1 }
        my @piece_pairs;
        for my $pair (@{ $chain->{pairs} }) {
            push @piece_pairs, $pair if $in_piece{ $pair->{g1} } && $in_piece{ $pair->{g2} };
        }

        my @pairs_by_pct = sort by_tiling_pct_descending (@piece_pairs ? @piece_pairs : @{ $chain->{pairs} });
        my $best_pair  = $pairs_by_pct[0];
        my $tiling_ref = $best_pair->{tiling_ref};

        my $max_tiling = $best_pair->{num_tiling};
        if (@piece_pairs) {
            my @by_tiling = sort by_num_tiling_descending @piece_pairs;
            $max_tiling = $by_tiling[0]{num_tiling};
        }

        my $gene_list  = join(",", @ordered);
        my (@coords, @descs, @evals);
        for my $gene (@ordered) {
            push @coords, $locs{$gene}{coord};
            push @descs,  $desc{$gene} // "?";
            push @evals,  $subjects{$tiling_ref}{$gene}{evalue} // "?";
        }
        my $coord_list = join(",", @coords);
        my $desc_list  = join(",", @descs);
        my $eval_list  = join(",", @evals);

        # junction tiling hits and genomic distances
        my @jhits    = get_junction_hits(@ordered);
        my @jdists   = get_junction_dists(@ordered);
        my $junc_str  = @ordered > 2
            ? join("|", unknown_as_question(@jhits))
            : "NA";
        my $jdist_str = @ordered > 1
            ? join("|", unknown_as_question(@jdists))
            : "NA";

        # skipped genes — those sitting between non-adjacent junction pairs
        my @skipped     = get_skipped_genes(@ordered);
        my $skipped_str = @skipped ? join(",", @skipped) : "none";

        my $tiling_pct  = $best_pair->{tiling_pct};
        my $tiling_gap  = $best_pair->{tiling_gap};
        my $overlap_ref = pieces_overlap_ref(@ordered);
        $n_overlap++ if defined $overlap_ref;
        my $flag = compute_flags(\@ordered, $tiling_pct, $overlap_ref);
        $flag = 'CLEAN' if $flag eq "";

        print MERGE join("\t",
            "merge_$merge_id",
            scalar @ordered,
            $gene_list,
            $coord_list,
            $ref,
            $desc_list,
            $tiling_ref,
            $desc{$tiling_ref} // "?",
            sprintf("%.1f", $tiling_pct),
            $tiling_gap,
            $eval_list,
            $junc_str,
            $jdist_str,
            $skipped_str,
            $max_tiling,
            scalar @ordered - 1,
            $flag,
            (@ranked_ids ? join(",", @ranked_ids) : $tiling_ref),
            $overlap_ref // "."
        ), "\n";

        $merge_id++;
    }
}

close MERGE;

print STDERR "Tiling tests skipped (pieces in the wrong order for the strand): $n_not_collinear\n";
print STDERR "Gene pairs reported: ",  scalar @all_pairs,       "\n";
print STDERR "Chains before refine: ", scalar @chains,          "\n";
print STDERR "Chains trimmed:       ", $n_trimmed,              "\n";
print STDERR "Chains split:         ", $n_split,                "\n";
print STDERR "Merge rows written:   ", $merge_id - 1,           "\n";
print STDERR "PIECES_OVERLAP rows:  ", $n_overlap,              "\n";
