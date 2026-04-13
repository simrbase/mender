#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## find_split_genes.pl
##
## Identifies chameleon gene annotations that may have been erroneously split
## during genome annotation. The core idea: if two or more adjacent chameleon
## genes each cover a non-overlapping portion of the same reference protein,
## they are candidates to be merged into a single gene.
##
## Evidence comes from diamond blastp of chameleon proteins vs a well-annotated
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
##                                  diamond.out genes.gff subject.fa query.fa
##
## Output files:
##
## split_genes_summary.txt
##   One row per gene pair that passes all filters. Best supporting hit only.
##   Columns:
##     g1, g2             : the two chameleon gene IDs (genomic order)
##     g1_coord, g2_coord : chromosome.start.end for each gene
##     genomic_dist       : number of rank positions between g1 and g2
##                          (1 = directly adjacent, 2 = one gene between them, etc)
##     g1_desc, g2_desc   : descriptions from query fasta headers
##     best_hit           : reference protein with highest combined coverage
##                          pct for this specific gene pair. Selected purely
##                          by combined span / subject length — no description
##                          preference applied here.
##     hit_desc           : description of the best_hit protein
##     combined_cov_pct   : % of best_hit covered by g1+g2 alignments together
##                          (span from min sstart to max send / subject length)
##     tiling_gap         : aa between the two alignments on best_hit
##                          positive = gap between alignments (they don't quite meet)
##                          zero     = alignments meet perfectly end to end
##                          negative = alignments overlap slightly
##     g1_subj, g2_subj   : alignment coordinates on the subject (sstart-send)
##     g1_cov, g2_cov     : individual scovhsp% for each gene vs best_hit
##     g1_evalue,g2_evalue: diamond evalues for each gene vs best_hit
##     num_tiling_hits    : how many different reference proteins support this
##                          pair as tiling fragments — higher = more confident
##
## split_genes_detail.txt
##   One row per supporting tiling hit per gene pair, sorted by coverage desc.
##   Same columns as summary but for every tiling hit, not just the best.
##   Use this to see the full evidence for any pair from the summary.
##   Columns: g1, g2, genomic_dist, hit, hit_desc, combined_cov_pct,
##            tiling_gap, g1_subj, g2_subj, g1_cov, g2_cov,
##            g1_evalue, g2_evalue
##
## split_genes_merge.txt
##   One row per merge candidate. Adjacent pairs sharing a gene are chained
##   into multi-gene groups (e.g. A-B + B-C becomes A,B,C in one row).
##   Chains are trimmed and split based on junction tiling hit asymmetry
##   before writing — see CHAIN TRIMMING AND SPLITTING below.
##   The best_hit is chosen as the reference protein with the most chain genes
##   having valid entries against it, breaking ties by combined coverage span,
##   with preference for hits that have a real description over bare Ensembl IDs.
##   Columns:
##     merge_id              : unique identifier
##     num_genes             : how many chameleon genes to merge
##     genes_in_order        : comma-sep gene IDs in genomic order
##     coords_in_order       : comma-sep coordinates in genomic order
##     ref                   : chromosome
##     gene_descs            : comma-sep descriptions in genomic order
##     best_hit              : reference protein best representing the whole chain
##     hit_desc              : description of best_hit
##     best_combined_cov_pct : combined coverage of best_hit by the best pair
##     best_gap              : tiling gap for the best pair
##     evalues               : comma-sep evalues in genomic order vs best_hit
##                             '?' means that gene has no hit vs best_hit
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
##
## =============================================================================
## KEY CONCEPTS
## =============================================================================
##
## COMBINED_COV_PCT
##   Brief: What fraction of the reference protein is spanned by the two genes
##          together.
##
##   Detail: This is NOT the sum of the two individual coverages. It is the
##           span from the leftmost alignment start to the rightmost alignment
##           end on the reference protein, divided by the reference protein
##           length. Any gap between the two alignments is included in the span,
##           so combined_cov_pct slightly overstates true coverage when there
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
##     combined_cov_pct = 471 / 500 * 100 = 94.2%
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
##     best_hit = whichever of A/B/C gives highest combined_cov_pct
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
##          all vs the same chain best_hit reference protein.
##
##   Detail: These come from your diamond blastp search of chameleon proteins
##           vs the reference proteome — NOT from the evalue strings embedded
##           in your fasta headers (those are from a different search).
##           A '?' means that gene had no diamond hit vs the chain best_hit
##           that passed the evalue/coverage filters. This can happen when the
##           chain best_hit was chosen because it covers most of the chain but
##           one terminal gene tiles against a slightly different reference
##           protein. Check split_genes_detail.txt for that gene pair to see
##           what it does tile against.
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
##     LOW_COV            : combined_cov_pct < 60%. The genes together cover
##                          less than 60% of the reference protein, suggesting
##                          domain sharing rather than a genuine split gene.
##                          Recommended to skip_flags in production runs.
##
##     SKIPPED_GENE       : a non-adjacent gene sits inside the merge locus
##                          (genomic_dist > 1 for at least one junction). The
##                          skipped gene may be an additional split fragment that
##                          failed filters, a different gene, or an artifact.
##                          Check skipped_genes column and GFF before merging.
##                          Recommended to skip_flags and review manually.
##
##     MULTI_ISOFORM_JOIN : at least one source gene has more than one annotated
##                          transcript. The merged gene will be built by cross-
##                          product (N*M transcript combinations). Not all
##                          cross-product transcripts are necessarily real —
##                          they require transcript-level evidence to validate.
##                          Review isoform structure in a genome browser before
##                          relying on isoform-level annotations from this merge.
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
##                          Recommended action: inspect hit_desc and pident. If
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
##                          "?" in the evalues column vs the chain best_hit, that
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
##     LOW_COV            : best_combined_cov_pct < 60%
##
##     SKIPPED_GENE       : at least one junction spans a non-adjacent gene.
##
##   Interpreting combinations:
##     STRONG                         -> act on this merge
##     STRONG,SKIPPED_GENE            -> strong evidence but check skipped gene;
##                                       it may need to be added to the merge
##     STRONG,TRANSITIVE_JOIN         -> inspect the ? junction carefully;
##                                       require IsoSeq if possible
##     SINGLE_HIT                     -> inspect hit_desc and pident; see above
##     SINGLE_HIT,TRANSITIVE_JOIN     -> high risk; require IsoSeq FULL_SPAN
##     WEAK_END,LOW_COV               -> low confidence; review carefully
##     STRONG,LOW_COV                 -> protein evidence consistent but incomplete;
##                                       may be domain sharing, not a split gene
##     LARGE_SPAN,SINGLE_HIT          -> high FP risk; require IsoSeq
##     LARGE_SPAN_EXTREME             -> add to skip_flags unless IsoSeq confirms
##     MULTI_ISOFORM_JOIN             -> merge proceeds; verify transcript structure
##
## =============================================================================

my $no_asym_trim       = 0;
my $large_span_warn    = 500000;   # bp: flag LARGE_SPAN above this
my $large_span_extreme = 2000000;  # bp: flag LARGE_SPAN_EXTREME above this

GetOptions(
    "no_asym_trim"         => \$no_asym_trim,
    "large_span_warn=i"    => \$large_span_warn,
    "large_span_extreme=i" => \$large_span_extreme,
) or die "usage: $0 [--no_asym_trim] [--large_span_warn N] [--large_span_extreme N] blast gff subject_fa query_fa\n";

my $usage_str = "usage: $0 [--no_asym_trim] [--large_span_warn N] [--large_span_extreme N] blast gff subject_fa query_fa\n";
my $blast      = shift or die $usage_str;
my $gff        = shift or die $usage_str;
my $subject_fa = shift or die $usage_str;
my $query_fa   = shift or die $usage_str;

my $wiggle         = 15;  # aa wiggle room for tiling junction
my $max_dist       =  4;  # max genomic rank distance between gene pair
my $asym_threshold =  6;  # min hits on strong side to trigger trim/split
my $low_cov_thresh = 60;  # combined_cov_pct below this gets LOW_COV flag
my $strong_thresh  =  3;  # min junction hits for STRONG flag
my $frag_scov_max  = 85;  # scovhsp <= this is always a fragment candidate
my $frag_qcov_min  = 90;  # qcovhsp >= this enables the soft scovhsp extension
my $frag_scov_soft = 90;  # with qcovhsp >= frag_qcov_min, keep if scovhsp < this

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
    my $joined   = join(" | ", grep { defined $_ } $symbol, $d);
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

# Sort all genes by chromosome then start position
my @by_start = sort {
    $locs{$a}{ref}   cmp $locs{$b}{ref} ||
    $locs{$a}{start} <=> $locs{$b}{start}
} keys %locs;

# ---------------------------------------------------------------------------
# Identify nested genes: genes whose coordinates fall entirely within another
# gene. These are excluded from the rank array because they would inflate the
# apparent genomic distance between flanking genes.
# Uses a sliding window lookback of 10 positions since nesting is always local.
# ---------------------------------------------------------------------------
my %nested;  # nested_gene -> parent_gene
for my $i (0 .. $#by_start) {
    my $gene     = $by_start[$i];
    my $lookback = $i < 10 ? 0 : $i - 10;
    for my $j ($lookback .. $i-1) {
        my $parent = $by_start[$j];
        next unless $locs{$gene}{ref} eq $locs{$parent}{ref};
        next if $locs{$parent}{end} < $locs{$gene}{start};
        if ($locs{$gene}{end} <= $locs{$parent}{end}) {
            $nested{$gene} = $parent;
            last;
        }
    }
}

# Build rank array and rank lookup excluding nested genes
my @all_sorted = grep { !exists $nested{$_} } @by_start;
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
# ---------------------------------------------------------------------------
my %subjects;
my %gene_hits;
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

    if ($evalue <= 1e-20 && $is_fragment) {
        # keep only the best-coverage HSP per gene-subject pair
        if (!exists $subjects{$subject}{$gene_id} ||
             $subjects{$subject}{$gene_id}{cov} < $scovhsp){
            $subjects{$subject}{$gene_id}{cov}    = $scovhsp;
            $subjects{$subject}{$gene_id}{sstart} = $sstart;
            $subjects{$subject}{$gene_id}{send}   = $send;
            $subjects{$subject}{$gene_id}{evalue} = $evalue;
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
    "best_hit", "hit_desc",
    "combined_cov_pct",
    "tiling_gap",
    "g1_subj", "g2_subj",
    "g1_cov", "g2_cov",
    "g1_evalue", "g2_evalue",
    "num_tiling_hits"
), "\n";

print DETAIL join("\t",
    "g1", "g2",
    "genomic_dist",
    "hit", "hit_desc",
    "combined_cov_pct",
    "tiling_gap",
    "g1_subj", "g2_subj",
    "g1_cov", "g2_cov",
    "g1_evalue", "g2_evalue"
), "\n";

# Collect all reported pairs for downstream chaining into merge candidates
# Also build a lookup: "g1\tg2" -> {num_tiling_hits, genomic_dist}
my @all_pairs;
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

        my @shared = grep { exists $big->{$_} } keys %$small;
        next unless @shared;

        my $best_hit      = undef;
        my $best_pct      = 0;
        my $best_gap      = undef;
        my $best_g1_range = undef;
        my $best_g2_range = undef;
        my $best_g1_cov   = undef;
        my $best_g2_cov   = undef;
        my $best_g1_eval  = undef;
        my $best_g2_eval  = undef;
        my $num_tiling    = 0;
        my @tiling_hits;  # all hits that pass tiling filter, for detail output

        for my $hit (@shared) {
            my $s1 = $subjects{$hit}{$g1}{sstart};
            my $e1 = $subjects{$hit}{$g1}{send};
            my $s2 = $subjects{$hit}{$g2}{sstart};
            my $e2 = $subjects{$hit}{$g2}{send};

            # Tiling gap: how far apart are the two alignments on the subject?
            # Positive = gap between them, negative = they overlap slightly.
            # We allow abs(gap) <= wiggle to account for alignment fuzziness.
            my $gap = ($s1 < $s2) ? $s2 - $e1 : $s1 - $e2;
            next if abs($gap) > $wiggle;

            $num_tiling++;

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
                gap    => $gap,
                g1_r   => "$s1-$e1",
                g2_r   => "$s2-$e2",
                g1_cov => $subjects{$hit}{$g1}{cov},
                g2_cov => $subjects{$hit}{$g2}{cov},
                g1_ev  => $subjects{$hit}{$g1}{evalue},
                g2_ev  => $subjects{$hit}{$g2}{evalue},
            };

            if ($pct > $best_pct) {
                $best_pct      = $pct;
                $best_hit      = $hit;
                $best_gap      = $gap;
                $best_g1_range = "$s1-$e1";
                $best_g2_range = "$s2-$e2";
                $best_g1_cov   = $subjects{$hit}{$g1}{cov};
                $best_g2_cov   = $subjects{$hit}{$g2}{cov};
                $best_g1_eval  = $subjects{$hit}{$g1}{evalue};
                $best_g2_eval  = $subjects{$hit}{$g2}{evalue};
            }
        }

        next unless defined $best_hit;

        my $gdist    = $j - $i;
        my $g1_desc  = $desc{$g1} // "?";
        my $g2_desc  = $desc{$g2} // "?";
        my $hit_desc = $desc{$best_hit} // "?";

        print SUMMARY join("\t",
            $g1, $g2,
            $locs{$g1}{coord}, $locs{$g2}{coord},
            $gdist,
            $g1_desc, $g2_desc,
            $best_hit, $hit_desc,
            sprintf("%.1f", $best_pct),
            $best_gap,
            $best_g1_range, $best_g2_range,
            $best_g1_cov, $best_g2_cov,
            $best_g1_eval, $best_g2_eval,
            $num_tiling
        ), "\n";

        # Detail: all tiling hits sorted by combined coverage descending
        for my $th (sort { $b->{pct} <=> $a->{pct} } @tiling_hits) {
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
            best_hit   => $best_hit,
            hit_desc   => $hit_desc,
            best_pct   => $best_pct,
            best_gap   => $best_gap,
            g1_cov     => $best_g1_cov,
            g2_cov     => $best_g2_cov,
            g1_eval    => $best_g1_eval,
            g2_eval    => $best_g2_eval,
            num_tiling => $num_tiling,
        };
        $pair_tiling{"$g1\t$g2"} = $num_tiling;
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

    my @matching_idx = grep {
        $chains[$_]{genes}{$g1} || $chains[$_]{genes}{$g2}
    } 0 .. $#chains;

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
        my $max_hits = (sort { $b <=> $a } grep { defined $_ } @hits)[0] // 0;

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
    return grep { scalar @$_ >= 2 } @pieces;
}

# ---------------------------------------------------------------------------
# Helper: select best_hit for a chain — reference protein with most chain
# genes covered, preferring hits with real descriptions over bare Ensembl IDs,
# breaking ties by combined span
# ---------------------------------------------------------------------------
sub select_chain_best_hit {
    my ($ordered_ref) = @_;
    my @ordered = @$ordered_ref;

    my %candidate_hits;
    for my $gene (@ordered) {
        for my $hit (keys %{ $gene_hits{$gene} }) {
            $candidate_hits{$hit} = 1;
        }
    }

    my $best_hit      = undef;
    my $best_n        = 0;
    my $best_span     = 0;
    my $best_has_desc = 0;

    for my $hit (keys %candidate_hits) {
        my @with_hit = grep { exists $subjects{$hit}{$_} } @ordered;
        my $n = scalar @with_hit;
        next unless $n > 1;

        my $min_s = (sort { $a <=> $b } map { $subjects{$hit}{$_}{sstart} } @with_hit)[0];
        my $max_e = (sort { $b <=> $a } map { $subjects{$hit}{$_}{send}   } @with_hit)[0];
        my $span  = $max_e - $min_s + 1;
        my $has_desc = ($desc{$hit} // "") =~ /[\s\|]/ ? 1 : 0;

        if (!defined $best_hit) {
            $best_hit = $hit; $best_n = $n; $best_span = $span; $best_has_desc = $has_desc;
        } elsif ($has_desc && !$best_has_desc) {
            $best_hit = $hit; $best_n = $n; $best_span = $span; $best_has_desc = $has_desc;
        } elsif (!$has_desc && $best_has_desc) {
            next;
        } elsif ($n > $best_n || ($n == $best_n && $span > $best_span)) {
            $best_hit = $hit; $best_n = $n; $best_span = $span; $best_has_desc = $has_desc;
        }
    }
    return $best_hit;
}

# ---------------------------------------------------------------------------
# Helper: compute flags for a merge row
# ---------------------------------------------------------------------------
sub compute_flags {
    my ($ordered_ref, $best_pct) = @_;
    my @ordered = @$ordered_ref;
    my @flags;

    # LARGE_SPAN / LARGE_SPAN_EXTREME: total genomic footprint of the merged locus.
    # Calculated from the outermost coordinates of all genes in the chain.
    # A large span is biologically plausible for some gene families but should be
    # reviewed — very large spans combined with weak evidence are high-risk merges.
    my $chain_start = (sort { $a <=> $b } map { $locs{$_}{start} } @ordered)[0];
    my $chain_end   = (sort { $b <=> $a } map { $locs{$_}{end}   } @ordered)[0];
    my $chain_span  = $chain_end - $chain_start + 1;
    if    ($chain_span > $large_span_extreme) { push @flags, "LARGE_SPAN_EXTREME" }
    elsif ($chain_span > $large_span_warn)    { push @flags, "LARGE_SPAN"         }

    # LOW_COV applies to all chain sizes
    push @flags, "LOW_COV" if $best_pct < $low_cov_thresh;

    # SKIPPED_GENE: any junction with genomic dist > 1
    my @skipped = get_skipped_genes(@ordered);
    push @flags, "SKIPPED_GENE" if @skipped;

    # MULTI_ISOFORM_JOIN: any source gene has more than one annotated transcript.
    # The merge script builds transcripts by cross-product (N*M combinations).
    # Not all cross-product transcripts are necessarily biologically real —
    # they require transcript-level evidence (IsoSeq isoform data) to validate.
    push @flags, "MULTI_ISOFORM_JOIN"
        if grep { ($gene_n_transcripts{$_} // 1) > 1 } @ordered;

    # 2-gene pairs: only SINGLE_HIT possible beyond the above
    if (@ordered == 2) {
        my $hits = $pair_tiling{"$ordered[0]\t$ordered[1]"} // 0;
        push @flags, "SINGLE_HIT" if $hits == 1;
        return @flags ? join(",", @flags) : "";
    }

    # 3+ gene chains
    my @hits    = get_junction_hits(@ordered);
    my @numeric = grep { defined $_ } @hits;

    # TRANSITIVE_JOIN: one or more consecutive genes in this chain have no direct
    # pairwise tiling evidence between them (junction = "?" in junction_tiling_hits).
    # They were placed in the same chain only because each is independently connected
    # to another gene via a different reference protein — not because they themselves
    # share a common reference where their alignments tile end-to-end.
    # See FLAG DEFINITIONS for full description and biological implications.
    push @flags, "TRANSITIVE_JOIN" if grep { !defined $_ } @hits;

    unless (@numeric) {
        return @flags ? join(",", @flags) : "";
    }

    my $min_hits = (sort { $a <=> $b } @numeric)[0];
    my $max_hits = (sort { $b <=> $a } @numeric)[0];

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
open MERGE, ">split_genes_merge.txt" or die "cant open merge output $!\n";
print MERGE join("\t",
    "merge_id",
    "num_genes",
    "genes_in_order",
    "coords_in_order",
    "ref",
    "gene_descs",
    "best_hit", "hit_desc",
    "best_combined_cov_pct",
    "best_gap",
    "evalues",
    "junction_tiling_hits",
    "junction_genomic_dist",
    "skipped_genes",
    "max_tiling_hits",
    "num_pairs_in_chain",
    "flag"
), "\n";

my $merge_id   = 1;
my $n_trimmed  = 0;
my $n_split    = 0;

for my $chain (sort {
    my $a_gene = (sort { $locs{$a}{start} <=> $locs{$b}{start} } keys %{$a->{genes}})[0];
    my $b_gene = (sort { $locs{$a}{start} <=> $locs{$b}{start} } keys %{$b->{genes}})[0];
    $locs{$a_gene}{ref} cmp $locs{$b_gene}{ref} ||
    $locs{$a_gene}{start} <=> $locs{$b_gene}{start}
} @chains) {

    # genes in genomic order before refinement
    my @raw_ordered = sort { $locs{$a}{start} <=> $locs{$b}{start} }
                      keys %{ $chain->{genes} };

    # apply trimming and splitting — may return multiple pieces
    my @pieces = refine_chain(@raw_ordered);
    $n_trimmed++ if scalar @pieces == 1 && scalar @{$pieces[0]} < scalar @raw_ordered;
    $n_split++   if scalar @pieces > 1;

    for my $piece (@pieces) {
        my @ordered = @$piece;
        my $ref     = $locs{ $ordered[0] }{ref};

        # select best hit for this piece
        my $chain_best_hit = select_chain_best_hit(\@ordered);

        # fall back to best pairwise hit if needed
        unless (defined $chain_best_hit) {
            my $best_pair = (sort { $b->{best_pct} <=> $a->{best_pct} }
                             @{ $chain->{pairs} })[0];
            $chain_best_hit = $best_pair->{best_hit};
        }

        # find best pairwise hit within this piece for coverage/gap reporting
        my @piece_pairs = grep {
            my $g1 = $_->{g1}; my $g2 = $_->{g2};
            grep { $_ eq $g1 } @ordered and grep { $_ eq $g2 } @ordered
        } @{ $chain->{pairs} };

        my $best_pair = @piece_pairs
            ? (sort { $b->{best_pct} <=> $a->{best_pct} } @piece_pairs)[0]
            : (sort { $b->{best_pct} <=> $a->{best_pct} } @{ $chain->{pairs} })[0];

        my $max_tiling = @piece_pairs
            ? (sort { $b->{num_tiling} <=> $a->{num_tiling} } @piece_pairs)[0]{num_tiling}
            : $best_pair->{num_tiling};

        my $gene_list  = join(",", @ordered);
        my $coord_list = join(",", map { $locs{$_}{coord}  } @ordered);
        my $desc_list  = join(",", map { $desc{$_} // "?"  } @ordered);
        my $eval_list  = join(",", map {
            $subjects{$chain_best_hit}{$_}{evalue} // "?"
        } @ordered);

        # junction tiling hits and genomic distances
        my @jhits    = get_junction_hits(@ordered);
        my @jdists   = get_junction_dists(@ordered);
        my $junc_str  = @ordered > 2
            ? join("|", map { $_ // "?" } @jhits)
            : "NA";
        #my $jdist_str = @ordered > 2
        #    ? join("|", map { $_ // "?" } @jdists)
        #    : "NA";
        my $jdist_str = @ordered > 1
            ? join("|", map { $_ // "?" } @jdists)
            : "NA";

        # skipped genes — those sitting between non-adjacent junction pairs
        my @skipped     = get_skipped_genes(@ordered);
        my $skipped_str = @skipped ? join(",", @skipped) : "none";

        my $best_pct = $best_pair->{best_pct};
        my $best_gap = $best_pair->{best_gap};
        my $flag     = compute_flags(\@ordered, $best_pct) ;
        $flag = 'CLEAN' if $flag eq "";
      

        print MERGE join("\t",
            "merge_$merge_id",
            scalar @ordered,
            $gene_list,
            $coord_list,
            $ref,
            $desc_list,
            $chain_best_hit,
            $desc{$chain_best_hit} // "?",
            sprintf("%.1f", $best_pct),
            $best_gap,
            $eval_list,
            $junc_str,
            $jdist_str,
            $skipped_str,
            $max_tiling,
            scalar @ordered - 1,
            $flag
        ), "\n";

        $merge_id++;
    }
}

close MERGE;

print STDERR "Gene pairs reported: ",  scalar @all_pairs,       "\n";
print STDERR "Chains before refine: ", scalar @chains,          "\n";
print STDERR "Chains trimmed:       ", $n_trimmed,              "\n";
print STDERR "Chains split:         ", $n_split,                "\n";
print STDERR "Merge rows written:   ", $merge_id - 1,           "\n";
