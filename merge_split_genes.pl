#!/usr/bin/perl
use strict;
use warnings;

## merge_split_genes.pl
##
## Merges erroneously split gene annotations in a GFF3 file based on the
## output of find_split_genes.pl + validate_with_isoseq.pl.
##
## For each merge candidate, the source genes are replaced by a single new
## gene whose transcripts are built by joining the exon/CDS/UTR features
## of the source transcripts in genomic order. The source gene entries are
## removed from the output GFF.
##
## Transcript joining logic:
##   Each source gene contributes one or more transcripts. All transcripts
##   from all source genes are combined by cross-product — every transcript
##   from gene A is joined with every transcript from gene B (and C, etc).
##   For the common case of 1 transcript per gene, this produces exactly
##   one merged transcript. For multi-isoform genes it produces N*M transcripts.
##   All sub-features (exon, CDS, five_prime_UTR, three_prime_UTR) are
##   collected from all source transcripts in a pair and sorted by coordinate.
##   CDS phase values are recalculated from scratch after joining.
##
## New gene ID scheme:
##   New genes are assigned IDs starting at CCA3g100001000.1, incrementing
##   by 1000 (CCA3g100002000.1, CCA3g100003000.1, etc). The script scans
##   the input GFF first to find the highest existing CCA3g1XXXXX000.1 ID
##   and starts the counter above that.
##
## Usage:
##   perl merge_split_genes.pl [options] isoseq_validated.txt input.gff output.gff
##
## Arguments:
##   isoseq_validated.txt  output of validate_with_isoseq.pl (or merge_candidates.txt)
##   input.gff            full GFF3 annotation (must be parent-before-child sorted)
##   output.gff           name for the new merged GFF3 output file
##
## Options:
##   --flags FLAGLIST   comma-sep list of flags to process (default: all)
##                      e.g. --flags STRONG,CLEAN
##                      Use 'all' to process every row regardless of flag
##   --min_tiling N     only process merges with max_tiling_hits >= N (default: 1)
##   --min_cov N        only process merges with combined_cov_pct >= N (default: 0)
##   --skip_isoseq F    skip merges with isoseq_flag matching F
##                      e.g. --skip_isoseq PARTIAL_SPAN
##   --removed FILE     filename for removed genes GFF (default: removed_genes.gff)
##   --dry_run          report what would be merged without writing any files
##   --date YYYY-MM-DD  date string for merge_date attribute (default: today)
##
## Output files:
##   output.gff           full annotation with merged genes replacing source genes
##   removed_genes.gff    all removed source genes and their children as audit trail
##                        each gene block is preceded by a comment showing which
##                        merge_id it came from and what it was merged with
##   STDERR               summary statistics
##
## New gene attributes:
##   ID=CCA3g100001000.1
##   Name=CCA3g100001000.1
##   merged_from=CCA3g007692000.1;CCA3g007693000.1
##   merge_source=computational
##   merge_date=YYYY-MM-DD
##
## New transcript attributes:
##   ID=CCA3t100001001.1
##   Parent=CCA3g100001000.1
##   Name=CCA3t100001001.1
##   merged_from=CCA3t007692001.1;CCA3t007693001.1  (source transcripts joined)
##   merge_source=computational
##
## Sub-feature renaming:
##   exon:          CCA3t100001001.1.exon.1, .exon.2, ...
##   CDS:           CCA3t100001001.1.CDS.1,  .CDS.2,  ...
##   five_prime_UTR:  CCA3t100001001.1.five_prime_UTR.1, ...
##   three_prime_UTR: CCA3t100001001.1.three_prime_UTR.1, ...
##
## CDS phase recalculation:
##   After joining CDS features from multiple source transcripts, phases are
##   recalculated from the start of the CDS (phase 0 for first CDS feature,
##   then tracking frame across features). This is necessary because the
##   original phase values assumed the CDS started at that gene's fragment.
##   NOTE: this recalculation assumes the merged transcript is in-frame end
##   to end. If the split happened mid-codon the recalculated phases may be
##   off by 1-2. Manual review of merged CDS structure is recommended.
##
## Notes on input GFF:
##   The GFF must be sorted so that parent features appear before their
##   children (gene before mRNA, mRNA before exon/CDS/UTR). This is standard
##   GFF3. If your GFF is not sorted, use one of:
##     gt gff3 -sort -tidy -retainids your.gff > sorted.gff
##     agat_convert_sp_gxf2gxf.pl -g your.gff -o sorted.gff
##   The script will exit with an error if a child feature is encountered
##   before its parent.
##
## =============================================================================
## FILTER FLAGS AND ARGUMENTS
## =============================================================================
##
## --skip_flags accepts any comma-separated combination of these values:
##
##   STRONG               all junctions have >= 3 tiling hits (you would not skip this)
##   SINGLE_HIT           all junctions have exactly 1 tiling hit.
##                        This is a REVIEW flag, not a default skip. For single-copy
##                        genes in the reference proteome, 1 hit is the expected maximum.
##                        Inspect hit_desc and pident before deciding. See FLAG DEFINITIONS
##                        in find_split_genes.pl for full guidance.
##   WEAK_END             terminal junction has 1-2 hits but survived trimming
##   WEAK_INTERNAL        internal junction has 1-2 hits but survived splitting
##   LOW_COV              combined coverage < 60% — likely domain sharing, not split gene
##   SKIPPED_GENE         a non-adjacent gene sits inside the merge locus — review manually
##   TRANSITIVE_JOIN      one or more consecutive genes in the chain have no direct
##                        pairwise tiling evidence; connection is inferred transitively.
##                        Compounds other flags: STRONG,TRANSITIVE_JOIN should be treated
##                        with the same caution as WEAK_END. Always check IsoSeq.
##   MULTI_ISOFORM_JOIN   at least one source gene has multiple transcripts; merged gene
##                        will contain cross-product isoform combinations.
##   LARGE_SPAN           merged locus exceeds large_span_warn genomic span (default 500kb)
##   LARGE_SPAN_EXTREME   merged locus exceeds large_span_extreme (default 2Mb) —
##                        strongly recommended to add to skip_flags unless IsoSeq confirms
##   CLEAN                2+ tiling hits, coverage >= 60%, directly adjacent — no concerns flagged
##
## --isoseq_flag values (for reference, used in isoseq_validated.txt col 20):
##   FULL_SPAN       at least one read spans all genes — strong confirmation
##   PARTIAL_SPAN    reads exist but don't reach a terminal gene — handled by --fix_partial
##   NO_SPANNERS     no spanning reads — may reflect expression stage not gene structure
##
## --gene_template / --trans_template bracket notation:
##   [GCOUNT:N]   gene counter, zero-padded to N digits, incremented per new gene
##   [TCOUNT:N]   transcript counter, zero-padded to N digits, incremented per transcript
##   :N is optional — omit for no zero-padding
##   All other text is literal
##
##   Examples:
##     CCA3 style:
##       --gene_template  "CCA3g1[GCOUNT:6]000"
##       --trans_template "CCA3t1[GCOUNT:6][TCOUNT:3]"
##       generates: CCA3g100001000 / CCA3t100001001
##
##     Helixer bat genome (ACI1_HiC_scaffold_1_000001):
##       --gene_template  "ACI1_HiC_scaffold_1_[GCOUNT:6]"
##       --trans_template "ACI1_HiC_scaffold_1_[GCOUNT:6].[TCOUNT]"
##       generates: ACI1_HiC_scaffold_1_000001 / ACI1_HiC_scaffold_1_000001.1
##
##     Medicago (MS.gene00001 / MS.gene00001.t1):
##       --gene_template  "MS.gene[GCOUNT:5]"
##       --trans_template "MS.gene[GCOUNT:5].t[TCOUNT]"
##       generates: MS.gene00001 / MS.gene00001.t1
##
##     Trinity (TRINITY_DN1000_c115_g5 / TRINITY_DN1000_c115_g5_i1):
##       --gene_template  "TRINITY_DN1000_c115_g[GCOUNT]"
##       --trans_template "TRINITY_DN1000_c115_g[GCOUNT]_i[TCOUNT]"
##       generates: TRINITY_DN1000_c115_g1 / TRINITY_DN1000_c115_g1_i1
##
##   Default (date-stamped, safe for any organism):
##       gene:        MERGE20260410g000001000
##       transcript:  MERGE20260410t000001001
##
##   Note: version suffix (e.g. .1) is intentionally excluded from templates.
##   New genes are always version 1. Versions are incremented only when
##   manually editing exon boundaries after the initial merge.
##
## =============================================================================
## EXAMPLE USAGE
## =============================================================================
##
## Dry run — see what would be merged without writing any files:
##
##   perl merge_split_genes.pl --fix_partial --skip_flags SKIPPED_GENE,LOW_COV \
##       --dry_run \
##       isoseq_validated.txt \
##       input.gff \
##       output.gff
##
## Recommended run — all candidates except SKIPPED_GENE and LOW_COV,
## with PARTIAL_SPAN terminal trimming, using CCA3 ID format:
##
##   perl merge_split_genes.pl \
##       --fix_partial \
##       --skip_flags SKIPPED_GENE,LOW_COV \
##       --gene_template  "CCA3g1[GCOUNT:6]000" \
##       --trans_template "CCA3t1[GCOUNT:6][TCOUNT:3]" \
##       isoseq_validated.txt \
##       input.gff \
##       output.gff
##
## Conservative run — only IsoSeq-confirmed merges, STRONG protein evidence,
## skip PARTIAL_SPAN, SKIPPED_GENE, LOW_COV:
##
##   perl merge_split_genes.pl \
##       --flags STRONG \
##       --skip_flags SKIPPED_GENE,LOW_COV,SINGLE_HIT,WEAK_END \
##       --require_isoseq FULL_SPAN \
##       --isoseq_min_spanning 3 \
##       --gene_template  "CCA3g1[GCOUNT:6]000" \
##       --trans_template "CCA3t1[GCOUNT:6][TCOUNT:3]" \
##       isoseq_validated.txt \
##       input.gff \
##       output.gff
##
## Run without IsoSeq validation (protein evidence only):
##
##   perl merge_split_genes.pl \
##       --flags STRONG \
##       --skip_flags SKIPPED_GENE,LOW_COV \
##       --gene_template  "ACI1_HiC_scaffold_1_[GCOUNT:6]" \
##       --trans_template "ACI1_HiC_scaffold_1_[GCOUNT:6].[TCOUNT]" \
##       merge_candidates.txt \
##       input.gff \
##       output.gff
##
## Run using default date-stamped IDs (safe for any organism):
##
##   perl merge_split_genes.pl \
##       --fix_partial \
##       --skip_flags SKIPPED_GENE,LOW_COV \
##       isoseq_validated.txt \
##       input.gff \
##       output.gff
##
## After any run, sanity check gene counts:
##
##   grep -c $'\tgene\t' input.gff
##   grep -c $'\tgene\t' output.gff
##   grep -c $'\tgene\t' removed_genes.gff
##   # output = input - removed + merged_count
##
##   # find all Mender-created features
##   grep -P '\tMend\t' output.gff | grep -P '\tgene\t' | wc -l
##
## =============================================================================


use POSIX qw(strftime);
use feature "state";
use Getopt::Long;

my $flags_opt         = "all";
my $skip_flags_opt    = "";
my $min_tiling        = 1;
my $min_cov           = 0;
my $require_isoseq    = "";
my $isoseq_min_span   = 0;
my $fix_partial       = 0;
my $removed_file      = "removed_genes.gff";
my $dry_run           = 0;
my $date              = strftime("%Y-%m-%d", localtime);

# New gene ID template options
# Templates use bracket notation to define ID structure:
#   [GCOUNT:N]  - gene counter, zero-padded to N digits, incremented per new gene
#   [TCOUNT:N]  - transcript counter, zero-padded to N digits, incremented per transcript
#   All other text is literal and included as-is in the ID
#
# Default templates use today's date so IDs are unique and datestamped:
#   gene:        MERGE20260410g000001000
#   transcript:  MERGE20260410t000001001
#
# Examples for common annotation styles:
#   CCA3 style:
#     --gene_template    "CCA3g1[GCOUNT:6]000"
#     --trans_template   "CCA3t1[GCOUNT:6][TCOUNT:3]"
#
#   Helixer bat genome style (e.g. ACI1_HiC_scaffold_1_000001):
#     --gene_template    "ACI1_HiC_scaffold_1_[GCOUNT:6]"
#     --trans_template   "ACI1_HiC_scaffold_1_[GCOUNT:6].[TCOUNT]"
#
#   Medicago style (MS.gene00001 / MS.gene00001.t1):
#     --gene_template    "MS.gene[GCOUNT:5]"
#     --trans_template   "MS.gene[GCOUNT:5].t[TCOUNT]"
#
#   Trinity style (TRINITY_DN1000_c115_g5 / TRINITY_DN1000_c115_g5_i1):
#     --gene_template    "TRINITY_DN1000_c115_g[GCOUNT]"
#     --trans_template   "TRINITY_DN1000_c115_g[GCOUNT]_i[TCOUNT]"
#
# The script scans the input GFF for existing IDs matching the gene template
# and starts the counter above the highest found, so new IDs never collide
# with existing ones.
#
# Note: version suffix (e.g. .1) is intentionally excluded from templates.
# New genes are always version 1 — versions are incremented only when
# manually editing exon boundaries after the initial merge.

my $default_date_tag  = strftime("%Y%m%d", localtime);
my $gene_template     = "MERGE${default_date_tag}g[GCOUNT:6]000";
my $trans_template    = "MERGE${default_date_tag}t[GCOUNT:6][TCOUNT:3]";

GetOptions(
    "flags=s"              => \$flags_opt,
    "skip_flags=s"         => \$skip_flags_opt,
    "min_tiling=i"         => \$min_tiling,
    "min_cov=f"            => \$min_cov,
    "require_isoseq=s"     => \$require_isoseq,
    "isoseq_min_spanning=i"=> \$isoseq_min_span,
    "fix_partial"          => \$fix_partial,
    "removed=s"            => \$removed_file,
    "dry_run"              => \$dry_run,
    "date=s"               => \$date,
    "gene_template=s"      => \$gene_template,
    "trans_template=s"     => \$trans_template,
) or die "Error in options\n";

my $usage = "usage: $0 [options] isoseq_validated.txt input.gff output.gff\n";
my $merge_file  = shift or die $usage;
my $gff_file    = shift or die $usage;
my $output_file = shift;
die $usage unless defined $output_file || $dry_run;
$output_file //= "/dev/null";  # dry_run: never written, just needs a value

# parse allowed flags (include filter)
my %allowed_flags;
if ($flags_opt eq "all") {
    $allowed_flags{all} = 1;
} else {
    for my $flag_item (split /,/, $flags_opt) { $allowed_flags{$flag_item} = 1; }
}

# parse skip flags (exclude filter) — any row whose flag contains one of these is skipped
my @skip_flag_list;
if ($skip_flags_opt) {
    @skip_flag_list = split /,/, $skip_flags_opt;
}

# ---------------------------------------------------------------------------
# Parse ID templates into components and build generator functions
# ---------------------------------------------------------------------------

# parse_id_template: takes a template string like "CCA3g1[GCOUNT:6]000"
# returns hashref with:
#   prefix      - literal text before [GCOUNT]
#   gcount_len  - zero-pad width for gene counter (0 = no padding)
#   middle      - literal text between [GCOUNT] and [TCOUNT] (if any)
#   tcount_len  - zero-pad width for transcript counter (0 = no padding)
#   suffix      - literal text after last bracket block
#   has_tcount  - 1 if [TCOUNT] present
sub parse_id_template {
    my ($tmpl) = @_;
    my %t = (prefix => "", gcount_len => 6, middle => "",
             tcount_len => 3, suffix => "", has_tcount => 0);

    # extract [GCOUNT:N] or [GCOUNT]
    if ($tmpl =~ /^(.*?)\[GCOUNT(?::(\d+))?\](.*)$/) {
        $t{prefix}     = $1;
        $t{gcount_len} = defined $2 ? $2 + 0 : 0;
        my $rest       = $3;

        # extract [TCOUNT:N] or [TCOUNT] from rest
        if ($rest =~ /^(.*?)\[TCOUNT(?::(\d+))?\](.*)$/) {
            $t{middle}     = $1;
            $t{tcount_len} = defined $2 ? $2 + 0 : 0;
            $t{suffix}     = $3;
            $t{has_tcount} = 1;
        } else {
            $t{suffix} = $rest;
        }
    } else {
        # no brackets found — use template as-is with appended counter
        $t{prefix} = $tmpl;
    }
    return \%t;
}

# build_gene_id: given parsed template and gene counter, return gene ID string
sub build_gene_id {
    my ($t, $gn) = @_;
    my $gc = $t->{gcount_len} > 0
        ? sprintf("%0$t->{gcount_len}d", $gn)
        : $gn;
    return "$t->{prefix}$gc$t->{suffix}";
}

# build_trans_id: given parsed gene template, parsed trans template,
# gene counter and transcript counter, return transcript ID string
sub build_trans_id {
    my ($t, $gn, $tn) = @_;
    my $gc = $t->{gcount_len} > 0
        ? sprintf("%0$t->{gcount_len}d", $gn)
        : $gn;
    if ($t->{has_tcount}) {
        my $tc = $t->{tcount_len} > 0
            ? sprintf("%0$t->{tcount_len}d", $tn)
            : $tn;
        return "$t->{prefix}$gc$t->{middle}$tc$t->{suffix}";
    } else {
        return "$t->{prefix}$gc$t->{suffix}";
    }
}

# parse the templates
my $gene_tmpl  = parse_id_template($gene_template);
my $trans_tmpl = parse_id_template($trans_template);

# ---------------------------------------------------------------------------
# Pass 1: scan GFF to find highest existing gene counter matching template
# so we can start our counter above it and avoid ID collisions
# ---------------------------------------------------------------------------

# build scan pattern from gene template prefix and suffix
my $scan_prefix  = quotemeta($gene_tmpl->{prefix});
my $scan_gclen   = $gene_tmpl->{gcount_len};
my $scan_suffix  = quotemeta($gene_tmpl->{suffix});
my $scan_pattern = $scan_gclen > 0
    ? qr/ID=${scan_prefix}(\d{$scan_gclen})${scan_suffix}/
    : qr/ID=${scan_prefix}(\d+)${scan_suffix}/;

my $max_counter = 0;
open GFF1, $gff_file or die "cant open gff: $gff_file $!\n";
while (my $scan_line = <GFF1>) {
    next unless $scan_line =~ $scan_pattern;
    my $n = $1 + 0;
    $max_counter = $n if $n > $max_counter;
}
close GFF1;

my $gene_counter = $max_counter + 1;
print STDERR "Gene ID template:          $gene_template\n";
print STDERR "Transcript ID template:    $trans_template\n";
print STDERR "Starting gene counter at:  ", build_gene_id($gene_tmpl, $gene_counter), "\n";

# ---------------------------------------------------------------------------
# Load merge table: collect merge groups to process
# ---------------------------------------------------------------------------

# column indices (0-based) for the isoseq_validated.txt
# works for both old (18-col) and new (20-col) format
# we key on genes_in_order (col 2) and flag (col 14 old / col 16 new)
# detect format from header
my ($col_merge_id, $col_genes, $col_flag, $col_tiling, $col_cov,
    $col_isoseq_flag, $col_isoseq_count, $col_isoseq_detail);

open MERGE, $merge_file or die "cant open merge file: $merge_file $!\n";
my $hdr = <MERGE>;
chomp $hdr;
my @hdr_cols = split "\t", $hdr;
for my $i (0 .. $#hdr_cols) {
    $col_merge_id       = $i if $hdr_cols[$i] eq "merge_id";
    $col_genes          = $i if $hdr_cols[$i] eq "genes_in_order";
    $col_flag           = $i if $hdr_cols[$i] eq "flag";
    $col_tiling         = $i if $hdr_cols[$i] eq "max_tiling_hits";
    $col_cov            = $i if $hdr_cols[$i] eq "best_combined_cov_pct";
    $col_isoseq_flag    = $i if $hdr_cols[$i] eq "isoseq_flag";
    $col_isoseq_count   = $i if $hdr_cols[$i] eq "spanning_isoseq_count";
    $col_isoseq_detail  = $i if $hdr_cols[$i] eq "spanning_isoseq_detail";
}
die "Could not find required columns in merge file header\n"
    unless defined $col_genes && defined $col_flag;

my @merge_groups;   # list of arrayrefs of gene IDs to merge
my %gene_in_merge;  # gene_id -> 1 if it is part of any merge group

my $n_total    = 0;
my $n_skipped  = 0;
my $n_accepted = 0;

while (my $mline = <MERGE>) {
    chomp $mline;
    my @f = split "\t", $mline;
    $n_total++;

    my $flag          = $f[$col_flag]        // "";
    my $tiling        = $f[$col_tiling]      // 0;
    my $cov           = $f[$col_cov]         // 0;
    my $isoseq_flag    = defined $col_isoseq_flag   ? ($f[$col_isoseq_flag]   // "") : "";
    my $spanning_count = defined $col_isoseq_count  ? ($f[$col_isoseq_count]  // 0)  : 0;

    # include filter: --flags
    unless (exists $allowed_flags{all}) {
        my $ok = 0;
        for my $af (keys %allowed_flags) {
            $ok = 1 if $flag =~ /\b\Q$af\E\b/;
        }
        unless ($ok) { $n_skipped++; next; }
    }

    # exclude filter: --skip_flags — skip if flag contains any skip term
    my $skip = 0;
    for my $sf (@skip_flag_list) {
        if ($flag =~ /\b\Q$sf\E\b/) { $skip = 1; last; }
    }
    if ($skip) { $n_skipped++; next; }

    # numeric filters
    if ($tiling < $min_tiling) { $n_skipped++; next; }
    if ($cov    < $min_cov)    { $n_skipped++; next; }

    # isoseq filters
    if ($require_isoseq && $isoseq_flag ne $require_isoseq) { $n_skipped++; next; }
    if ($isoseq_min_span > 0 && $spanning_count < $isoseq_min_span) { $n_skipped++; next; }

    my @genes = split ",", $f[$col_genes];
    next unless @genes >= 2;

    # --fix_partial: if PARTIAL_SPAN, trim unsupported terminal genes using
    # the isoseq spanning detail column, then merge the remainder
    if ($fix_partial && $isoseq_flag eq "PARTIAL_SPAN") {
        if (defined $col_isoseq_detail) {
            my $detail = $f[$col_isoseq_detail] // "none";
            # collect which genes are covered by any spanning read
            my %covered_by_iso;
            for my $entry (split /;/, $detail) {
                if ($entry =~ /\(\d+of\d+\):(.+)/) {
                    for my $gid (split /\+/, $1) {
                        $covered_by_iso{$gid} = 1;
                    }
                }
            }
            # trim uncovered terminal genes iteratively
            while (@genes >= 2) {
                last if exists $covered_by_iso{ $genes[0] };
                shift @genes;
            }
            while (@genes >= 2) {
                last if exists $covered_by_iso{ $genes[-1] };
                pop @genes;
            }
            # if trimming left fewer than 2 genes, skip this merge
            if (@genes < 2) {
                my $orig_genes = $f[$col_genes] // "";
                my $mid_val    = defined $col_merge_id ? ($f[$col_merge_id] // "?") : "?";
                print STDERR "PARTIAL_SPAN trim left <2 genes, skipping $mid_val: $orig_genes\n";
                $n_skipped++;
                next;
            }
        }
    }

    my $merge_id_val  = defined $col_merge_id ? ($f[$col_merge_id] // "?") : "?";
    my $new_gene_id_preview = build_gene_id($gene_tmpl, $gene_counter + $n_accepted);
    push @merge_groups, {
        genes       => \@genes,
        merge_id    => $merge_id_val,
        new_gene_id => $new_gene_id_preview,
    };
    for my $g (@genes) { $gene_in_merge{$g} = 1; }
    $n_accepted++;
}
close MERGE;

print STDERR "Merge candidates in file:  $n_total\n";
print STDERR "Skipped by filters:        $n_skipped\n";
print STDERR "Accepted for merging:      $n_accepted\n";
print STDERR "Source genes to replace:   ", scalar keys %gene_in_merge, "\n";

if ($dry_run) {
    print STDERR "\n-- DRY RUN: merge groups that would be processed --\n";
    print STDERR "Output would be written to: $output_file\n";
    print STDERR "Removed genes would go to:  $removed_file\n";
    print STDERR "fix_partial:                ", ($fix_partial ? "yes" : "no"), "\n";
    print STDERR "skip_flags:                 ", ($skip_flags_opt || "none"), "\n\n";
    for my $entry (@merge_groups) {
        print STDERR $entry->{merge_id}, "\t",
                     $entry->{new_gene_id}, "\t",
                     join(",", @{ $entry->{genes} }), "\n";
    }
    exit 0;
}

# open output files
open OUT,     ">$output_file"  or die "cant open output file: $output_file $!\n";
open REMOVED, ">$removed_file" or die "cant open removed file: $removed_file $!\n";
# reconstruct the command line for provenance
my $cmdline = join(" ", $0,
    ($fix_partial       ? "--fix_partial"                      : ()),
    ($skip_flags_opt    ? "--skip_flags $skip_flags_opt"       : ()),
    ($flags_opt ne "all"? "--flags $flags_opt"                 : ()),
    ($min_tiling > 1    ? "--min_tiling $min_tiling"           : ()),
    ($min_cov > 0       ? "--min_cov $min_cov"                 : ()),
    ($require_isoseq    ? "--require_isoseq $require_isoseq"   : ()),
    ($isoseq_min_span>0 ? "--isoseq_min_spanning $isoseq_min_span" : ()),
    ($removed_file ne "removed_genes.gff" ? "--removed $removed_file" : ()),
    $merge_file, $gff_file, $output_file
);

my $provenance_comment = "# Created by: $cmdline";
my $date_comment       = "# Date: $date";

print REMOVED "##gff-version 3\n";
print REMOVED "# Genes removed by merge_split_genes.pl on $date\n";
print REMOVED "# Command: $cmdline\n";
print REMOVED "# Each gene block is preceded by a comment showing its merge group\n";

# open log file named after output gff
(my $log_file = $output_file) =~ s/\.gff$/.log/i;
$log_file .= ".log" unless $log_file =~ /\.log$/;
open LOG, ">$log_file" or die "cant open log file: $log_file $!\n";
print LOG "## Mender - merge_split_genes.pl run log\n";
print LOG "## date:           $date\n";
print LOG "## command:        $cmdline\n";
print LOG "## input_gff:      $gff_file\n";
print LOG "## merge_table:    $merge_file\n";
print LOG "## output_gff:     $output_file\n";
print LOG "## removed_gff:    $removed_file\n";
print LOG "## fix_partial:    ", ($fix_partial ? "yes" : "no"), "\n";
print LOG "## skip_flags:     ", ($skip_flags_opt || "none"), "\n";
print LOG "## flags:          $flags_opt\n";
print LOG "## min_tiling:     $min_tiling\n";
print LOG "## min_cov:        $min_cov\n";
print LOG "##\n";
print LOG "## merge_id\tnew_gene_id\tsource_genes\n";

# ---------------------------------------------------------------------------
# Pass 2: load the full GFF into memory, grouped by gene
# We need:
#   %gene_features  : gene_id -> [GFF lines for the gene feature itself]
#   %gene_children  : gene_id -> {transcript_id -> [GFF lines]}
#   %transcript_children : transcript_id -> [GFF lines for exon/CDS/UTR]
#   @gff_order      : all gene IDs in original order for output
# Lines not belonging to a gene model are stored separately (comments, etc)
# ---------------------------------------------------------------------------
my %gene_line;           # gene_id -> the gene GFF line
my %gene_transcripts;    # gene_id -> [transcript_ids in order of appearance]
my %transcript_line;     # transcript_id -> the mRNA GFF line
my %transcript_features; # transcript_id -> [child feature lines]
my %transcript_parent;   # transcript_id -> gene_id
my @gff_order;           # all lines in original order as [type, key, line]
                         # type: 'header', 'gene', 'passthrough'
my @header_lines;        # ##gff-version etc

my $gff_line_num = 0;

open GFF2, $gff_file or die "cant open gff: $gff_file $!\n";
while (my $line = <GFF2>) {
    chomp $line;
    $gff_line_num++;

    # preserve header/comment lines
    if ($line =~ /^#/) {
        push @header_lines, $line;
        next;
    }

    # skip blank lines
    next unless $line =~ /\S/;

    my @f = split "\t", $line;
    if (@f < 9) {
        print STDERR "WARNING: skipping malformed line $gff_line_num: $line\n";
        next;
    }

    my $type  = $f[2];
    my $attrs = $f[8];

    if ($type eq "gene") {
        my ($id) = $attrs =~ /ID=([^;]+)/;
        unless (defined $id) {
            print STDERR "WARNING: gene feature at line $gff_line_num has no ID, skipping\n";
            next;
        }
        $gene_line{$id} = $line;
        push @gff_order, ['gene', $id, $line];
    }
    elsif ($type eq "mRNA") {
        my ($id)     = $attrs =~ /ID=([^;]+)/;
        my ($parent) = $attrs =~ /Parent=([^;]+)/;
        unless (defined $id && defined $parent) {
            print STDERR "WARNING: mRNA at line $gff_line_num missing ID or Parent, skipping\n";
            next;
        }
        # strict: parent gene must already exist
        unless (exists $gene_line{$parent}) {
            die "ERROR: mRNA '$id' at line $gff_line_num has Parent='$parent' " .
                "which has not been seen yet.\n" .
                "Your GFF is not sorted (parent must appear before child).\n" .
                "Sort with: gt gff3 -sort -tidy -retainids your.gff > sorted.gff\n" .
                "       or: agat_convert_sp_gxf2gxf.pl -g your.gff -o sorted.gff\n";
        }
        $transcript_line{$id}   = $line;
        $transcript_parent{$id} = $parent;
        push @{ $gene_transcripts{$parent} }, $id;
    }
    elsif (defined $attrs && $attrs =~ /Parent=/) {
        # exon, CDS, UTR, or any other child feature
        my ($parent) = $attrs =~ /Parent=([^;]+)/;
        # strict: parent transcript must already exist
        unless (defined $parent && exists $transcript_line{$parent}) {
            # check if parent is a gene (some GFFs have features directly
            # under gene with no mRNA level) — pass through as-is
            if (defined $parent && exists $gene_line{$parent}) {
                # direct child of gene — treat as pass-through attached to gene
                push @{ $transcript_features{"__gene__$parent"} }, $line;
            } else {
                my $parent_str = $parent // "undef";
                die "ERROR: feature '$type' at line $gff_line_num has Parent='$parent_str' " .
                    "which has not been seen yet.\n" .
                    "Your GFF is not sorted (parent must appear before child).\n" .
                    "Sort with: gt gff3 -sort -tidy -retainids your.gff > sorted.gff\n" .
                    "       or: agat_convert_sp_gxf2gxf.pl -g your.gff -o sorted.gff\n";
            }
        } else {
            push @{ $transcript_features{$parent} }, $line;
        }
    }
    else {
        # no Parent attribute and not a gene — pass-through line
        # (repeat regions, sequence-region pragmas, etc)
        push @gff_order, ['passthrough', undef, $line];
    }
}
close GFF2;

print STDERR "GFF lines read: $gff_line_num\n";
print STDERR "Genes loaded:   ", scalar keys %gene_line, "\n";
print STDERR "Transcripts:    ", scalar keys %transcript_line, "\n";

# ---------------------------------------------------------------------------
# Helper: parse GFF attributes into ordered list of [key,value] pairs
# ---------------------------------------------------------------------------
sub parse_attrs {
    my ($attrs) = @_;
    my @pairs;
    for my $kv (split /;/, $attrs) {
        $kv =~ s/^\s+|\s+$//g;
        next unless $kv =~ /=/;
        my ($k, $v) = split /=/, $kv, 2;
        push @pairs, [$k, $v];
    }
    return @pairs;
}

sub attrs_to_string {
    my @pairs = @_;
    return join(";", map { my $p = $_; "$p->[0]=$p->[1]" } @pairs);
}

# ---------------------------------------------------------------------------
# Helper: get genomic coordinates from a GFF line
# ---------------------------------------------------------------------------
sub line_coords {
    my ($line) = @_;
    my @f = split "\t", $line;
    return ($f[3]+0, $f[4]+0);
}

# ---------------------------------------------------------------------------
# Helper: recalculate CDS phases for a sorted list of CDS lines on + strand
# or - strand. Phase = (3 - (cumulative_length % 3)) % 3
# For - strand, CDS features are sorted descending by coordinate.
# ---------------------------------------------------------------------------
sub recalc_phases {
    my ($cds_lines_ref, $strand) = @_;
    my @cds = @$cds_lines_ref;

    # sort by coordinate: ascending for +, descending for -
    if ($strand eq "+") {
        @cds = sort { (split "\t", $a)[3] <=> (split "\t", $b)[3] } @cds;
    } else {
        @cds = sort { (split "\t", $b)[4] <=> (split "\t", $a)[4] } @cds;
    }

    my $cumulative = 0;
    my @result;
    for my $line (@cds) {
        my @f = split "\t", $line;
        my $len   = $f[4] - $f[3] + 1;
        my $phase = (3 - ($cumulative % 3)) % 3;
        $f[7] = $phase;
        push @result, join("\t", @f);
        $cumulative += $len - $phase;
    }
    return @result;
}

# ---------------------------------------------------------------------------
# Build merge operations: for each merge group, construct the new gene
# ---------------------------------------------------------------------------
my %genes_to_remove;  # gene_id -> 1
my @new_gene_blocks;  # each element is a list of GFF lines for one new gene

my $n_merged       = 0;
my $n_transcripts  = 0;

for my $group (@merge_groups) {
    my @genes    = @{ $group->{genes} };
    my $merge_id = $group->{merge_id};

    # verify all genes exist in GFF
    my @missing = grep { !exists $gene_line{$_} } @genes;
    if (@missing) {
        print STDERR "WARNING: genes not found in GFF, skipping merge [",
                     join(",", @genes), "]: missing=", join(",", @missing), "\n";
        next;
    }

    # verify all genes are on the same strand
    my @strands = map { (split "\t", $gene_line{$_})[6] } @genes;
    my %seen_strands = map { $_ => 1 } @strands;
    if (keys %seen_strands > 1) {
        print STDERR "WARNING: mixed-strand chain, skipping merge [$merge_id]: ",
                     join(",", map { "$genes[$_]($strands[$_])" } 0..$#genes), "\n";
        next;
    }

    # get strand and chromosome from first gene
    my @gf    = split "\t", $gene_line{ $genes[0] };
    my $chr   = $gf[0];
    my $strand = $gf[6];

    # new gene coordinates: min start to max end across all genes
    my $new_start = (sort { $a <=> $b }
                     map { (split "\t", $gene_line{$_})[3] } @genes)[0];
    my $new_end   = (sort { $b <=> $a }
                     map { (split "\t", $gene_line{$_})[4] } @genes)[0];

    # assign new gene ID using template
    my $new_gene_id = build_gene_id($gene_tmpl, $gene_counter++);

    # write to log (guard against duplicate writes)
    state %logged;
    unless ($logged{$merge_id}++) {
        print LOG join("\t", $merge_id, $new_gene_id, join(",", @genes)), "\n";
    }

    # build merged_from list
    my $merged_from = join(",", @genes);

    # detect multi-isoform sources: any gene with > 1 transcript
    # (the MULTI_ISOFORM_JOIN flag in the merge table comes from find_split_genes.pl;
    # here we re-check independently and record it as a GFF attribute so the
    # information is preserved in the output annotation)
    my $multi_isoform = 0;
    for my $g (@genes) {
        $multi_isoform = 1 if scalar(@{ $gene_transcripts{$g} // [] }) > 1;
    }
    if ($multi_isoform) {
        my @multi = grep { scalar(@{ $gene_transcripts{$_} // [] }) > 1 } @genes;
        print STDERR "INFO: MULTI_ISOFORM_JOIN [$merge_id]: cross-product transcripts created; ",
                     "verify isoform structure. Multi-isoform source genes: ",
                     join(",", @multi), "\n";
    }

    # new gene GFF line
    my @gene_attr_parts = (
        "ID=$new_gene_id",
        "Name=$new_gene_id",
        "merged_from=$merged_from",
        "merge_source=Mender",
        "merge_id=$merge_id",
        "merge_date=$date",
    );
    push @gene_attr_parts, "multi_isoform_join=1" if $multi_isoform;
    my $new_gene_attrs = join(";", @gene_attr_parts);
    my $new_gene_line = join("\t",
        $chr, "Mender", "gene",
        $new_start, $new_end,
        ".", $strand, ".",
        $new_gene_attrs
    );

    # collect all transcripts from all genes in genomic order
    # for each gene, get its transcripts sorted by start coordinate
    my @all_source_transcripts;  # ordered list of transcript IDs
    for my $gene (@genes) {
        my @tids = @{ $gene_transcripts{$gene} // [] };
        # sort transcripts by their start coordinate
        @tids = sort {
            (split "\t", $transcript_line{$a})[3]
            <=>
            (split "\t", $transcript_line{$b})[3]
        } @tids;
        push @all_source_transcripts, @tids;
    }

    # build cross-product of transcripts across genes
    # group transcripts by their parent gene
    my %tids_by_gene;
    for my $tid (@all_source_transcripts) {
        push @{ $tids_by_gene{ $transcript_parent{$tid} } }, $tid;
    }

    # generate all combinations: one transcript from each gene
    # start with [[]] and extend
    my @combos = ([]);
    for my $gene (@genes) {
        my @tids = @{ $tids_by_gene{$gene} // [] };
        if (!@tids) {
            # gene has no transcripts — skip this gene's contribution
            next;
        }
        my @new_combos;
        for my $combo (@combos) {
            for my $tid (@tids) {
                push @new_combos, [@$combo, $tid];
            }
        }
        @combos = @new_combos;
    }

    # assign transcript numbers sequentially
    my $t_num = 1;
    my @new_transcript_blocks;

    for my $combo (@combos) {
        my @source_tids = @$combo;
        my $new_tid = build_trans_id($trans_tmpl,
                              $gene_counter - 1,  # same numeric as gene
                              $t_num);
        $t_num++;

        # merged_from for transcript
        my $t_merged_from = join(",", @source_tids);

        # new transcript GFF line
        # get coords: min start to max end across all source transcripts
        my $t_start = (sort { $a <=> $b }
                       map { (split "\t", $transcript_line{$_})[3] }
                       @source_tids)[0];
        my $t_end   = (sort { $b <=> $a }
                       map { (split "\t", $transcript_line{$_})[4] }
                       @source_tids)[0];

        my $t_attrs = join(";",
            "ID=$new_tid",
            "Parent=$new_gene_id",
            "Name=$new_tid",
            "merged_from=$t_merged_from",
            "merge_source=Mender",
            "merge_id=$merge_id"
        );
        my $t_line = join("\t",
            $chr, "Mender", "mRNA",
            $t_start, $t_end,
            ".", $strand, ".",
            $t_attrs
        );

        # collect sub-features by type from all source transcripts
        my %by_type;  # type -> [lines]
        for my $tid (@source_tids) {
            for my $feat_line (@{ $transcript_features{$tid} // [] }) {
                my @ff = split "\t", $feat_line;
                push @{ $by_type{ $ff[2] } }, $feat_line;
            }
        }

        # sort each feature type by coordinate (ascending for +, desc for -)
        my @new_feat_lines;
        my %type_counter;

        # determine which source transcript is most upstream and downstream
        # for UTR filtering — sort source_tids by transcript start coordinate
        my @sorted_source_tids = sort {
            (split "\t", $transcript_line{$a})[3]
            <=>
            (split "\t", $transcript_line{$b})[3]
        } @source_tids;
        my $most_upstream_tid   = $sorted_source_tids[0];
        my $most_downstream_tid = $sorted_source_tids[-1];

        # process in a defined order
        for my $ftype ("exon", "CDS", "five_prime_UTR", "three_prime_UTR") {
            my @feats = @{ $by_type{$ftype} // [] };
            next unless @feats;

            # UTR filtering: only keep UTRs from the appropriate terminal transcript
            # five_prime_UTR: only from the most upstream source transcript
            # three_prime_UTR: only from the most downstream source transcript
            # This prevents internal UTRs from split gene fragments appearing
            # mid-transcript after merging
            if ($ftype eq "five_prime_UTR") {
                @feats = grep {
                    my $fl = $_;
                    my ($parent) = (split "\t", $fl)[8] =~ /Parent=([^;]+)/;
                    defined $parent && $parent eq $most_upstream_tid;
                } @feats;
                next unless @feats;
            }
            elsif ($ftype eq "three_prime_UTR") {
                @feats = grep {
                    my $fl = $_;
                    my ($parent) = (split "\t", $fl)[8] =~ /Parent=([^;]+)/;
                    defined $parent && $parent eq $most_downstream_tid;
                } @feats;
                next unless @feats;
            }

            # sort by start coordinate ascending
            @feats = sort {
                (split "\t", $a)[3] <=> (split "\t", $b)[3]
            } @feats;

            # recalculate CDS phases after joining
            if ($ftype eq "CDS") {
                @feats = recalc_phases(\@feats, $strand);
            }

            $type_counter{$ftype} = 0;
            for my $fl (@feats) {
                $type_counter{$ftype}++;
                my @ff = split "\t", $fl;
                my $new_feat_id = sprintf("%s.%s.%d",
                    $new_tid, $ftype, $type_counter{$ftype});
                $ff[1] = "Mender";
                $ff[8] = "ID=$new_feat_id;Parent=$new_tid";
                push @new_feat_lines, join("\t", @ff);
            }
        }

        # handle any other feature types not in the list above
        for my $ftype (sort keys %by_type) {
            next if grep { $_ eq $ftype } qw(exon CDS five_prime_UTR three_prime_UTR);
            my @feats = sort {
                (split "\t", $a)[3] <=> (split "\t", $b)[3]
            } @{ $by_type{$ftype} };
            $type_counter{$ftype} = 0;
            for my $fl (@feats) {
                $type_counter{$ftype}++;
                my @ff = split "\t", $fl;
                my $new_feat_id = sprintf("%s.%s.%d",
                    $new_tid, $ftype, $type_counter{$ftype});
                $ff[1] = "Mender";
                $ff[8] = "ID=$new_feat_id;Parent=$new_tid";
                push @new_feat_lines, join("\t", @ff);
            }
        }

        push @new_transcript_blocks, $t_line, @new_feat_lines;
        $n_transcripts++;
    }

    # mark source genes for removal
    for my $g (@genes) { $genes_to_remove{$g} = 1; }

    push @new_gene_blocks, [$new_gene_line, @new_transcript_blocks];
    $n_merged++;
}

# ---------------------------------------------------------------------------
# Pass 3: write output GFF
# Print header, then all genes in original order, replacing merged ones
# with new gene blocks. New gene blocks are inserted at the position of
# the first source gene in the merge group.
# ---------------------------------------------------------------------------

# build a set of which positions in @gff_order are the first gene of a merge
my %first_gene_of_merge;  # gene_id -> new_gene_block index
my %new_block_by_first;   # first_gene_id -> block
{
    my $block_idx = 0;
    for my $block (@new_gene_blocks) {
        # extract merged_from from the gene line (first element of block)
        my ($merged_from) = $block->[0] =~ /merged_from=([^;]+)/;
        next unless defined $merged_from;
        my @source_genes = split /,/, $merged_from;
        # find which source gene appears first in @gff_order
        my %rank;
        for my $i (0 .. $#gff_order) {
            my $entry = $gff_order[$i];
            $rank{ $entry->[1] } = $i if defined $entry->[1];
        }
        my $first = (sort { ($rank{$a}//9999) <=> ($rank{$b}//9999) }
                     @source_genes)[0];
        $new_block_by_first{$first} = $block;
        $block_idx++;
    }
}

# build lookup: gene_id -> merge_id and partner genes for removed.gff comments
my %gene_merge_comment;  # gene_id -> comment string
for my $group (@merge_groups) {
    # find the merge_id for this group from the merge table
    # we stored groups in order so find by gene set
    my $partners = join(",", @{ $group->{genes} });
    my $mid      = $group->{merge_id};
    for my $g (@{ $group->{genes} }) {
        $gene_merge_comment{$g} = "$mid: $partners";
    }
}

# print header to output file
# if the first header line is a ##gff-version directive, insert provenance
# at index 1 (immediately after ##gff-version), otherwise prepend it
my @out_header = @header_lines;
my $provenance1 = "# Created by: $cmdline";
my $provenance2 = "# Date: $date";
if (@out_header && $out_header[0] =~ /^##gff-version/) {
    splice @out_header, 1, 0, $provenance1, $provenance2;
} else {
    unshift @out_header, $provenance2, $provenance1;
}
for my $hline (@out_header) { print OUT "$hline\n"; }

# print all features in original order, replacing merged genes with new blocks
for my $entry (@gff_order) {
    my ($etype, $key, $line) = @$entry;

    if ($etype eq "passthrough") {
        print OUT "$line\n";
        next;
    }

    # gene entry
    my $gene_id = $key;

    if (exists $genes_to_remove{$gene_id}) {
        # write this gene and all its children to removed.gff
        my $comment = $gene_merge_comment{$gene_id} // "unknown";
        print REMOVED "# removed gene: $gene_id | merged_with: $comment\n";
        print REMOVED $gene_line{$gene_id}, "\n";
        for my $feat (@{ $transcript_features{"__gene__$gene_id"} // [] }) {
            print REMOVED "$feat\n";
        }
        for my $tid (@{ $gene_transcripts{$gene_id} // [] }) {
            print REMOVED $transcript_line{$tid}, "\n";
            for my $feat_line (@{ $transcript_features{$tid} // [] }) { print REMOVED "$feat_line\n"; }
        }

        # if this is the first gene of a merge group, emit new block to output
        if (exists $new_block_by_first{$gene_id}) {
            for my $new_line (@{ $new_block_by_first{$gene_id} }) { print OUT "$new_line\n"; }
        }
        next;
    }

    # print original gene line to output
    print OUT $gene_line{$gene_id}, "\n";

    # print any direct-child features (not via mRNA)
    for my $feat (@{ $transcript_features{"__gene__$gene_id"} // [] }) {
        print OUT "$feat\n";
    }

    # print its transcripts and their children
    for my $tid (@{ $gene_transcripts{$gene_id} // [] }) {
        print OUT $transcript_line{$tid}, "\n";
        for my $feat_line (@{ $transcript_features{$tid} // [] }) { print OUT "$feat_line\n"; }
    }
}

close OUT;
close REMOVED;

print LOG "##\n";
print LOG "## total_merged:   $n_merged\n";
print LOG "## total_removed:  ", scalar keys %genes_to_remove, "\n";
print LOG "## total_transcripts: $n_transcripts\n";
close LOG;

print STDERR "Output GFF written to:     $output_file\n";
print STDERR "Removed genes written to:  $removed_file\n";
print STDERR "Log written to:            $log_file\n";
print STDERR "New merged genes written:  $n_merged\n";
print STDERR "New transcripts written:   $n_transcripts\n";
print STDERR "Source genes removed:      ", scalar keys %genes_to_remove, "\n";
