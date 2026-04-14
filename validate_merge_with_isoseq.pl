#!/usr/bin/perl
use strict;
use warnings;

## validate_with_isoseq.pl
##
## Validates split gene merge candidates by finding IsoSeq reads that span
## two or more genes in a merge group. A spanning IsoSeq read is independent
## molecular evidence that the split genes are transcribed as a single unit.
##
## Steps to prepare input:
##   grep -P "\tgene\t" genes.gff > models.gene.gff
##   grep -P "\tmRNA\t" CCA3C_isoseq.agat.mRNA.gff > isoseq.mrna.gff
##   bedtools intersect -a isoseq.mrna.gff -b models.gene.gff -wa -wb | \
##       perl -p -e 's/.+\sID=([^;]+).+\sID=([^;]+).*$/$1\t$2/' > overlaps
##
## Usage: perl validate_with_isoseq.pl merge_candidates.txt overlaps > isoseq_validated.txt
##
## Output columns (all columns from merge table plus):
##
##   spanning_isoseq_count : total number of IsoSeq reads that overlap 2+
##                           genes in this merge group
##
##   spanning_isoseq_detail: one entry per spanning read, sorted by how many
##                           group genes it covers (most first), then by read ID.
##                           Format per entry:
##                             read_id(N_covered of M_total):gene1+gene2+...
##                           Example:
##                             PB.5678.2(3of3):CCA3g000013+CCA3g000014+CCA3g000015
##                           Multiple entries separated by semicolons.
##                           'none' if no spanning reads found.
##
##                           The (NofM) notation tells you how much of the merge
##                           group each read covers:
##                             (2of2) = read spans both genes in a pair — full support
##                             (2of3) = read spans 2 of 3 genes — partial support
##                             (3of3) = read spans all 3 genes — full support
##                           Genes are listed in genomic order (same order as
##                           genes_in_order column in merge table).
##
##   isoseq_flag           : structural signal from IsoSeq spanning pattern.
##                           Note: absence of spanning reads is NOT flagged —
##                           IsoSeq came from specific developmental stages so
##                           zero spanning reads may simply reflect expression
##                           pattern, not gene structure.
##
##                           Flag values:
##                             FULL_SPAN    : at least one read spans all genes
##                                           in the group — strong confirmation
##                             PARTIAL_SPAN : spanning reads exist BUT none of
##                                           them reach one or both terminal
##                                           genes — structural signal that the
##                                           terminal gene may not belong in the
##                                           merge. This is meaningful even with
##                                           stage-specific IsoSeq because reads
##                                           long enough to span genes 1+2 should
##                                           also cross gene 3 if they are truly
##                                           one transcriptional unit.
##                             none         : no spanning reads found (not a flag,
##                                           just absence of evidence)
##
##                           Example of PARTIAL_SPAN:
##                             3-gene merge A,B,C
##                             Reads span A+B consistently, none reach C
##                             -> PARTIAL_SPAN: gene C may not belong
##                             This was observed for JMJD1C where all 7 reads
##                             spanned genes 1+2 but none crossed into gene 3.
##
##
## Note: validate_merge_with_isoseq.pl is optional. This script accepts either
##   merge_candidates.txt or isoseq_validated.txt as input — columns are
##   detected by header name. If IsoSeq columns are absent, options that
##   depend on them (--fix_partial, --require_isoseq, --isoseq_min_spanning)
##   are silently ignored and filtering falls back to protein evidence only.
##
##   Without IsoSeq, recommended usage:
##
##   perl merge_split_genes.pl \
##       --flags STRONG \
##       --skip_flags SKIPPED_GENE,LOW_COV \
##       merge_candidates.txt \
##       input.gff \
##       output.gff
##

my $merge_file    = shift or die "usage: $0 merge_candidates.txt overlaps\n";
my $overlaps_file = shift or die "usage: $0 merge_candidates.txt overlaps\n";

# ---------------------------------------------------------------------------
# Load overlaps: isoseq_read -> {gene -> 1} and gene -> {isoseq_read -> 1}
# Each line is one IsoSeq-gene overlap from bedtools intersect
# ---------------------------------------------------------------------------
my %iso_to_genes;  # iso_read_id -> {gene_id -> 1}
my %gene_to_isos;  # gene_id     -> {iso_read_id -> 1}

open OVERLAPS, $overlaps_file or die "cant open overlaps file: $overlaps_file $!\n";
while (my $line = <OVERLAPS>){
    chomp $line;
    my ($isoseq, $gene_id) = split "\t", $line;
    next unless defined $isoseq && defined $gene_id;
    $iso_to_genes{$isoseq}{$gene_id} = 1;
    $gene_to_isos{$gene_id}{$isoseq} = 1;
}
close OVERLAPS;

print STDERR "IsoSeq reads in overlaps: ",  scalar keys %iso_to_genes, "\n";
print STDERR "Genes with IsoSeq overlap: ", scalar keys %gene_to_isos, "\n";

# ---------------------------------------------------------------------------
# Process merge table, adding spanning read evidence to each row
# ---------------------------------------------------------------------------
open MERGE, $merge_file or die "cant open merge file: $merge_file $!\n";

my $header = <MERGE>;
chomp $header;
print $header, "\t", "spanning_isoseq_count", "\t", "spanning_isoseq_detail",
               "\t", "isoseq_flag", "\n";

my $total         = 0;
my $with_spanning = 0;
my $full_span     = 0;
my $partial_span  = 0;

while (my $line = <MERGE>){
    chomp $line;
    my @fields = split "\t", $line;

    # gene IDs are in column 3 (0-indexed: col 2), comma separated
    # these are in genomic order from the merge table
    my @genes   = split ",", $fields[2];
    my $n_genes = scalar @genes;

    # collect all IsoSeq reads that overlap any gene in this group
    # for each read, track which group genes it overlaps
    my %candidate_isos;  # iso_read_id -> {gene_id -> 1}
    for my $gene (@genes){
        if (exists $gene_to_isos{$gene}){
            for my $iso (keys %{ $gene_to_isos{$gene} }){
                $candidate_isos{$iso}{$gene} = 1;
            }
        }
    }

    # a spanning read overlaps 2+ genes in this group
    my @spanning_details;
    for my $iso (sort keys %candidate_isos){
        # get the group genes this read covers, in genomic order
        my @covered = grep { exists $candidate_isos{$iso}{$_} } @genes;
        my $n_covered = scalar @covered;
        next unless $n_covered >= 2;

        # format: read_id(N_covered of N_total):gene1+gene2+...
        my $detail = sprintf "%s(%dof%d):%s",
            $iso,
            $n_covered,
            $n_genes,
            join("+", @covered);

        push @spanning_details, [$n_covered, $iso, $detail, [@covered]];
    }

    # sort by genes covered descending, then alphabetically by read ID
    my @sorted = sort { $b->[0] <=> $a->[0] || $a->[1] cmp $b->[1] }
                 @spanning_details;

    my $spanning_count  = scalar @sorted;
    my $spanning_detail = $spanning_count > 0
        ? join(";", map { $_->[2] } @sorted)
        : "none";

    # ---------------------------------------------------------------------------
    # Determine isoseq_flag
    # FULL_SPAN: at least one read covers all genes in the group
    # PARTIAL_SPAN: reads exist but no read reaches one or both terminal genes
    #   — check if any terminal gene is absent from ALL spanning reads
    #   — this is a structural signal regardless of expression stage
    # ---------------------------------------------------------------------------
    my $isoseq_flag = "NO_SPANNERS";

    if ($spanning_count > 0) {
        # check if any read spans all genes
        my $has_full = grep { $_->[0] == $n_genes } @sorted;

        if ($has_full) {
            $isoseq_flag = "FULL_SPAN";
        } else {
            # check which terminal genes are covered by at least one spanning read
            my $first_gene = $genes[0];
            my $last_gene  = $genes[-1];

            my $first_covered = grep {
                grep { $_ eq $first_gene } @{ $_->[3] }
            } @sorted;

            my $last_covered = grep {
                grep { $_ eq $last_gene } @{ $_->[3] }
            } @sorted;

            # PARTIAL_SPAN if a terminal gene is not reached by any spanning read
            # Only meaningful for chains with 3+ genes — for 2-gene pairs
            # any spanning read is by definition covering both terminals
            if ($n_genes >= 3 && (!$first_covered || !$last_covered)) {
                $isoseq_flag = "PARTIAL_SPAN";
                $partial_span++;
            } else {
                # spanning reads exist but none reach all genes and
                # both terminals are covered — just partial internal coverage
                $isoseq_flag = "FULL_SPAN" if $n_genes == 2;
            }
        }
        $full_span++ if $isoseq_flag eq "FULL_SPAN";
    }

    print join("\t", @fields, $spanning_count, $spanning_detail, $isoseq_flag), "\n";

    $total++;
    $with_spanning++ if $spanning_count > 0;
}

close MERGE;

print STDERR "Merge candidates processed: $total\n";
print STDERR "With spanning IsoSeq reads: $with_spanning\n";
print STDERR "  FULL_SPAN:                $full_span\n";
print STDERR "  PARTIAL_SPAN:             $partial_span\n";
print STDERR "Without spanning IsoSeq:    ", $total - $with_spanning, "\n";
