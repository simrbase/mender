#!/usr/bin/perl
use strict;
use warnings;

## validate_merge_translation.pl
##
## Step 8 of the Mender pipeline. Post-merge translation quality validator.
## Validates merged gene models from merge_split_genes.pl (and, when run,
## rebuild_with_miniprot.pl) by:
##
##   A  Parse merged GFF — collect new gene IDs, source gene lists, GFF blocks,
##      and the miniprot result recorded on each gene
##   C  Translate merged proteins (gffread -y); pick representative per gene
##   E  Build or verify diamond database from ref_fa
##   F  Diamond blastp merged proteins vs reference proteome; compute coverage metrics
##   I  Load optional merge table for flag annotations
##   K  Assign PASS / FAIL / REVIEW per merge
##   L  Write outputs (TSV report, pass/fail/review GFF3s, pass proteins FASTA)
##
## (Steps B, D, G, H and J -- the MSA junction scoring and the source-protein,
## reference-sequence and SwissProt loading that fed it -- were removed on
## 2026-09-28: in a negative control of 100 deliberately false merges every one
## scored GOOD_MSA, the same range as real merges. Protein coverage and the
## miniprot rebuild separate them; see README.)
##
## Usage:
##   perl validate_merge_translation.pl \
##       --merged_gff merged.gff3 \
##       --genome_fa  genome.fa \
##       --ref_fa     reference_proteome.fa \
##       --out        validation_run1 \
##       [options]
##
## Required:
##   --merged_gff FILE    Output of merge_split_genes.pl (or rebuild_with_miniprot.pl)
##   --genome_fa FILE     Genome sequence (for gffread -y)
##   --ref_fa FILE        Reference proteome FASTA (same as diamond step)
##   --out PREFIX         Prefix for all output files
##
## Options:
##   --diamond_db FILE    Pre-built diamond DB from ref_fa; built if absent
##   --merge_table FILE   TSV from find_split_genes.pl / validate_merge_with_isoseq.pl
##                        Used to annotate source merge flags in the report
##   --threads N          Diamond threads (default: 4)
##   --evalue E           Diamond evalue cutoff (default: 1e-10)
##   --scratch DIR        Directory for intermediate/scratch files (default: same as --out)
##   --min_merged_cov F   Min merged protein coverage by best ref hit (default: 0.80)
##   --min_ref_cov F      Min ref protein coverage by merged protein (default: 0.50)
##   --gffread_bin PATH   Path to gffread (default: gffread in PATH)
##   --diamond_bin PATH   Path to diamond (default: diamond in PATH)
##
## Outputs:
##   <out>/report.tsv               Per-merge validation results (TSV, 17 columns)
##   <out>/pass.gff3                Merged genes that PASS all checks
##   <out>/fail.gff3                Merged genes that FAIL (clear evidence of bad merge)
##   <out>/review.gff3              Merged genes that need REVIEW (borderline)
##   <out>/pass_proteins.fa         Translated proteins for PASS merges
##
## Note: <out> is used as a directory name, not a filename prefix.

use Getopt::Long;
use File::Path qw(make_path);
use List::Util qw(max min sum);

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
my $merged_gff     = undef;
my $genome_fa      = undef;
my $ref_fa         = undef;
my $diamond_db     = undef;
my $merge_table    = undef;
my $out            = undef;
my $threads        = 4;
my $evalue         = "1e-10";
my $min_merged_cov = 0.80;
my $min_ref_cov    = 0.50;
my $gffread_bin    = "gffread";
my $diamond_bin    = "diamond";
my $scratch        = undef;   # directory for intermediate/scratch files

GetOptions(
    "merged_gff=s"      => \$merged_gff,
    "genome_fa=s"       => \$genome_fa,
    "ref_fa=s"          => \$ref_fa,
    "diamond_db=s"      => \$diamond_db,
    "merge_table=s"     => \$merge_table,
    "out=s"             => \$out,
    "threads=i"         => \$threads,
    "evalue=s"          => \$evalue,
    "min_merged_cov=f"  => \$min_merged_cov,
    "min_ref_cov=f"     => \$min_ref_cov,
    "gffread_bin=s"     => \$gffread_bin,
    "diamond_bin=s"     => \$diamond_bin,
    "scratch=s"         => \$scratch,
) or die "Run: perl $0 --help for usage\n";

die "ERROR: --merged_gff is required\n" unless defined $merged_gff;
die "ERROR: --genome_fa is required\n"  unless defined $genome_fa;
die "ERROR: --ref_fa is required\n"     unless defined $ref_fa;
die "ERROR: --out is required\n"        unless defined $out;

for my $pair ([$merged_gff, "merged_gff"], [$genome_fa,  "genome_fa"],  [$ref_fa,   "ref_fa"]) {
    die "ERROR: --$pair->[1] file not found: $pair->[0]\n" unless -f $pair->[0];
}
if (defined $diamond_db) {
    die "ERROR: --diamond_db file not found: $diamond_db\n"
        unless -f $diamond_db || -f "${diamond_db}.dmnd";
}

# $out is the run directory; $scratch is for intermediate files
my $out_dir     = $out;
my $scratch_dir = defined($scratch) ? $scratch : $out;

# Ensure both directories exist
make_path($out_dir)     unless -d $out_dir;
make_path($scratch_dir) unless -d $scratch_dir;

print "=" x 60, "\n";
print "VALIDATE MERGE TRANSLATION\n";
print "=" x 60, "\n";
printf "merged_gff:     %s\n", $merged_gff;
printf "genome_fa:      %s\n", $genome_fa;
printf "ref_fa:         %s\n", $ref_fa;
printf "out_dir:        %s\n", $out_dir;
printf "scratch_dir:    %s\n", $scratch_dir;
printf "min_merged_cov: %.2f\n", $min_merged_cov;
printf "min_ref_cov:    %.2f\n", $min_ref_cov;
printf "threads:        %d\n",   $threads;
print "-" x 60, "\n";

# ---------------------------------------------------------------------------
# Helper: run a shell command and die on failure
# ---------------------------------------------------------------------------
sub run_cmd {
    my ($cmd, $desc) = @_;
    print "CMD [$desc]: $cmd\n";
    my $ret = system($cmd);
    die "ERROR: command failed (exit $ret): $desc\n" if $ret != 0;
}

# ---------------------------------------------------------------------------
# Helper: read a FASTA file into a hash {id => sequence}
# Returns (\%seqs, \@order)
# ---------------------------------------------------------------------------
sub read_fasta {
    my ($file) = @_;
    my (%seqs, @order);
    my $id;
    open my $fh, "<", $file or die "Cannot open FASTA $file: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>(\S+)/) {
            $id = $1;
            unless (exists $seqs{$id}) {
                push @order, $id;
                $seqs{$id} = "";
            }
        } elsif (defined $id) {
            $seqs{$id} .= $line;
        }
    }
    close $fh;
    return (\%seqs, \@order);
}

# ---------------------------------------------------------------------------
# Helper: write a FASTA file from a list of (id, seq) pairs
# ---------------------------------------------------------------------------
sub write_fasta {
    my ($file, @entries) = @_;  # @entries = ([id, seq], ...)
    open my $fh, ">", $file or die "Cannot write FASTA $file: $!\n";
    for my $e (@entries) {
        print $fh ">$e->[0]\n$e->[1]\n";
    }
    close $fh;
}

# ---------------------------------------------------------------------------
# Helper: parse GFF3 attribute string → hash
# ---------------------------------------------------------------------------
sub parse_attrs {
    my ($attr_str) = @_;
    my %attrs;
    for my $kv (split /;/, $attr_str) {
        if ($kv =~ /^(\w+)=(.*)$/) {
            $attrs{$1} = $2;
        }
    }
    return %attrs;
}

# ---------------------------------------------------------------------------
# Helper: map unaligned aa position (1-based) to alignment column index
# ---------------------------------------------------------------------------
sub unaligned_to_col {
    my ($seq, $target) = @_;
    my $count = 0;
    for my $i (0 .. length($seq) - 1) {
        my $c = substr($seq, $i, 1);
        unless ($c eq '-' || $c eq '.') {
            $count++;
            return $i if $count == $target;
        }
    }
    return length($seq) - 1;  # fallback: last column
}

# ---------------------------------------------------------------------------
# STEP A: Parse the merged GFF
#
# Collect:
#   %merge_info{new_gene_id}  = { merge_id, source_genes[], chrom, strand }
#   %mrgd_tid_to_gene{tid}    = new_gene_id
#   %gene_blocks{new_gene_id} = [gff lines]
#   @gff_gene_order           = [new_gene_ids in file order]
#   @gff_header_lines         = pragma/comment lines before first gene
# ---------------------------------------------------------------------------
print "\nSTEP A: Parsing merged GFF...\n";

my %merge_info;          # new_gene_id -> { merge_id, source_genes, chrom, strand }
my %mrgd_tid_to_gene;    # merged transcript ID -> new gene ID
my %gene_blocks;         # new_gene_id -> [array of GFF lines]
my @gff_gene_order;      # gene IDs in file order
my @gff_header_lines;    # ##pragma and comment lines before first gene block

{
    my $cur_gene_id = undef;
    open my $fh, "<", $merged_gff or die "Cannot open $merged_gff: $!\n";
    while (my $line = <$fh>) {
        chomp $line;

        # header lines (pragma, comments, blank) before first gene block
        if (!defined $cur_gene_id && ($line =~ /^#/ || $line =~ /^\s*$/)) {
            push @gff_header_lines, $line;
            next;
        }
        next if $line =~ /^#/;   # skip comments inside body too
        next if $line =~ /^\s*$/;

        my @f = split /\t/, $line;
        next unless @f >= 9;

        my $feat = $f[2];
        my %attrs = parse_attrs($f[8]);

        if ($feat eq "gene") {
            my $gid = $attrs{ID} // next;
            $cur_gene_id = $gid;
            push @gff_gene_order, $gid;
            push @{ $gene_blocks{$gid} }, $line;

            # only process Mender-merged genes
            if (defined $attrs{merged_from}) {
                my @src = split /,/, $attrs{merged_from};
                $merge_info{$gid} = {
                    merge_id     => $attrs{merge_id} // ".",
                    source_genes => \@src,
                    chrom        => $f[0],
                    strand       => $f[6],
                    frame_break  => $attrs{frame_break} // "",   # merge_split_genes.pl: junctions out of frame
                    miniprot     => $attrs{miniprot}    // "",   # rebuild_with_miniprot.pl: SUPPORTED | NOT_SUPPORTED:<reason>
                };
            }
        } elsif (defined $cur_gene_id) {
            push @{ $gene_blocks{$cur_gene_id} }, $line;
            if ($feat eq "mRNA") {
                my $tid = $attrs{ID} // next;
                $mrgd_tid_to_gene{$tid} = $cur_gene_id;
            }
        }
    }
    close $fh;
}

my $n_merged = scalar keys %merge_info;
printf "  Merged genes found:       %d\n", $n_merged;
printf "  Total genes in merged GFF: %d\n", scalar @gff_gene_order;
die "ERROR: no merged genes (merged_from attribute) found in $merged_gff\n" unless $n_merged > 0;

# ---------------------------------------------------------------------------
# STEP C: Translate merged proteins with gffread
# ---------------------------------------------------------------------------
print "\nSTEP C: Translating merged proteins (gffread)...\n";

my $merged_prot_fa = "$scratch_dir/merged_proteins.fa";
my $gffread_log    = "$scratch_dir/gffread.log";

# Truncate gffread log; steps C and D both append to it.
open(my $grl, ">", $gffread_log) or warn "Cannot create $gffread_log: $!\n"; close $grl;

run_cmd("$gffread_bin $merged_gff -g $genome_fa -y $merged_prot_fa >> $gffread_log 2>&1",
        "gffread translate merged GFF");

# Parse merged proteins — pick representative per gene (longest without internal stop)
my ($mrgd_prot_ref, $mrgd_prot_order) = read_fasta($merged_prot_fa);

my %merged_prot;           # new_gene_id -> protein sequence (representative)
my %has_internal_stop;     # new_gene_id -> 0/1
my %tid_internal_stop;     # merged transcript id -> 0/1 (every transcript, not just the representative)
my %internal_stop_pos;     # new_gene_id -> comma-sep positions or "."
my %merged_prot_len;       # new_gene_id -> aa length (after stripping terminal *)

for my $tid (@$mrgd_prot_order) {
    my $gid = $mrgd_tid_to_gene{$tid};
    next unless defined $gid && exists $merge_info{$gid};

    my $raw = $mrgd_prot_ref->{$tid};
    # gffread emits either '*' or '.' for a stop codon depending on version/
    # invocation (see the tr/./X/ cleanup below) -- match both here too.
    (my $seq = $raw) =~ s/[*.]$//;  # strip terminal stop

    # detect internal stops
    my @ipos;
    while ($seq =~ /[*.]/g) { push @ipos, pos($seq) }
    my $has_stop = @ipos ? 1 : 0;
    $tid_internal_stop{$tid} = $has_stop;

    # pick representative: prefer no internal stops, then longest
    if (!exists $merged_prot{$gid}) {
        $merged_prot{$gid}        = $seq;
        $has_internal_stop{$gid}  = $has_stop;
        $internal_stop_pos{$gid}  = @ipos ? join(",", @ipos) : ".";
        $merged_prot_len{$gid}    = length($seq);
    } else {
        # prefer stop-free over stop-containing; break ties by length
        my $cur_stop = $has_internal_stop{$gid};
        my $cur_len  = $merged_prot_len{$gid};
        if (($cur_stop && !$has_stop) ||
            (!$cur_stop && !$has_stop && length($seq) > $cur_len) ||
            ($cur_stop  &&  $has_stop && length($seq) > $cur_len)) {
            $merged_prot{$gid}       = $seq;
            $has_internal_stop{$gid} = $has_stop;
            $internal_stop_pos{$gid} = @ipos ? join(",", @ipos) : ".";
            $merged_prot_len{$gid}   = length($seq);
        }
    }
}

my $n_prot_found = scalar keys %merged_prot;
printf "  Merged proteins extracted: %d / %d\n", $n_prot_found, $n_merged;
printf "  With internal stops:       %d\n", scalar grep { $has_internal_stop{$_} } keys %has_internal_stop;

# Write a clean per-gene FASTA for diamond:
#   - one sequence per merged gene (representative)
#   - internal stops (*) removed
#   - unknown AA (.) replaced with X  (gffread emits . for N-containing codons;
#     diamond rejects . but accepts X)
#   - header is the gene ID so diamond hits come back keyed by gene ID directly
my $diamond_query_fa = "$scratch_dir/merged_proteins_clean.fa";
{
    open my $fh, ">", $diamond_query_fa
        or die "Cannot write $diamond_query_fa: $!\n";
    for my $gid (sort keys %merged_prot) {
        (my $clean = $merged_prot{$gid}) =~ s/\*//g;
        $clean =~ tr/./X/;
        print $fh ">$gid\n$clean\n";
    }
    close $fh;
}

# ---------------------------------------------------------------------------
# STEP E: Build or verify diamond database
# ---------------------------------------------------------------------------
print "\nSTEP E: Diamond database...\n";

unless (defined $diamond_db) {
    $diamond_db = "$scratch_dir/ref.dmnd";
}
unless (-f $diamond_db || -f "${diamond_db}.dmnd") {
    run_cmd("$diamond_bin makedb --in $ref_fa --db $diamond_db --quiet",
            "Build diamond database from ref_fa");
}
print "  Diamond db: $diamond_db\n";

# ---------------------------------------------------------------------------
# STEP F: Diamond blastp — merged proteins vs reference
# ---------------------------------------------------------------------------
print "\nSTEP F: Diamond blastp (merged proteins vs reference)...\n";

my $diamond_out = "$scratch_dir/merged_vs_ref.tsv";
run_cmd(
    "$diamond_bin blastp " .
    "--db $diamond_db " .
    "--query $diamond_query_fa " .
    "--out $diamond_out " .
    "--outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore scovhsp slen qcovhsp qlen " .
    "--evalue $evalue " .
    "--max-target-seqs 5 " .
    "--threads $threads " .
    "--quiet",
    "Diamond blastp merged proteins"
);

# Parse diamond output — query IDs are now gene IDs (from the clean FASTA)
# outfmt 6 fields (0-indexed):
#   0:qseqid  1:sseqid  2:pident  3:length  4:qstart  5:qend
#   6:sstart  7:send    8:evalue  9:bitscore 10:scovhsp 11:slen
#   12:qcovhsp 13:qlen
my %gene_diamond_hits;  # new_gene_id -> [{...}, ...]  sorted by bitscore desc
my %gene_cov_metrics;   # new_gene_id -> { merged_cov_by_ref, ... }

{
    my %raw_hits;  # gene_id -> [hits]
    open my $dh, "<", $diamond_out or die "Cannot open $diamond_out: $!\n";
    while (my $line = <$dh>) {
        chomp $line;
        my @f = split /\t/, $line;
        next unless @f >= 14;
        my $gid = $f[0];
        push @{ $raw_hits{$gid} }, {
            sseqid   => $f[1],
            pident   => $f[2]+0,
            qstart   => $f[4]+0,
            qend     => $f[5]+0,
            sstart   => $f[6]+0,
            send     => $f[7]+0,
            evalue   => $f[8],
            bitscore => $f[9]+0,
            scovhsp  => $f[10]+0,
            slen     => $f[11]+0,
            qcovhsp  => $f[12]+0,
            qlen     => $f[13]+0,
        };
    }
    close $dh;

    # Sort hits by bitscore, deduplicate by sseqid, keep top 5
    for my $gid (keys %raw_hits) {
        my @sorted = sort { $b->{bitscore} <=> $a->{bitscore} } @{ $raw_hits{$gid} };
        my %seen;
        my @top;
        for my $h (@sorted) {
            next if $seen{$h->{sseqid}}++;
            push @top, $h;
            last if @top >= 5;
        }
        $gene_diamond_hits{$gid} = \@top;
    }
}

for my $gid (keys %merge_info) {
    my @top_hits = @{ $gene_diamond_hits{$gid} // [] };

    if (@top_hits) {
        my $best = $top_hits[0];
        my $qlen = $best->{qlen} || $merged_prot_len{$gid} || 1;
        $gene_cov_metrics{$gid} = {
            merged_cov_by_ref  => ($best->{qend} - $best->{qstart} + 1) / $qlen,
            ref_cov_by_merged  => $best->{scovhsp} / 100,
            best_ref_hit       => $best->{sseqid},
            best_ref_pident    => $best->{pident},
            ref_len            => $best->{slen},
            n_ref_hits         => scalar @top_hits,
        };
    } else {
        $gene_cov_metrics{$gid} = {
            merged_cov_by_ref  => 0,
            ref_cov_by_merged  => 0,
            best_ref_hit       => ".",
            best_ref_pident    => 0,
            ref_len            => 0,
            n_ref_hits         => 0,
        };
    }
}

my $n_hits = grep { scalar(@{ $gene_diamond_hits{$_} }) > 0 } keys %gene_diamond_hits;
printf "  Merged genes with diamond hits: %d / %d\n", $n_hits, $n_merged;

# ---------------------------------------------------------------------------
# STEP I: Load optional merge table for flag annotations
# ---------------------------------------------------------------------------
my %merge_flags;  # merge_id -> flag string

if (defined $merge_table && -f $merge_table) {
    print "\nLoading merge table for flag annotations...\n";
    open my $mfh, "<", $merge_table or die "Cannot open $merge_table: $!\n";
    my $hdr = <$mfh>; chomp $hdr;
    my @cols = split /\t/, $hdr;
    my %ci;
    $ci{$cols[$_]} = $_ for 0..$#cols;
    my $c_mid  = $ci{merge_id} // $ci{merge_id};
    my $c_flag = $ci{flag};
    while (my $line = <$mfh>) {
        chomp $line;
        my @f = split /\t/, $line;
        next unless defined $c_mid && defined $c_flag;
        my $mid  = $f[$c_mid]  // next;
        my $flag = $f[$c_flag] // ".";
        $merge_flags{$mid} = $flag;
    }
    close $mfh;
    printf "  Merge flags loaded: %d entries\n", scalar keys %merge_flags;
}

# ---------------------------------------------------------------------------
# STEP K: PASS / FAIL / REVIEW assignment
# ---------------------------------------------------------------------------
print "\nSTEP K: Assigning PASS / FAIL / REVIEW...\n";

my %results;  # new_gene_id -> { overall_result, fail_reasons }

for my $gid (keys %merge_info) {
    my $transl_flag  = $has_internal_stop{$gid} ? "FRAMESHIFT_DETECTED" : "OK";
    my $has_hit      = @{ $gene_diamond_hits{$gid} // [] } ? 1 : 0;
    my $cov          = $gene_cov_metrics{$gid};
    my $mrgd_cov     = $cov->{merged_cov_by_ref};
    my $ref_cov      = $cov->{ref_cov_by_merged};
    my @fail_reasons;

    my $result;

    if ($transl_flag eq "FRAMESHIFT_DETECTED") {
        push @fail_reasons, "INTERNAL_STOP";
    }
    # merge_split_genes.pl found a junction where the downstream fragment is out of
    # frame: its translation is wrong whether or not a stop happens to appear
    if ($merge_info{$gid}{frame_break}) {
        push @fail_reasons, "FRAME_BREAK";
    }
    # rebuild_with_miniprot.pl could not align the tiling reference protein across the
    # merged locus as one gene (it covers only some fragments, has a frameshift or
    # in-frame stop, or too little of the protein): not supported as one gene
    if ($merge_info{$gid}{miniprot} =~ /^NOT_SUPPORTED/) {
        push @fail_reasons, "NO_PROTEIN_SPAN";
    }
    if (!$has_hit) {
        push @fail_reasons, "NO_HIT";
    } elsif ($mrgd_cov < $min_merged_cov) {
        # One reference protein must cover most of the merged protein. This is the
        # check that separates true from false merges; the MSA junction score does
        # not (100 deliberately false merges of two distinct adjacent genes: all
        # GOOD_MSA, junction scores 0.59-1.00; merged_cov_by_ref median 0.57 vs 0.98
        # for real merges; >= 0.80 kept 89/93 real, 12/100 false). A GOOD_MSA no
        # longer overrides this.
        push @fail_reasons, "LOW_MERGED_COV";
    }

    if (@fail_reasons) {
        $result = "FAIL";
    } elsif ($transl_flag eq "OK"
          && $mrgd_cov >= $min_merged_cov
          && $ref_cov  >= $min_ref_cov) {
        $result = "PASS";
    } else {
        $result = "REVIEW";
    }

    $results{$gid} = {
        overall_result   => $result,
        translation_flag => $transl_flag,
        fail_reasons     => \@fail_reasons,
    };
}

my %counts = (PASS => 0, FAIL => 0, REVIEW => 0);
$counts{ $results{$_}{overall_result} }++ for keys %results;
printf "  PASS:   %d\n  FAIL:   %d\n  REVIEW: %d\n",
    $counts{PASS}, $counts{FAIL}, $counts{REVIEW};

# ---------------------------------------------------------------------------
# STEP L: Write outputs
# ---------------------------------------------------------------------------
print "\nSTEP L: Writing outputs...\n";

# --- L1: TSV report ---
my $report_file = "$out_dir/report.tsv";
open my $rep, ">", $report_file or die "Cannot write $report_file: $!\n";
print $rep join("\t",
    "merge_id",
    "new_gene_id",
    "source_genes",
    "merged_protein_len",
    "has_internal_stop",
    "internal_stop_pos",
    "merged_cov_by_ref",
    "ref_cov_by_merged",
    "best_ref_hit",
    "best_ref_pident",
    "ref_len",
    "n_ref_hits",
    "miniprot",
    "translation_flag",
    "source_flags",
    "fail_reasons",
    "overall_result",
), "\n";

for my $gid (sort keys %merge_info) {
    my $info  = $merge_info{$gid};
    my $cov   = $gene_cov_metrics{$gid};
    my $res   = $results{$gid};
    my $fail_str = @{ $res->{fail_reasons} } ? join("|", @{ $res->{fail_reasons} }) : ".";
    my $src_flag = $merge_flags{ $info->{merge_id} } // ".";

    print $rep join("\t",
        $info->{merge_id},
        $gid,
        join(",", @{ $info->{source_genes} }),
        $merged_prot_len{$gid}           // 0,
        $has_internal_stop{$gid}         // 0,
        $internal_stop_pos{$gid}         // ".",
        sprintf("%.4f", $cov->{merged_cov_by_ref}),
        sprintf("%.4f", $cov->{ref_cov_by_merged}),
        $cov->{best_ref_hit},
        sprintf("%.2f", $cov->{best_ref_pident}),
        $cov->{ref_len},
        $cov->{n_ref_hits},
        ($info->{miniprot} ne "" ? $info->{miniprot} : "not_run"),
        $res->{translation_flag},
        $src_flag,
        $fail_str,
        $res->{overall_result},
    ), "\n";
}
close $rep;
print "  Report:      $report_file\n";

# --- L2: Split GFF3 outputs ---
my $pass_gff   = "$out_dir/pass.gff3";
my $fail_gff   = "$out_dir/fail.gff3";
my $review_gff = "$out_dir/review.gff3";

open my $pfh, ">", $pass_gff   or die "Cannot write $pass_gff: $!\n";
open my $ffh, ">", $fail_gff   or die "Cannot write $fail_gff: $!\n";
open my $rfh, ">", $review_gff or die "Cannot write $review_gff: $!\n";

# Write GFF3 version header to all output files
my @gff_hdr = grep { /^##gff-version/ || /^##genome-build/ } @gff_header_lines;
push @gff_hdr, "##gff-version 3" unless @gff_hdr;
for my $fh ($pfh, $ffh, $rfh) {
    print $fh "$_\n" for @gff_hdr;
}

my %result_fh = (PASS => $pfh, FAIL => $ffh, REVIEW => $rfh);

for my $gid (@gff_gene_order) {
    next unless exists $merge_info{$gid};   # skip non-merged genes
    my $res = $results{$gid}{overall_result} // next;
    my $fh  = $result_fh{$res} // next;
    my @block = @{ $gene_blocks{$gid} // [] };

    # A merged gene is judged by its best transcript, but a multi-isoform merge
    # carries N x M cross-product transcripts: those with an internal stop are not
    # written to the PASS file (with their exons/CDS/UTRs); the gene records them.
    my @dropped;
    if ($res eq "PASS") {
        my %drop = map { $_ => 1 } grep { $tid_internal_stop{$_} && ($mrgd_tid_to_gene{$_} // "") eq $gid } keys %tid_internal_stop;
        if (%drop) {
            @block = grep {
                my @f = split "\t", $_;
                my %a = @f >= 9 ? parse_attrs($f[8]) : ();
                !(($f[2] // "") eq "mRNA" && $drop{$a{ID} // ""}) && !$drop{$a{Parent} // ""};
            } @block;
            @dropped = sort keys %drop;
        }
    }

    # Add translation validation flags to the gene feature line (first line of block)
    if (@block) {
        $block[0] .= ";transl_result=$res";
        $block[0] .= ";dropped_transcripts_internal_stop=" . join(",", @dropped) if @dropped;
    }
    print $fh "$_\n" for @block;
}
close $pfh; close $ffh; close $rfh;
printf "  Pass GFF:   %s  (%d genes)\n", $pass_gff,   $counts{PASS};
printf "  Fail GFF:   %s  (%d genes)\n", $fail_gff,   $counts{FAIL};
printf "  Review GFF: %s  (%d genes)\n", $review_gff, $counts{REVIEW};

# --- L3: Pass proteins FASTA ---
my $pass_prot_fa = "$out_dir/pass_proteins.fa";
open my $ppfh, ">", $pass_prot_fa or die "Cannot write $pass_prot_fa: $!\n";
for my $gid (sort keys %results) {
    next unless $results{$gid}{overall_result} eq "PASS";
    next unless defined $merged_prot{$gid} && $merged_prot{$gid} ne "";
    my $merge_id = $merge_info{$gid}{merge_id};
    my $src_str  = join(",", @{ $merge_info{$gid}{source_genes} });
    print $ppfh ">$gid merge_id=$merge_id merged_from=$src_str\n$merged_prot{$gid}\n";
}
close $ppfh;
printf "  Pass proteins: %s\n", $pass_prot_fa;

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
print "\n", "=" x 60, "\n";
print "VALIDATE MERGE TRANSLATION COMPLETE\n";
printf "  PASS:   %d\n  FAIL:   %d\n  REVIEW: %d\n",
    $counts{PASS}, $counts{FAIL}, $counts{REVIEW};
print "=" x 60, "\n";
