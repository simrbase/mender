#!/usr/bin/perl
use strict;
use warnings;

## validate_merge_translation.pl
##
## Step 8 of the Mender pipeline. Post-merge translation quality validator.
## Validates merged gene models from merge_split_genes.pl by:
##
##   A  Parse merged GFF — collect new gene IDs, source gene lists, GFF blocks
##   B  Parse original GFF — get representative transcript and CDS length per source gene
##   C  Translate merged proteins (gffread -y); pick representative per gene
##   D  Translate source gene proteins (gffread -y)
##   E  Build or verify diamond database from ref_fa
##   F  Diamond blastp merged proteins vs reference proteome; compute coverage metrics
##   G  Load reference proteome sequences into memory
##   H  SwissProt diamond blast (if --swissprot_fa); load hit sequences for MSA
##   I  Load optional merge table for flag annotations
##   J  MSA junction scoring — align merged + source + ref (SwissProt + ref_fa) proteins; score each junction
##   K  Assign PASS / FAIL / REVIEW per merge
##   L  Write outputs (TSV report, pass/fail/review GFF3s, pass proteins FASTA)
##
## Usage:
##   perl validate_merge_translation.pl \
##       --merged_gff merged.gff3 \
##       --orig_gff   original.gff3 \
##       --genome_fa  genome.fa \
##       --ref_fa     reference_proteome.fa \
##       --out        validation_run1 \
##       [options]
##
## Required:
##   --merged_gff FILE    Output of merge_split_genes.pl
##   --orig_gff FILE      Pre-merge annotation (for source gene proteins + CDS lengths)
##   --genome_fa FILE     Genome sequence (for gffread -y)
##   --ref_fa FILE        Reference proteome FASTA (same as diamond step)
##   --out PREFIX         Prefix for all output files
##
## Options:
##   --diamond_db FILE    Pre-built diamond DB from ref_fa; built if absent
##   --merge_table FILE   TSV from find_split_genes.pl / validate_merge_with_isoseq.pl
##                        Used to annotate source merge flags in the report
##   --threads N          Diamond/MAFFT threads (default: 4)
##   --evalue E           Diamond evalue cutoff (default: 1e-10)
##   --aligner STR        MSA tool: mafft_fast (default) | mafft_auto | kalign3
##                        kalign3 is ~8x faster and recommended for routine runs
##   --kalign_bin PATH    Path to kalign binary (default: kalign in PATH)
##   --no_msa             Skip MSA; report translation + diamond only
##   --min_junction F     Min per-junction MSA score to PASS (default: 0.5)
##   --min_merged_cov F   Min merged protein coverage by best ref hit (default: 0.60)
##   --min_ref_cov F      Min ref protein coverage by merged protein (default: 0.50)
##   --max_msa_refs N     Max ref hits to include per MSA (default: 3)
##   --keep_msa           Write per-merge MSA files to <out>_msa/
##   --w_conservation F   Weight for conservation sub-score (default: 0.3)
##   --w_continuity F     Weight for ref continuity sub-score (default: 0.3)
##   --w_gap F            Weight for gap pattern sub-score (default: 0.4)
##   --min_msa_refs N     Min refs in MSA for full-confidence flag (default: 2)
##   --gffread_bin PATH   Path to gffread (default: gffread in PATH)
##   --diamond_bin PATH   Path to diamond (default: diamond in PATH)
##   --mafft_bin PATH     Path to mafft (default: mafft in PATH)
##
## SwissProt MSA reference (recommended):
##   --swissprot_fa FILE  UniProt/SwissProt FASTA for MSA reference sequences.
##                        Manually curated proteins give better junction scoring
##                        than the species reference proteome (no partial seqs,
##                        no annotation errors, no '.' characters).
##                        Falls back to --ref_fa sequences if absent.
##   --swissprot_db FILE  Pre-built diamond DB from swissprot_fa; built if absent.
##
## To download SwissProt:
##   wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
##   gunzip uniprot_sprot.fasta.gz
##   diamond makedb --in uniprot_sprot.fasta --db uniprot_sprot.dmnd
##
## Outputs:
##   <out>_report.tsv               Per-merge validation results
##   <out>_pass.gff3                Merged genes that PASS all checks
##   <out>_fail.gff3                Merged genes that FAIL (clear evidence of bad merge)
##   <out>_review.gff3              Merged genes that need REVIEW (borderline)
##   <out>_pass_proteins.fa         Translated proteins for PASS merges
##   <out>_msa/                     Per-merge alignment files (if --keep_msa)

use Getopt::Long;
use POSIX qw(floor);
use File::Path qw(make_path);
use List::Util qw(max min sum);

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
my $merged_gff     = undef;
my $orig_gff       = undef;
my $genome_fa      = undef;
my $ref_fa         = undef;
my $diamond_db     = undef;
my $merge_table    = undef;
my $out            = undef;
my $threads        = 4;
my $evalue         = "1e-10";
my $no_msa         = 0;
my $min_junction   = 0.5;
my $min_merged_cov = 0.60;
my $min_ref_cov    = 0.50;
my $max_msa_refs   = 3;
my $keep_msa       = 0;
my $w_conservation = 0.3;
my $w_continuity   = 0.3;
my $w_gap          = 0.4;
my $junction_window = 5;
my $min_msa_refs   = 2;   # warn + flag GOOD_MSA_LOW_REF when fewer refs in MSA
my $gffread_bin    = "gffread";
my $diamond_bin    = "diamond";
my $mafft_bin      = "mafft";
my $kalign_bin     = "kalign";
my $aligner        = "mafft_fast";  # mafft_fast | mafft_auto | kalign3
my $swissprot_fa   = undef;
my $swissprot_db   = undef;
my $scratch        = undef;   # directory for intermediate/scratch files

GetOptions(
    "merged_gff=s"      => \$merged_gff,
    "orig_gff=s"        => \$orig_gff,
    "genome_fa=s"       => \$genome_fa,
    "ref_fa=s"          => \$ref_fa,
    "diamond_db=s"      => \$diamond_db,
    "merge_table=s"     => \$merge_table,
    "out=s"             => \$out,
    "threads=i"         => \$threads,
    "evalue=s"          => \$evalue,
    "no_msa"            => \$no_msa,
    "aligner=s"         => \$aligner,
    "min_junction=f"    => \$min_junction,
    "min_merged_cov=f"  => \$min_merged_cov,
    "min_ref_cov=f"     => \$min_ref_cov,
    "max_msa_refs=i"    => \$max_msa_refs,
    "keep_msa"          => \$keep_msa,
    "w_conservation=f"  => \$w_conservation,
    "w_continuity=f"    => \$w_continuity,
    "w_gap=f"           => \$w_gap,
    "min_msa_refs=i"    => \$min_msa_refs,
    "gffread_bin=s"     => \$gffread_bin,
    "diamond_bin=s"     => \$diamond_bin,
    "mafft_bin=s"       => \$mafft_bin,
    "kalign_bin=s"      => \$kalign_bin,
    "swissprot_fa=s"    => \$swissprot_fa,
    "swissprot_db=s"    => \$swissprot_db,
    "scratch=s"         => \$scratch,
) or die "Run: perl $0 --help for usage\n";

die "ERROR: --aligner must be one of: mafft_fast, mafft_auto, kalign3 (got: $aligner)\n"
    unless $aligner =~ /^(mafft_fast|mafft_auto|kalign3)$/;

die "ERROR: --merged_gff is required\n" unless defined $merged_gff;
die "ERROR: --orig_gff is required\n"   unless defined $orig_gff;
die "ERROR: --genome_fa is required\n"  unless defined $genome_fa;
die "ERROR: --ref_fa is required\n"     unless defined $ref_fa;
die "ERROR: --out is required\n"        unless defined $out;

for my $pair ([$merged_gff, "merged_gff"], [$orig_gff, "orig_gff"],
              [$genome_fa,  "genome_fa"],  [$ref_fa,   "ref_fa"]) {
    die "ERROR: --$pair->[1] file not found: $pair->[0]\n" unless -f $pair->[0];
}
if (defined $diamond_db) {
    die "ERROR: --diamond_db file not found: $diamond_db\n"
        unless -f $diamond_db || -f "${diamond_db}.dmnd";
}
if (defined $swissprot_fa && !-f $swissprot_fa) {
    print "WARNING: --swissprot_fa file not found: $swissprot_fa\n";
    print "         MSA will fall back to --ref_fa sequences\n";
    $swissprot_fa = undef;
}

my $wsum = $w_conservation + $w_continuity + $w_gap;
die sprintf("ERROR: --w_conservation + --w_continuity + --w_gap must sum to 1.0 (got %.3f)\n", $wsum)
    unless abs($wsum - 1.0) < 0.001;

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
printf "orig_gff:       %s\n", $orig_gff;
printf "genome_fa:      %s\n", $genome_fa;
printf "ref_fa:         %s\n", $ref_fa;
printf "swissprot_fa:   %s\n", ($swissprot_fa // "(not set — MSA will use ref_fa sequences)");
printf "out_dir:        %s\n", $out_dir;
printf "scratch_dir:    %s\n", $scratch_dir;
printf "no_msa:         %s\n", ($no_msa ? "yes" : "no");
printf "aligner:        %s\n", ($no_msa ? "n/a" : $aligner);
printf "min_junction:   %.2f\n", $min_junction;
printf "min_merged_cov: %.2f\n", $min_merged_cov;
printf "min_ref_cov:    %.2f\n", $min_ref_cov;
printf "w_conservation: %.2f  w_continuity: %.2f  w_gap: %.2f\n",
    $w_conservation, $w_continuity, $w_gap;
printf "min_msa_refs:   %d\n", $min_msa_refs;
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
# STEP B: Parse the original GFF
#
# For each gene, determine the representative transcript (longest total CDS)
# and compute its total CDS bp length.
#
# Collect:
#   %orig_tid_to_gene{tid}   = gene_id
#   %orig_gene_cds_len{gid}  = total CDS bp of representative transcript
#   %orig_gene_rep_tid{gid}  = transcript ID of representative transcript
# ---------------------------------------------------------------------------
print "\nSTEP B: Parsing original GFF for source gene CDS lengths...\n";

my %orig_tid_to_gene;
my %orig_gene_cds_len;    # gene_id -> CDS bp of representative transcript
my %orig_gene_rep_tid;    # gene_id -> rep transcript ID

{
    # First pass: build gene -> transcripts map and transcript -> gene map
    my %gene_transcripts;    # gene_id -> [tid, ...]
    my %tid_cds_len;         # tid -> cumulative CDS bp

    open my $fh, "<", $orig_gff or die "Cannot open $orig_gff: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^#/ || $line =~ /^\s*$/;
        my @f = split /\t/, $line;
        next unless @f >= 9;
        my $feat = $f[2];
        my %attrs = parse_attrs($f[8]);

        if ($feat eq "mRNA") {
            my $tid = $attrs{ID}     // next;
            my $gid = $attrs{Parent} // next;
            $orig_tid_to_gene{$tid} = $gid;
            push @{ $gene_transcripts{$gid} }, $tid;
            $tid_cds_len{$tid} //= 0;
        } elsif ($feat eq "CDS") {
            my $pid = $attrs{Parent} // next;
            # Parent of CDS can be a transcript ID
            my $cds_bp = $f[4] - $f[3] + 1;
            $tid_cds_len{$pid} += $cds_bp;
        }
    }
    close $fh;

    # For each gene, pick representative transcript (longest CDS)
    for my $gid (keys %gene_transcripts) {
        my @tids = @{ $gene_transcripts{$gid} };
        my $best_tid = $tids[0];
        my $best_len = $tid_cds_len{$tids[0]} // 0;
        for my $tid (@tids) {
            my $len = $tid_cds_len{$tid} // 0;
            if ($len > $best_len) {
                $best_len = $len;
                $best_tid = $tid;
            }
        }
        $orig_gene_cds_len{$gid} = $best_len;
        $orig_gene_rep_tid{$gid} = $best_tid;
    }
}

# Build set of all source gene IDs we need
my %need_source_gene;
for my $gid (keys %merge_info) {
    $need_source_gene{$_} = 1 for @{ $merge_info{$gid}{source_genes} };
}
my $n_src_found = grep { defined $orig_gene_cds_len{$_} } keys %need_source_gene;
printf "  Source genes needed:  %d\n", scalar keys %need_source_gene;
printf "  Found in orig GFF:    %d\n", $n_src_found;

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
my %internal_stop_pos;     # new_gene_id -> comma-sep positions or "."
my %merged_prot_len;       # new_gene_id -> aa length (after stripping terminal *)

for my $tid (@$mrgd_prot_order) {
    my $gid = $mrgd_tid_to_gene{$tid};
    next unless defined $gid && exists $merge_info{$gid};

    my $raw = $mrgd_prot_ref->{$tid};
    (my $seq = $raw) =~ s/\*$//;  # strip terminal stop

    # detect internal stops
    my @ipos;
    while ($seq =~ /\*/g) { push @ipos, pos($seq) }
    my $has_stop = @ipos ? 1 : 0;

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
# STEP D: Translate source gene proteins with gffread
# ---------------------------------------------------------------------------
print "\nSTEP D: Translating source proteins (gffread)...\n";

my $source_prot_fa = "$scratch_dir/source_proteins.fa";
run_cmd("$gffread_bin $orig_gff -g $genome_fa -y $source_prot_fa >> $gffread_log 2>&1",
        "gffread translate original GFF");

my ($src_prot_ref, $src_prot_order) = read_fasta($source_prot_fa);

# Map source gene IDs to their representative protein sequence
my %source_prot;  # gene_id -> protein sequence
for my $tid (@$src_prot_order) {
    my $gid = $orig_tid_to_gene{$tid} // next;
    next unless $need_source_gene{$gid};
    my $rep_tid = $orig_gene_rep_tid{$gid} // next;
    next unless $tid eq $rep_tid;
    (my $seq = $src_prot_ref->{$tid}) =~ s/\*$//;
    $source_prot{$gid} = $seq;
}
printf "  Source proteins extracted: %d / %d\n",
    scalar keys %source_prot, scalar keys %need_source_gene;

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
# STEP G: Load reference proteome sequences (for MSA)
# ---------------------------------------------------------------------------
print "\nSTEP G: Loading reference proteome sequences...\n";

my %ref_prot;  # hit_id -> sequence
{
    my ($rp, undef) = read_fasta($ref_fa);
    %ref_prot = %$rp;
}
printf "  Reference sequences loaded: %d\n", scalar keys %ref_prot;

# ---------------------------------------------------------------------------
# STEP H: SwissProt diamond blast + load sequences (for MSA only)
#
# SwissProt (manually curated) gives better junction scoring than the species
# reference proteome: no partial sequences, no annotation errors, no '.'
# characters from N-containing codons.
# Falls back to ref_fa sequences in the MAFFT loop when not set.
# ---------------------------------------------------------------------------
my %sprot_hits;  # new_gene_id -> [{sseqid, bitscore, ...}]
my %sprot_prot;  # hit_id -> sequence
my $use_sprot = 0;

if (defined $swissprot_fa) {
    print "\nSTEP H: SwissProt diamond blast (for MSA reference sequences)...\n";

    # Build/verify SwissProt diamond DB
    unless (defined $swissprot_db) {
        # fall back: build in scratch dir (but normally passed via --swissprot_db from mender_db)
        $swissprot_db = "$scratch_dir/swissprot.dmnd";
    }
    unless (-f $swissprot_db || -f "${swissprot_db}.dmnd") {
        run_cmd("$diamond_bin makedb --in $swissprot_fa --db $swissprot_db --quiet",
                "Build diamond database from swissprot_fa");
    }
    print "  SwissProt db: $swissprot_db\n";

    my $sprot_out = "$scratch_dir/merged_vs_swissprot.tsv";
    run_cmd(
        "$diamond_bin blastp " .
        "--db $swissprot_db " .
        "--query $diamond_query_fa " .
        "--out $sprot_out " .
        "--outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore scovhsp slen qcovhsp qlen " .
        "--evalue $evalue " .
        "--max-target-seqs 5 " .
        "--threads $threads " .
        "--quiet",
        "Diamond blastp merged proteins vs SwissProt"
    );

    # Parse SwissProt hits
    {
        my %sprot_raw;
        open my $sdh, "<", $sprot_out or die "Cannot open $sprot_out: $!\n";
        while (my $line = <$sdh>) {
            chomp $line;
            my @f = split /\t/, $line;
            next unless @f >= 14;
            push @{ $sprot_raw{$f[0]} }, {
                sseqid   => $f[1],
                pident   => $f[2]+0,
                bitscore => $f[9]+0,
                scovhsp  => $f[10]+0,
                slen     => $f[11]+0,
            };
        }
        close $sdh;

        for my $gid (keys %sprot_raw) {
            my @sorted = sort { $b->{bitscore} <=> $a->{bitscore} } @{ $sprot_raw{$gid} };
            my %seen;
            my @top;
            for my $h (@sorted) {
                next if $seen{$h->{sseqid}}++;
                push @top, $h;
                last if @top >= $max_msa_refs;
            }
            $sprot_hits{$gid} = \@top;
        }
    }

    # Collect all SwissProt hit IDs we need, then load just those sequences
    my %need_sprot_ids;
    for my $gid (keys %sprot_hits) {
        $need_sprot_ids{$_->{sseqid}} = 1 for @{ $sprot_hits{$gid} };
    }

    if (%need_sprot_ids) {
        open my $sfh, "<", $swissprot_fa or die "Cannot open $swissprot_fa: $!\n";
        my $sid;
        while (my $line = <$sfh>) {
            chomp $line;
            if ($line =~ /^>(\S+)/) {
                $sid = $1;
                $sid = undef unless $need_sprot_ids{$sid};
            } elsif (defined $sid) {
                $sprot_prot{$sid} .= $line;
            }
        }
        close $sfh;
    }

    my $n_sprot_hits = grep { @{ $sprot_hits{$_} } > 0 } keys %sprot_hits;
    printf "  Merged genes with SwissProt hits: %d / %d\n", $n_sprot_hits, $n_merged;
    printf "  SwissProt sequences loaded: %d\n", scalar keys %sprot_prot;
    $use_sprot = 1;
} else {
    print "\nSTEP H: Skipped (no --swissprot_fa — MSA will use ref_fa sequences)\n";
}

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
# STEP J: MAFFT junction scoring
# ---------------------------------------------------------------------------
my $msa_dir    = "$scratch_dir/msa";
my $kalign_log = "$scratch_dir/kalign.log";

make_path($msa_dir) if $keep_msa && !$no_msa;

# Truncate kalign log at run start; each merge call appends to it.
if (!$no_msa && $aligner eq "kalign3") {
    open(my $fh, ">", $kalign_log) or warn "Cannot create $kalign_log: $!\n";
    close $fh;
}

printf "\nSTEP J: MSA junction scoring (aligner: %s)...\n", $aligner unless $no_msa;
printf "  kalign log: %s\n", $kalign_log if !$no_msa && $aligner eq "kalign3";
print  "\nSTEP J: Skipping MSA (--no_msa)\n" if $no_msa;

my %junction_scores;  # new_gene_id -> [score1, score2, ...]
my %msa_flag;         # new_gene_id -> flag string

my $pid = $$;  # for temp file naming

my $mafft_done = 0;
my @merge_ids_sorted = sort keys %merge_info;
my $mafft_total = scalar @merge_ids_sorted;

for my $gid (@merge_ids_sorted) {
    $mafft_done++;
    if ($mafft_done == 1 || $mafft_done % 50 == 0 || $mafft_done == $mafft_total) {
        printf "  MSA: %d / %d (%.0f%%)\n",
            $mafft_done, $mafft_total, 100 * $mafft_done / $mafft_total;
    }
    my @src_genes   = @{ $merge_info{$gid}{source_genes} };
    my $merge_id    = $merge_info{$gid}{merge_id};
    my @top_hits    = @{ $gene_diamond_hits{$gid} };  # ref hits (for NO_HIT check)
    my @msa_hits    = $use_sprot                       # hits to use for MSA ref seqs
                    ? @{ $sprot_hits{$gid} // [] }
                    : @top_hits;
    my $msa_prot_hr = $use_sprot ? \%sprot_prot : \%ref_prot;
    my $n_junctions = scalar(@src_genes) - 1;

    # --- No diamond hit: skip MSA entirely ---
    if (!@top_hits) {
        $msa_flag{$gid}      = "NO_HIT";
        $junction_scores{$gid} = [];
        next;
    }

    # --- --no_msa flag ---
    if ($no_msa) {
        $msa_flag{$gid}      = "SKIPPED";
        $junction_scores{$gid} = [];
        next;
    }

    # --- Build input FASTA for MAFFT ---
    my @fasta_entries;
    my $merged_label = "MERGED_${gid}";
    $merged_label =~ s/[^A-Za-z0-9._-]/_/g;

    unless (defined $merged_prot{$gid} && $merged_prot{$gid} ne "") {
        $msa_flag{$gid}      = "WEAK_JUNCTION";
        $junction_scores{$gid} = [];
        next;
    }
    # Clean for MSA: remove internal stops, replace N-codon dots with X
    (my $clean_merged = $merged_prot{$gid}) =~ s/\*//g;
    $clean_merged =~ tr/./X/;
    push @fasta_entries, [$merged_label, $clean_merged];

    for my $sg (@src_genes) {
        if (defined $source_prot{$sg} && $source_prot{$sg} ne "") {
            (my $lbl = "SOURCE_${sg}") =~ s/[^A-Za-z0-9._-]/_/g;
            (my $clean_src = $source_prot{$sg}) =~ s/\*//g;
            $clean_src =~ tr/./X/;
            push @fasta_entries, [$lbl, $clean_src];
        } else {
            print STDERR "WARN: no source protein for $sg in merge $merge_id\n";
        }
    }

    # Top N ref sequences — use SwissProt hits if available, else ref_fa hits
    my $n_refs_added = 0;
    for my $hit (@msa_hits) {
        last if $n_refs_added >= $max_msa_refs;
        my $rid = $hit->{sseqid};
        if (defined $msa_prot_hr->{$rid} && $msa_prot_hr->{$rid} ne "") {
            (my $lbl = "REF_${rid}") =~ s/[^A-Za-z0-9._-]/_/g;
            push @fasta_entries, [$lbl, $msa_prot_hr->{$rid}];
            $n_refs_added++;
        }
    }

    # Also add the best ref_fa hit (used for merge detection) when SwissProt is
    # providing the MSA refs — it's the protein the source fragments tiled against,
    # so seeing it in the alignment helps interpret junction quality.
    if ($use_sprot && @top_hits) {
        my $rid = $top_hits[0]{sseqid};
        if (defined $ref_prot{$rid} && $ref_prot{$rid} ne "") {
            (my $lbl = "REFPROT_${rid}") =~ s/[^A-Za-z0-9._-]/_/g;
            push @fasta_entries, [$lbl, $ref_prot{$rid}];
        }
    }

    # Need at least 2 sequences for MAFFT
    if (@fasta_entries < 2) {
        $msa_flag{$gid}      = "WEAK_JUNCTION";
        $junction_scores{$gid} = [];
        next;
    }

    # Write input FASTA
    my $input_fa  = "${msa_dir}/${merge_id}_input.fa";
    my $aligned_fa = "${msa_dir}/${merge_id}_aligned.fa";
    unless ($keep_msa) {
        $input_fa   = "/tmp/vmt_${pid}_${merge_id}_input.fa";
        $aligned_fa = "/tmp/vmt_${pid}_${merge_id}_aligned.fa";
    }
    write_fasta($input_fa, @fasta_entries);

    # Run aligner
    my $merged_len = length($merged_prot{$gid});
    my $aln_cmd;
    if ($aligner eq "kalign3") {
        $aln_cmd = "$kalign_bin -i $input_fa -f fasta -o $aligned_fa >> $kalign_log 2>&1";
    } elsif ($aligner eq "mafft_auto") {
        # --auto with speed fallback for very long proteins
        my $opts = $merged_len > 5000
            ? "--auto --retree 1 --reorder --quiet"
            : "--auto --reorder --quiet";
        $aln_cmd = "$mafft_bin $opts $input_fa > $aligned_fa 2>/dev/null";
    } else {
        # mafft_fast (default): FFT-NS-2, no iterative refinement
        $aln_cmd = "$mafft_bin --retree 2 --maxiterate 0 --reorder --quiet $input_fa > $aligned_fa 2>/dev/null";
    }
    my $aln_ret = system($aln_cmd);
    if ($aln_ret != 0) {
        print STDERR "WARN: aligner ($aligner) failed for $merge_id (exit $aln_ret); flagging WEAK_JUNCTION\n";
        $msa_flag{$gid}      = "WEAK_JUNCTION";
        $junction_scores{$gid} = [];
        unlink $input_fa, $aligned_fa unless $keep_msa;
        next;
    }

    # Parse aligned FASTA
    my ($aln_seqs, $aln_order) = read_fasta($aligned_fa);

    # Find MERGED sequence in alignment (MAFFT may reorder sequences)
    my ($aln_merged_id) = grep { /^MERGED_/ } @$aln_order;
    unless (defined $aln_merged_id && defined $aln_seqs->{$aln_merged_id}) {
        print STDERR "WARN: MERGED sequence missing from MAFFT output for $merge_id\n";
        $msa_flag{$gid}      = "WEAK_JUNCTION";
        $junction_scores{$gid} = [];
        unlink $input_fa, $aligned_fa unless $keep_msa;
        next;
    }
    my $aln_len = length($aln_seqs->{$aln_merged_id});
    my @ref_ids_in_aln = grep { /^REF_/ } @$aln_order;

    # Compute junction positions in unaligned merged protein aa coordinates
    # J_k = floor( sum(CDS_bp for source genes 1..k) / 3 )
    # Then map to alignment column.
    my @junc_aa;     # junction positions in unaligned aa (1-based)
    my $cumulative_bp = 0;
    for my $k (0 .. $n_junctions - 1) {
        my $sg = $src_genes[$k];
        my $cds_bp = $orig_gene_cds_len{$sg} // 0;
        if ($cds_bp == 0) {
            print STDERR "WARN: CDS length unknown for $sg — junction $k position may be inaccurate\n";
        }
        $cumulative_bp += $cds_bp;
        push @junc_aa, floor($cumulative_bp / 3);
    }

    # Score each junction
    my @scores;
    for my $junc_pos (@junc_aa) {
        next if $junc_pos <= 0 || $junc_pos > length($merged_prot{$gid});
        my $junc_col = unaligned_to_col($aln_seqs->{$aln_merged_id}, $junc_pos);
        my ($score, undef, undef, undef) = score_junction(
            $aln_seqs, $aln_order,
            $aln_merged_id, \@ref_ids_in_aln,
            $junc_col, $aln_len,
            $junction_window,
            $w_conservation, $w_continuity, $w_gap
        );
        push @scores, $score;
    }

    $junction_scores{$gid} = \@scores;

    # Warn when the MSA has fewer reference sequences than the minimum; scores
    # are less reliable with sparse reference support.
    if ($n_refs_added < $min_msa_refs) {
        my $warn = "WARN: $merge_id scored with only $n_refs_added reference sequence(s)"
                 . " (min_msa_refs=$min_msa_refs) — junction score is low-confidence\n";
        print STDERR $warn;
        if ($kalign_log) {
            open(my $fh, ">>", $kalign_log) or die "Cannot open $kalign_log: $!";
            print $fh $warn;
            close $fh;
        }
    }

    # Assign msa_flag
    if (!@scores) {
        # No scoreable junctions (e.g., single-gene source with degenerate input)
        $msa_flag{$gid} = "WEAK_JUNCTION";
    } else {
        my $min_score = min(@scores);
        # GOOD_MSA: min junction score meets threshold AND no internal stop.
        # GOOD_MSA_LOW_REF: same, but flagged as low-confidence due to sparse refs.
        if ($min_score >= $min_junction && !$has_internal_stop{$gid}) {
            $msa_flag{$gid} = $n_refs_added < $min_msa_refs ? "GOOD_MSA_LOW_REF" : "GOOD_MSA";
        } else {
            $msa_flag{$gid} = "WEAK_JUNCTION";
        }
    }

    # Cleanup temp files unless keeping
    unless ($keep_msa) {
        unlink $input_fa, $aligned_fa;
    }
}

# ---------------------------------------------------------------------------
# STEP J helper: junction scoring sub
# ---------------------------------------------------------------------------
sub score_junction {
    my ($aln_seqs, $aln_order, $merged_id, $ref_ids, $junc_col, $aln_len,
        $window, $wc, $wt, $wg) = @_;

    my $col_start = max(0, $junc_col - $window);
    my $col_end   = min($aln_len - 1, $junc_col + $window);
    my $n_cols    = $col_end - $col_start + 1;

    my @ref_ids = @$ref_ids;

    # Conservation uses merged + refs only; source fragments are excluded because
    # they are identical to the corresponding halves of the merged protein by
    # construction and would artificially inflate the consensus score.
    my @cons_ids = ($merged_id, @ref_ids);

    # --- Conservation score ---
    my ($con_sum, $con_cols) = (0, 0);
    for my $col ($col_start .. $col_end) {
        my @chars = map { substr($aln_seqs->{$_}, $col, 1) } @cons_ids;
        my @non_gap = grep { $_ ne '-' && $_ ne '.' } @chars;
        next unless @non_gap;
        my %cnt;
        $cnt{$_}++ for @non_gap;
        my $max_cnt = max(values %cnt);
        $con_sum += $max_cnt / scalar(@non_gap);
        $con_cols++;
    }
    my $conservation = $con_cols > 0 ? $con_sum / $con_cols : 0.0;

    # --- Gap pattern score ---
    # Penalise columns where MERGED has gap but REF sequences do not
    my ($merged_gap_cols, $valid_cols) = (0, 0);
    if (@ref_ids) {
        for my $col ($col_start .. $col_end) {
            my $mc = substr($aln_seqs->{$merged_id}, $col, 1);
            my @rc = map { substr($aln_seqs->{$_}, $col, 1) } @ref_ids;
            my @ref_non_gap = grep { $_ ne '-' && $_ ne '.' } @rc;
            $valid_cols++;
            $merged_gap_cols++ if ($mc eq '-' || $mc eq '.') && @ref_non_gap;
        }
    } else {
        $valid_cols = $n_cols;
    }
    my $gap_pattern = $valid_cols > 0 ? 1 - ($merged_gap_cols / $valid_cols) : 0.5;

    # --- Ref continuity score ---
    # Fraction of REF sequences gapless across junc_col ± window (same window
    # as conservation and gap pattern — previously used ±3, now unified to ±window).
    my $ref_continuity;
    if (@ref_ids) {
        my $gapless_refs = 0;
        for my $rid (@ref_ids) {
            my $gapless = 1;
            for my $col ($col_start .. $col_end) {
                my $c = substr($aln_seqs->{$rid}, $col, 1);
                if ($c eq '-' || $c eq '.') { $gapless = 0; last; }
            }
            $gapless_refs++ if $gapless;
        }
        $ref_continuity = $gapless_refs / scalar(@ref_ids);
    } else {
        $ref_continuity = 0.5;  # neutral when no ref available
    }

    my $score = $wc * $conservation + $wt * $ref_continuity + $wg * $gap_pattern;
    return ($score, $conservation, $ref_continuity, $gap_pattern);
}

# ---------------------------------------------------------------------------
# STEP K: PASS / FAIL / REVIEW assignment
# ---------------------------------------------------------------------------
print "\nSTEP K: Assigning PASS / FAIL / REVIEW...\n";

my %results;  # new_gene_id -> { overall_result, fail_reasons }

for my $gid (keys %merge_info) {
    my $transl_flag  = $has_internal_stop{$gid} ? "FRAMESHIFT_DETECTED" : "OK";
    my $mflag        = $msa_flag{$gid}      // ($no_msa ? "SKIPPED" : "NO_HIT");
    my $cov          = $gene_cov_metrics{$gid};
    my $mrgd_cov     = $cov->{merged_cov_by_ref};
    my $ref_cov      = $cov->{ref_cov_by_merged};
    my @fail_reasons;

    my $result;

    if ($transl_flag eq "FRAMESHIFT_DETECTED") {
        push @fail_reasons, "INTERNAL_STOP";
    }
    if ($mrgd_cov < $min_merged_cov && $mflag ne "GOOD_MSA") {
        push @fail_reasons, "LOW_MERGED_COV";
    }
    if ($mflag eq "NO_HIT" && $mrgd_cov < $min_merged_cov) {
        push @fail_reasons, "NO_HIT_LOW_COV" unless grep { $_ eq "LOW_MERGED_COV" } @fail_reasons;
    }

    if (@fail_reasons) {
        $result = "FAIL";
    } elsif ($transl_flag eq "OK"
          && $mrgd_cov >= $min_merged_cov
          && $ref_cov  >= $min_ref_cov
          && ($mflag eq "GOOD_MSA" || $mflag eq "GOOD_MSA_LOW_REF" || $no_msa)) {
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

# --- K1: TSV report ---
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
    "junction_scores",
    "min_junction_score",
    "mean_junction_score",
    "msa_flag",
    "translation_flag",
    "source_flags",
    "fail_reasons",
    "overall_result",
), "\n";

for my $gid (sort keys %merge_info) {
    my $info  = $merge_info{$gid};
    my $cov   = $gene_cov_metrics{$gid};
    my $res   = $results{$gid};
    my @jscores = @{ $junction_scores{$gid} // [] };

    my $min_js  = @jscores ? sprintf("%.4f", min(@jscores)) : ".";
    my $mean_js = @jscores ? sprintf("%.4f", sum(@jscores) / scalar(@jscores)) : ".";
    my $jscores_str = @jscores ? join("|", map { sprintf("%.4f", $_) } @jscores) : ".";

    my $fail_str = @{ $res->{fail_reasons} } ? join("|", @{ $res->{fail_reasons} }) : ".";
    my $mflag    = $msa_flag{$gid} // ($no_msa ? "SKIPPED" : "NO_HIT");
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
        $jscores_str,
        $min_js,
        $mean_js,
        $mflag,
        $res->{translation_flag},
        $src_flag,
        $fail_str,
        $res->{overall_result},
    ), "\n";
}
close $rep;
print "  Report:      $report_file\n";

# --- K2: Split GFF3 outputs ---
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

    # Add translation validation flags to the gene feature line (first line of block)
    if (@block) {
        my @jscores = @{ $junction_scores{$gid} // [] };
        my $min_js  = @jscores ? sprintf("%.4f", min(@jscores)) : ".";
        my $mflag   = $msa_flag{$gid} // ($no_msa ? "SKIPPED" : "NO_HIT");
        $block[0] .= ";transl_result=$res;msa_flag=$mflag;min_junction_score=$min_js";
    }
    print $fh "$_\n" for @block;
}
close $pfh; close $ffh; close $rfh;
printf "  Pass GFF:   %s  (%d genes)\n", $pass_gff,   $counts{PASS};
printf "  Fail GFF:   %s  (%d genes)\n", $fail_gff,   $counts{FAIL};
printf "  Review GFF: %s  (%d genes)\n", $review_gff, $counts{REVIEW};

# --- K3: Pass proteins FASTA ---
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
