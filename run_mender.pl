#!/usr/bin/perl
use strict;
use warnings;

## run_mender.pl — Mender pipeline wrapper
##
## Reads a mender.cfg config file and runs the split gene detection
## and merging pipeline end to end.
##
## Usage:
##   perl run_mender.pl --config mender.cfg
##   perl run_mender.pl --config mender.cfg --steps 1,2,3,4
##   perl run_mender.pl --config mender.cfg --steps 6
##   perl run_mender.pl --config mender.cfg --dry_run
##
## Steps:
##   1  prepare    clean protein fasta, extract GFF subsets
##   2  diamond    run diamond blastp
##   3  bedtools   run bedtools intersect for IsoSeq overlaps
##   4  find       run find_split_genes.pl
##   5  isoseq     run validate_merge_with_isoseq.pl
##   6  merge      run merge_split_genes.pl
##   7  check      run gt gff3validator and agat validation
##
## If --steps is not specified, all applicable steps are run in order.
## Steps 3 and 5 are skipped automatically if isoseq_gff is not set.

use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use Cwd qw(abs_path);

my $config_file = undef;
my $steps_opt   = undef;
my $dry_run     = 0;

GetOptions(
    "config=s" => \$config_file,
    "steps=s"  => \$steps_opt,
    "dry_run"  => \$dry_run,
) or die "usage: $0 --config mender.cfg [--steps 1,2,3] [--dry_run]\n";

die "usage: $0 --config mender.cfg\n" unless defined $config_file;
die "Config file not found: $config_file\n" unless -f $config_file;

# ---------------------------------------------------------------------------
# Parse config file
# ---------------------------------------------------------------------------
my %cfg;
my $current_section = "general";

open CFG, $config_file or die "Cannot open config: $config_file $!\n";
while (my $line = <CFG>) {
    chomp $line;
    $line =~ s/#.*//;        # strip comments
    $line =~ s/^\s+|\s+$//g; # trim whitespace
    next unless length $line;

    if ($line =~ /^\[(\w+)\]$/) {
        $current_section = lc($1);
        next;
    }

    if ($line =~ /^(\w+)\s*=\s*(.*)$/) {
        my ($key, $val) = ($1, $2);
        $val =~ s/\s+$//;
        $cfg{$current_section}{$key} = $val;
    }
}
close CFG;

# ---------------------------------------------------------------------------
# Resolve config values with defaults
# ---------------------------------------------------------------------------
sub cfg {
    my ($section, $key, $default) = @_;
    return $cfg{$section}{$key} // $default // "";
}

my $date_tag = strftime("%Y%m%d", localtime);

# input
my $gff           = cfg("input", "gff");
my $proteome_fa   = cfg("input", "proteome_fa");
my $subject_fa    = cfg("input", "subject_fa");
my $isoseq_gff    = cfg("input", "isoseq_gff");

# output
my $output_gff    = cfg("output", "output_gff",    "new_merges.gff");
my $removed_gff   = cfg("output", "removed_gff",   "removed_genes.gff");
my $workdir       = cfg("output", "workdir",        "mender_workdir");

# diamond
my $threads       = cfg("diamond", "threads",  "4");
my $evalue        = cfg("diamond", "evalue",   "1e-20");

# ids — substitute [DATE] token
my $gene_template = cfg("ids", "gene_template",
                        "MERGE${date_tag}g[GCOUNT:6]000");
my $trans_template = cfg("ids", "trans_template",
                         "MERGE${date_tag}t[GCOUNT:6][TCOUNT:3]");
$gene_template  =~ s/\[DATE\]/$date_tag/g;
$trans_template =~ s/\[DATE\]/$date_tag/g;

# merge filters
my $fix_partial        = cfg("merge_filters", "fix_partial",        "yes");
my $skip_flags         = cfg("merge_filters", "skip_flags",         "SKIPPED_GENE,LOW_COV");
my $flags              = cfg("merge_filters", "flags",              "all");
my $min_tiling         = cfg("merge_filters", "min_tiling",         "1");
my $min_cov            = cfg("merge_filters", "min_cov",            "0");
my $require_isoseq     = cfg("merge_filters", "require_isoseq",     "");
my $isoseq_min_span    = cfg("merge_filters", "isoseq_min_spanning","0");

# validation
my $run_gt    = cfg("validation", "run_gt",   "yes");
my $run_agat  = cfg("validation", "run_agat", "yes");

# paths
my $scripts_dir  = cfg("paths", "scripts_dir",  "") || dirname(abs_path($0));
my $diamond_bin  = cfg("paths", "diamond_bin",  "") || "diamond";
my $bedtools_bin = cfg("paths", "bedtools_bin", "") || "bedtools";
my $gt_bin       = cfg("paths", "gt_bin",       "") || "gt";
my $agat_bin     = cfg("paths", "agat_bin",     "") || "agat_convert_sp_gxf2gxf.pl";
my $perl_bin     = cfg("paths", "perl_bin",     "") || "perl";

# intermediate file paths
my $clean_fa         = "$workdir/protein.nostops.fa";
my $diamond_db       = "$workdir/subject.dmnd";
my $diamond_out      = "$workdir/diamond.out";
my $genes_gff        = "$workdir/models.gene.gff";
my $isoseq_mrna_gff  = "$workdir/isoseq.mrna.gff";
my $overlaps         = "$workdir/all.overlaps";
my $merge_table      = "$workdir/split_genes_merge.txt";
my $validated_table  = "$workdir/validated_merge.txt";
my $find_script      = "$scripts_dir/find_split_genes.pl";
my $isoseq_script    = "$scripts_dir/validate_merge_with_isoseq.pl";
my $merge_script     = "$scripts_dir/merge_split_genes.pl";

# ---------------------------------------------------------------------------
# Determine which steps to run
# ---------------------------------------------------------------------------
my %run_steps;
if (defined $steps_opt) {
    for my $s (split /,/, $steps_opt) {
        $run_steps{$s+0} = 1;
    }
} else {
    $run_steps{$_} = 1 for (1..7);
}

# skip IsoSeq steps if no isoseq_gff provided
if (!$isoseq_gff) {
    delete $run_steps{3};
    delete $run_steps{5};
    print "INFO: isoseq_gff not set — skipping steps 3 and 5\n";
    # merge step uses merge_table not validated_table
    $validated_table = $merge_table;
}

# ---------------------------------------------------------------------------
# Helper: run a command, print it, check exit code
# ---------------------------------------------------------------------------
sub run_cmd {
    my ($cmd, $desc) = @_;
    print "\n--- $desc ---\n";
    print "CMD: $cmd\n";
    return if $dry_run;
    my $ret = system($cmd);
    die "ERROR: command failed (exit $ret): $desc\n" if $ret != 0;
}

sub check_file {
    my ($file, $desc) = @_;
    die "ERROR: required file not found: $file ($desc)\n"
        unless $dry_run || -f $file;
}

# ---------------------------------------------------------------------------
# Print run summary
# ---------------------------------------------------------------------------
print "=" x 60, "\n";
print "MENDER PIPELINE\n";
print "=" x 60, "\n";
print "Config:          $config_file\n";
print "Date:            $date_tag\n";
print "Steps:           ", join(", ", sort { $a <=> $b } keys %run_steps), "\n";
print "Dry run:         ", ($dry_run ? "YES" : "no"), "\n";
print "-" x 60, "\n";
print "Input GFF:       $gff\n";
print "Proteome FA:     $proteome_fa\n";
print "Subject FA:      $subject_fa\n";
print "IsoSeq GFF:      ", ($isoseq_gff || "(not provided)"), "\n";
print "Output GFF:      $output_gff\n";
print "Workdir:         $workdir\n";
print "Gene template:   $gene_template\n";
print "Trans template:  $trans_template\n";
print "Skip flags:      $skip_flags\n";
print "Fix partial:     $fix_partial\n";
print "=" x 60, "\n";

# create workdir
unless ($dry_run) {
    mkdir $workdir unless -d $workdir;
    die "ERROR: could not create workdir: $workdir\n" unless -d $workdir;
}

# ---------------------------------------------------------------------------
# STEP 1: PREPARE
# ---------------------------------------------------------------------------
if ($run_steps{1}) {
    print "\n", "=" x 60, "\n";
    print "STEP 1: PREPARE INPUTS\n";
    print "=" x 60, "\n";

    check_file($proteome_fa, "proteome_fa");
    check_file($gff,         "gff");

    # remove internal stop codons from protein fasta
    run_cmd(
        "grep -v '\\*' $proteome_fa | " .
        "grep -v -P '^[^>].*\\.' > $clean_fa",
        "Remove internal stops from protein fasta"
    );

    # extract gene features for bedtools
    run_cmd(
        "grep -P '\\tgene\\t' $gff > $genes_gff",
        "Extract gene features from GFF"
    );

    if ($isoseq_gff) {
        check_file($isoseq_gff, "isoseq_gff");
        run_cmd(
            "grep -P '\\tmRNA\\t' $isoseq_gff > $isoseq_mrna_gff",
            "Extract IsoSeq mRNA features"
        );
    }
}

# ---------------------------------------------------------------------------
# STEP 2: DIAMOND
# ---------------------------------------------------------------------------
if ($run_steps{2}) {
    print "\n", "=" x 60, "\n";
    print "STEP 2: DIAMOND BLASTP\n";
    print "=" x 60, "\n";

    check_file($subject_fa, "subject_fa");
    check_file($clean_fa,   "cleaned proteome (run step 1 first)");

    run_cmd(
        "$diamond_bin makedb --in $subject_fa --db $diamond_db",
        "Build diamond database"
    );

    run_cmd(
        "$diamond_bin blastp " .
        "--db $diamond_db " .
        "--query $clean_fa " .
        "--outfmt 6 qseqid sseqid pident length mismatch gapopen " .
        "qstart qend sstart send evalue bitscore scovhsp slen " .
        "--threads $threads " .
        "--evalue $evalue " .
        "--out $diamond_out",
        "Run diamond blastp"
    );
}

# ---------------------------------------------------------------------------
# STEP 3: BEDTOOLS (IsoSeq overlaps)
# ---------------------------------------------------------------------------
if ($run_steps{3} && $isoseq_gff) {
    print "\n", "=" x 60, "\n";
    print "STEP 3: BEDTOOLS INTERSECT (IsoSeq overlaps)\n";
    print "=" x 60, "\n";

    check_file($isoseq_mrna_gff, "isoseq mRNA GFF (run step 1 first)");
    check_file($genes_gff,       "gene GFF (run step 1 first)");

    run_cmd(
        "$bedtools_bin intersect " .
        "-a $isoseq_mrna_gff -b $genes_gff -wa -wb | " .
        "$perl_bin -p -e " .
        "'s/.+\\sID=([^;]+).+\\sID=([^;]+).*\$/\$1\\t\$2/' " .
        "> $overlaps",
        "Bedtools intersect IsoSeq reads vs gene models"
    );
}

# ---------------------------------------------------------------------------
# STEP 4: FIND SPLIT GENES
# ---------------------------------------------------------------------------
if ($run_steps{4}) {
    print "\n", "=" x 60, "\n";
    print "STEP 4: FIND SPLIT GENES\n";
    print "=" x 60, "\n";

    check_file($diamond_out,  "diamond output (run step 2 first)");
    check_file($gff,          "gff");
    check_file($subject_fa,   "subject_fa");
    check_file($clean_fa,     "cleaned proteome (run step 1 first)");
    check_file($find_script,  "find_split_genes.pl");

    run_cmd(
        "$perl_bin $find_script " .
        "$diamond_out " .
        "$gff " .
        "$subject_fa " .
        "$clean_fa",
        "Find split gene candidates"
    );

    # move outputs to workdir
    run_cmd("mv split_genes_summary.txt split_genes_detail.txt " .
            "split_genes_merge.txt $workdir/",
            "Move find_split_genes outputs to workdir");
    $merge_table = "$workdir/split_genes_merge.txt";
}

# ---------------------------------------------------------------------------
# STEP 5: VALIDATE WITH ISOSEQ
# ---------------------------------------------------------------------------
if ($run_steps{5} && $isoseq_gff) {
    print "\n", "=" x 60, "\n";
    print "STEP 5: VALIDATE WITH ISOSEQ\n";
    print "=" x 60, "\n";

    check_file($merge_table,   "merge table (run step 4 first)");
    check_file($overlaps,      "overlaps (run step 3 first)");
    check_file($isoseq_script, "validate_merge_with_isoseq.pl");

    run_cmd(
        "$perl_bin $isoseq_script $merge_table $overlaps " .
        "> $validated_table",
        "Validate merge candidates with IsoSeq spanning reads"
    );
}

# ---------------------------------------------------------------------------
# STEP 6: MERGE
# ---------------------------------------------------------------------------
if ($run_steps{6}) {
    print "\n", "=" x 60, "\n";
    print "STEP 6: MERGE SPLIT GENES\n";
    print "=" x 60, "\n";

    check_file($validated_table, "merge/validated table (run steps 4-5 first)");
    check_file($gff,             "gff");
    check_file($merge_script,    "merge_split_genes.pl");

    # build merge command
    my $merge_cmd = "$perl_bin $merge_script";
    $merge_cmd .= " --fix_partial"              if $fix_partial eq "yes";
    $merge_cmd .= " --skip_flags $skip_flags"   if $skip_flags;
    $merge_cmd .= " --flags $flags"             if $flags && $flags ne "all";
    $merge_cmd .= " --min_tiling $min_tiling"   if $min_tiling > 1;
    $merge_cmd .= " --min_cov $min_cov"         if $min_cov > 0;
    $merge_cmd .= " --require_isoseq $require_isoseq"  if $require_isoseq;
    $merge_cmd .= " --isoseq_min_spanning $isoseq_min_span" if $isoseq_min_span > 0;
    $merge_cmd .= " --removed $removed_gff";
    $merge_cmd .= " --gene_template  \"$gene_template\"";
    $merge_cmd .= " --trans_template \"$trans_template\"";
    $merge_cmd .= " $validated_table $gff $output_gff";

    run_cmd($merge_cmd, "Merge split genes into annotation");
}

# ---------------------------------------------------------------------------
# STEP 7: VALIDATE OUTPUT
# ---------------------------------------------------------------------------
if ($run_steps{7}) {
    print "\n", "=" x 60, "\n";
    print "STEP 7: VALIDATE OUTPUT\n";
    print "=" x 60, "\n";

    check_file($output_gff, "output GFF (run step 6 first)");

    # gene count sanity check
    print "\n--- Gene count sanity check ---\n";
    unless ($dry_run) {
        my $before = `grep -c \$'\\tgene\\t' $gff`;
        my $after  = `grep -c \$'\\tgene\\t' $output_gff`;
        my $removed = `grep -c \$'\\tgene\\t' $removed_gff`;
        chomp ($before, $after, $removed);
        my $log_file = $output_gff;
        $log_file =~ s/\.gff$/.log/i;
        my $merged = `grep -c '^merge_' $log_file 2>/dev/null` || 0;
        chomp $merged;
        print "  Input genes:    $before\n";
        print "  Output genes:   $after\n";
        print "  Removed genes:  $removed\n";
        print "  Merged genes:   $merged\n";
        my $expected = $before - $removed + $merged;
        print "  Expected after: $expected\n";
        if ($after == $expected) {
            print "  Gene count: OK\n";
        } else {
            print "  WARNING: gene count mismatch! Expected $expected got $after\n";
        }
    }

    # gt gff3validator
    if ($run_gt eq "yes") {
        run_cmd("$gt_bin gff3validator $output_gff",
                "GFF3 spec validation (gt gff3validator)");
    }

    # agat biological validation
    if ($run_agat eq "yes") {
        run_cmd(
            "$agat_bin -g $output_gff -o /dev/null 2>&1 | " .
            "grep -i 'error\\|warn' | head -20",
            "Biological consistency validation (AGAT)"
        );
    }

    # Mender feature count
    print "\n--- Mender-created features ---\n";
    unless ($dry_run) {
        my $mender_genes = `grep -P '\\tMend\\t' $output_gff | grep -cP '\\tgene\\t'`;
        chomp $mender_genes;
        print "  New Mender genes: $mender_genes\n";
    }
}

print "\n", "=" x 60, "\n";
print $dry_run ? "DRY RUN COMPLETE\n" : "MENDER PIPELINE COMPLETE\n";
print "Output: $output_gff\n" unless $dry_run;
print "=" x 60, "\n";
