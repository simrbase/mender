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
##   7  gt         run gt gff3validator on merged GFF (fast; catches frame issues)
##   8  transl     run validate_merge_translation.pl → pass/fail/review GFFs
##   9  agat       run agat validation on pass GFF only (slow)
##
## If --steps is not specified, all applicable steps are run in order.
## Steps 3 and 5 are skipped automatically if isoseq_gff is not set.
## Step 8 is skipped automatically if run_translation_validation = no in config.
## Step 9 is skipped automatically if run_agat = no in config.

use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use Cwd qw(abs_path);

my @argv_orig   = @ARGV;   # capture before GetOptions consumes it
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

my $date_tag    = strftime("%Y%m%d",       localtime);
my $run_start   = time();
my $start_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime);

# input
my $gff           = cfg("input", "gff");
my $proteome_fa   = cfg("input", "proteome_fa");
my $subject_fa    = cfg("input", "subject_fa");
my $genome_fa     = cfg("input", "genome_fa");
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
my $asym_trim          = cfg("merge_filters", "asym_trim",          "yes");
my $large_span_warn    = cfg("merge_filters", "large_span_warn",    "500000");
my $large_span_hard    = cfg("merge_filters", "large_span_hard",    "2000000");

# validation (gt + agat)
my $run_gt    = cfg("validation", "run_gt",   "yes");
my $run_agat  = cfg("validation", "run_agat", "yes");

# translation validation (step 8)
my $run_transl_val  = cfg("validation", "run_translation_validation", "yes");
my $tv_min_junction = cfg("validation", "min_junction_score",         "0.5");
my $tv_min_mrgd_cov = cfg("validation", "min_merged_cov",             "0.60");
my $tv_min_ref_cov  = cfg("validation", "min_ref_cov",                "0.50");
my $tv_keep_msa     = cfg("validation", "keep_msa",                   "no");
my $tv_no_msa       = cfg("validation", "no_msa",                     "no");
my $tv_out_prefix   = cfg("validation", "transl_out_prefix",          "");
my $tv_swissprot_fa = cfg("validation", "swissprot_fa",               "");
my $tv_swissprot_db = cfg("validation", "swissprot_db",               "");
my $tv_aligner      = cfg("validation", "aligner",                    "mafft_fast");
my $tv_final_incl_review = cfg("validation", "final_gff_include_review", "yes");

# paths
my $scripts_dir  = cfg("paths", "scripts_dir",  "") || dirname(abs_path($0));
my $diamond_bin  = cfg("paths", "diamond_bin",  "") || "diamond";
my $bedtools_bin = cfg("paths", "bedtools_bin", "") || "bedtools";
my $gt_bin       = cfg("paths", "gt_bin",       "") || "gt";
my $agat_bin     = cfg("paths", "agat_bin",     "") || "agat_convert_sp_gxf2gxf.pl";
my $gffread_bin  = cfg("paths", "gffread_bin",  "") || "gffread";
my $mafft_bin    = cfg("paths", "mafft_bin",    "") || "mafft";
my $kalign_bin   = cfg("paths", "kalign_bin",   "") || "kalign";
my $perl_bin     = cfg("paths", "perl_bin",     "") || "perl";

# intermediate file paths
my $clean_fa         = "$workdir/protein.nostops.fa";
my $diamond_db       = "$workdir/subject.dmnd";
my $diamond_out      = "$workdir/diamond.out";
my $genes_gff        = "$workdir/models.gene.gff";
my $isoseq_mrna_gff  = "$workdir/isoseq.mrna.gff";
my $overlaps         = "$workdir/all.overlaps";
my $merge_table      = "$workdir/merge_candidates.txt";
my $validated_table  = "$workdir/isoseq_validated.txt";
my $find_script      = "$scripts_dir/find_split_genes.pl";
my $isoseq_script    = "$scripts_dir/validate_merge_with_isoseq.pl";
my $merge_script     = "$scripts_dir/merge_split_genes.pl";
my $transl_script    = "$scripts_dir/validate_merge_translation.pl";

# ---------------------------------------------------------------------------
# Determine which steps to run
# ---------------------------------------------------------------------------
my %run_steps;
my %step_skip;   # step -> reason string (for run report)
my @steps_requested;

if (defined $steps_opt) {
    for my $s (split /,/, $steps_opt) {
        $run_steps{$s+0} = 1;
    }
    @steps_requested = sort { $a <=> $b } keys %run_steps;
    # steps not in --steps are "not requested"
    for my $s (1..9) {
        $step_skip{$s} = "not in --steps" unless $run_steps{$s};
    }
} else {
    $run_steps{$_} = 1 for (1..9);
    @steps_requested = (1..9);
}

# skip IsoSeq steps if no isoseq_gff provided
if (!$isoseq_gff) {
    delete $run_steps{3};
    delete $run_steps{5};
    $step_skip{3} = "isoseq_gff not set";
    $step_skip{5} = "isoseq_gff not set";
    print "INFO: isoseq_gff not set — skipping steps 3 and 5\n";
    # merge step uses merge_table not validated_table
    $validated_table = $merge_table;
}

# skip translation validation if disabled in config
if ($run_transl_val =~ /^no$/i) {
    delete $run_steps{8};
    $step_skip{8} = "run_translation_validation = no";
    print "INFO: run_translation_validation = no — skipping step 8\n";
}

# skip agat if disabled in config
if ($run_agat =~ /^no$/i) {
    delete $run_steps{9};
    $step_skip{9} = "run_agat = no";
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
print "Large span warn: $large_span_warn bp\n";
print "Large span hard: $large_span_hard bp\n";
print "Transl val:      ", ($run_transl_val =~ /^no$/i ? "no" : "yes"), "\n";
if ($run_transl_val !~ /^no$/i) {
    print "  min_junction:  $tv_min_junction\n";
    print "  min_merged_cov:$tv_min_mrgd_cov\n";
    print "  min_ref_cov:   $tv_min_ref_cov\n";
    print "  no_msa:        ", ($tv_no_msa =~ /^yes$/i ? "yes" : "no"), "\n";
    print "  aligner:       $tv_aligner\n";
    print "  swissprot_fa:  ", ($tv_swissprot_fa || "(not set — MSA will use ref_fa)"), "\n";
}
print "=" x 60, "\n";

# ---------------------------------------------------------------------------
# Tool preflight check
# ---------------------------------------------------------------------------
{
    sub tool_ok {
        my ($bin) = @_;
        return -x $bin if $bin =~ m{/};   # full path: check executable
        return system("which \Q$bin\E > /dev/null 2>&1") == 0;
    }

    my @missing;
    my @checks;  # [tool_bin, label, condition]

    # diamond: always needed (step 2, step 8 db reuse)
    push @checks, [$diamond_bin,  "diamond",  1];

    # bedtools: only if IsoSeq GFF is set (step 3)
    push @checks, [$bedtools_bin, "bedtools", !!$isoseq_gff];

    # gt: only if run_gt = yes (step 7)
    push @checks, [$gt_bin,       "gt",       $run_gt  !~ /^no$/i];

    # gffread: only if translation validation is on (step 8)
    push @checks, [$gffread_bin,  "gffread",  $run_transl_val !~ /^no$/i];

    # mafft: only if translation validation is on, no_msa is off, and aligner is mafft_*
    push @checks, [$mafft_bin,    "mafft",    $run_transl_val !~ /^no$/i
                                              && $tv_no_msa =~ /^no$/i
                                              && $tv_aligner =~ /^mafft/];

    # kalign3: only if translation validation is on, no_msa is off, and aligner = kalign3
    push @checks, [$kalign_bin,   "kalign3",  $run_transl_val !~ /^no$/i
                                              && $tv_no_msa =~ /^no$/i
                                              && $tv_aligner eq "kalign3"];

    # agat: only if run_agat = yes (step 9)
    push @checks, [$agat_bin,     "agat",     $run_agat !~ /^no$/i];

    print "\n--- Tool preflight check ---\n";
    for my $c (@checks) {
        my ($bin, $label, $needed) = @$c;
        next unless $needed;
        if (tool_ok($bin)) {
            printf "  %-12s OK  (%s)\n", $label, $bin;
        } else {
            printf "  %-12s MISSING (%s)\n", $label, $bin;
            push @missing, "$label ($bin)";
        }
    }

    if (@missing) {
        print "\nERROR: required tools not found:\n";
        print "  $_\n" for @missing;
        print "Load the relevant modules or set bin paths in the [paths] config section.\n";
        exit 1;
    }
    print "  All required tools found.\n";

    # File existence check for swissprot_fa (warning only — falls back to ref_fa)
    if ($run_transl_val !~ /^no$/i && $tv_swissprot_fa) {
        if (-f $tv_swissprot_fa) {
            printf "  %-12s OK  (%s)\n", "swissprot_fa", $tv_swissprot_fa;
        } else {
            print "WARNING: swissprot_fa not found: $tv_swissprot_fa\n";
            print "         MSA will fall back to ref_fa sequences for junction scoring.\n";
            $tv_swissprot_fa = "";  # clear so it is not passed to the script
        }
    }
}

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
        "qstart qend sstart send evalue bitscore scovhsp slen qcovhsp qlen " .
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
        ($asym_trim =~ /^no$/i ? "--no_asym_trim " : "") .
        ($large_span_warn  ? "--large_span_warn $large_span_warn "       : "") .
        ($large_span_hard  ? "--large_span_extreme $large_span_hard "    : "") .
        "$diamond_out " .
        "$gff " .
        "$subject_fa " .
        "$clean_fa",
        "Find split gene candidates"
    );

    # move outputs to workdir
    run_cmd("mv split_genes_summary.txt split_genes_detail.txt $workdir/ && " .
            "mv split_genes_merge.txt $workdir/merge_candidates.txt",
            "Move find_split_genes outputs to workdir");
    $merge_table = "$workdir/merge_candidates.txt";
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
# STEP 7: GT GFF3VALIDATOR (fast; Mender-created genes only)
# ---------------------------------------------------------------------------
if ($run_steps{7}) {
    print "\n", "=" x 60, "\n";
    print "STEP 7: GT GFF3VALIDATOR (Mender genes only)\n";
    print "=" x 60, "\n";

    check_file($output_gff, "output GFF (run step 6 first)");

    # gene count sanity check
    print "\n--- Gene count sanity check ---\n";
    unless ($dry_run) {
        my $before  = `grep -c \$'\\tgene\\t' $gff`;
        my $after   = `grep -c \$'\\tgene\\t' $output_gff`;
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
        my $mender_genes = `grep -P '\\tMender\\t' $output_gff | grep -cP '\\tgene\\t'`;
        chomp $mender_genes;
        print "  New Mender genes: $mender_genes\n";
    }

    if ($run_gt eq "yes") {
        # Extract only Mender-created features into a temp GFF for validation.
        # This isolates merge-introduced errors from pre-existing annotation issues.
        my $mender_only_gff = "$workdir/mender_only.gff";
        print "\n--- GT GFF3 validation (Mender genes only) ---\n";
        unless ($dry_run) {
            open my $mo, ">", $mender_only_gff
                or die "Cannot write $mender_only_gff: $!\n";
            print $mo "##gff-version 3\n";
            open my $fh, "<", $output_gff
                or die "Cannot open $output_gff: $!\n";
            while (<$fh>) { print $mo $_ if /\tMender\t/ }
            close $fh;
            close $mo;

            my $n_mender_lines = `grep -c '' $mender_only_gff` + 0;
            chomp $n_mender_lines;
            print "CMD: $gt_bin gff3validator $mender_only_gff\n";
            print "NOTE: errors are reported for review only — pipeline continues\n";
            system("$gt_bin gff3validator $mender_only_gff");
            my $gt_exit = $? >> 8;
            if ($gt_exit == 0) {
                print "  gt gff3validator: OK — no errors in Mender genes\n";
            } else {
                print "  gt gff3validator: errors found in Mender genes (see above)\n";
            }
        }
    }
}

# ---------------------------------------------------------------------------
# STEP 8: TRANSLATION VALIDATION
# ---------------------------------------------------------------------------
if ($run_steps{8}) {
    print "\n", "=" x 60, "\n";
    print "STEP 8: TRANSLATION VALIDATION\n";
    print "=" x 60, "\n";

    check_file($output_gff,   "merged GFF (run step 6 first)");
    check_file($gff,          "original GFF");
    check_file($genome_fa,    "genome_fa");
    check_file($subject_fa,   "subject_fa (ref proteome)");
    check_file($transl_script, "validate_merge_translation.pl");

    # derive output prefix: use transl_out_prefix from config, else workdir/transl
    my $tv_out = $tv_out_prefix || "$workdir/transl";

    my $tv_cmd = "$perl_bin $transl_script"
        . " --merged_gff $output_gff"
        . " --orig_gff $gff"
        . " --genome_fa $genome_fa"
        . " --ref_fa $subject_fa"
        . " --diamond_db $diamond_db"
        . " --merge_table $validated_table"
        . " --out $tv_out"
        . " --threads $threads"
        . " --min_junction $tv_min_junction"
        . " --min_merged_cov $tv_min_mrgd_cov"
        . " --min_ref_cov $tv_min_ref_cov"
        . ($tv_keep_msa =~ /^yes$/i ? " --keep_msa"  : "")
        . ($tv_no_msa =~ /^yes$/i ? " --no_msa"  : "")
        . ($tv_swissprot_fa ? " --swissprot_fa $tv_swissprot_fa" : "")
        . ($tv_swissprot_db ? " --swissprot_db $tv_swissprot_db" : "")
        . " --aligner $tv_aligner"
        . " --gffread_bin $gffread_bin"
        . " --diamond_bin $diamond_bin"
        . " --mafft_bin $mafft_bin"
        . " --kalign_bin $kalign_bin";

    run_cmd($tv_cmd, "Translation validation of merged genes");

    unless ($dry_run) {
        print "\n--- Translation validation summary ---\n";
        my $pass_gff   = "${tv_out}_pass.gff3";
        my $report_tsv = "${tv_out}_report.tsv";
        if (-f $report_tsv) {
            my $pass   = `grep -c 'PASS'   $report_tsv 2>/dev/null` || 0;
            my $fail   = `grep -c 'FAIL'   $report_tsv 2>/dev/null` || 0;
            my $review = `grep -c 'REVIEW' $report_tsv 2>/dev/null` || 0;
            chomp ($pass, $fail, $review);
            # subtract header line matches if any
            print "  Report:   $report_tsv\n";
            print "  Pass GFF: ${tv_out}_pass.gff3\n";
            print "  Fail GFF: ${tv_out}_fail.gff3\n";
        }
    }
}

# ---------------------------------------------------------------------------
# STEP 9: AGAT BIOLOGICAL VALIDATION (on pass GFF only)
# ---------------------------------------------------------------------------
if ($run_steps{9}) {
    print "\n", "=" x 60, "\n";
    print "STEP 9: AGAT BIOLOGICAL VALIDATION\n";
    print "=" x 60, "\n";

    my $tv_out    = $tv_out_prefix || "$workdir/transl";
    my $pass_gff  = "${tv_out}_pass.gff3";

    # fall back to full merged GFF if translation validation was not run
    my $target_gff = (-f $pass_gff) ? $pass_gff : $output_gff;
    if ($target_gff eq $output_gff) {
        print "INFO: pass GFF not found — running AGAT on full merged GFF\n";
        print "      (run step 8 first to validate on pass merges only)\n";
    } else {
        print "  Running AGAT on: $pass_gff\n";
    }

    check_file($target_gff, "GFF for AGAT validation");

    run_cmd(
        "$agat_bin -g $target_gff -o /dev/null 2>&1 | " .
        "grep -i 'error\\|warn' | head -20",
        "Biological consistency validation (AGAT)"
    );
}

# ---------------------------------------------------------------------------
# FINAL VALIDATED GFF
# ---------------------------------------------------------------------------
# Builds a complete annotation GFF where FAIL merges (and optionally REVIEW)
# are replaced with their original source genes. Adds transl_result, msa_flag,
# and min_junction_score attributes to all retained Mender gene features.
# Generated automatically whenever translation validation outputs are present.
# ---------------------------------------------------------------------------
{
    my $tv_out     = $tv_out_prefix || "$workdir/transl";
    my $fail_gff   = "${tv_out}_fail.gff3";
    my $report_tsv = "${tv_out}_report.tsv";

    if (-f $fail_gff && !$dry_run) {
        print "\n", "=" x 60, "\n";
        print "FINAL VALIDATED GFF\n";
        print "=" x 60, "\n";

        # Load validation attributes from TSV report (result, msa_flag, min_junction_score)
        my %transl_attrs;  # new_gene_id -> hashref
        if (-f $report_tsv) {
            open my $fh, "<", $report_tsv or die "Cannot open $report_tsv: $!\n";
            <$fh>;  # skip header
            while (my $line = <$fh>) {
                chomp $line;
                my @f = split /\t/, $line;
                next unless @f >= 18;
                # cols: 1=new_gene_id, 13=min_junction_score, 15=msa_flag, 17=overall_result
                $transl_attrs{$f[1]} = {
                    result   => $f[17],
                    msa_flag => $f[15],
                    min_js   => $f[13],
                };
            }
            close $fh;
        }

        # Determine which Mender genes to exclude (replace with source genes)
        my %exclude;   # mender_gene_id -> [source_gene_ids]
        my $n_fail_excl   = 0;
        my $n_review_excl = 0;

        my $parse_excl_gff = sub {
            my ($gff_file) = @_;
            open my $fh, "<", $gff_file or return;
            while (<$fh>) {
                next unless /\tgene\t/;
                my $attr_str = (split /\t/)[8] // "";
                my ($gid)   = $attr_str =~ /\bID=([^;]+)/;
                my ($mfrom) = $attr_str =~ /\bmerged_from=([^;]+)/;
                next unless $gid && $mfrom;
                $exclude{$gid} = [split /,/, $mfrom];
            }
            close $fh;
        };

        $parse_excl_gff->($fail_gff);
        $n_fail_excl = scalar keys %exclude;

        my $review_gff = "${tv_out}_review.gff3";
        if ($tv_final_incl_review =~ /^no$/i && -f $review_gff) {
            my $before = scalar keys %exclude;
            $parse_excl_gff->($review_gff);
            $n_review_excl = (scalar keys %exclude) - $before;
        }

        # Load source gene blocks from removed_genes.gff
        my %removed_blocks;  # source_gene_id -> [gff lines]
        if (-f $removed_gff) {
            open my $fh, "<", $removed_gff or die "Cannot open $removed_gff: $!\n";
            my ($cur_id, @cur_block);
            while (my $line = <$fh>) {
                chomp $line;
                next if $line =~ /^#/ || $line =~ /^\s*$/;
                if ($line =~ /\tgene\t/) {
                    $removed_blocks{$cur_id} = [@cur_block] if $cur_id && @cur_block;
                    @cur_block = ($line);
                    ($cur_id) = $line =~ /\bID=([^;]+)/;
                } else {
                    push @cur_block, $line;
                }
            }
            $removed_blocks{$cur_id} = [@cur_block] if $cur_id && @cur_block;
            close $fh;
        }

        # Write final validated GFF
        (my $validated_gff = $output_gff) =~ s/(\.[^.]+)$/_validated$1/;
        {
            open my $in,  "<", $output_gff   or die "Cannot open $output_gff: $!\n";
            open my $out, ">", $validated_gff or die "Cannot write $validated_gff: $!\n";

            my ($cur_id, @cur_block, $cur_is_mender, $cur_is_excluded);

            my $flush = sub {
                return unless @cur_block;
                if ($cur_is_excluded) {
                    # Replace excluded Mender gene with its original source gene blocks
                    for my $src_id (@{ $exclude{$cur_id} }) {
                        print $out "$_\n" for @{ $removed_blocks{$src_id} // [] };
                    }
                } else {
                    # Annotate retained Mender genes with transl validation flags
                    if ($cur_is_mender && exists $transl_attrs{$cur_id}) {
                        my $ta = $transl_attrs{$cur_id};
                        $cur_block[0] .= ";transl_result=$ta->{result}" .
                                         ";msa_flag=$ta->{msa_flag}" .
                                         ";min_junction_score=$ta->{min_js}";
                    }
                    print $out "$_\n" for @cur_block;
                }
                @cur_block = ();
                ($cur_id, $cur_is_mender, $cur_is_excluded) = (undef, 0, 0);
            };

            while (my $line = <$in>) {
                chomp $line;
                if ($line =~ /^\s*$/ || $line =~ /^#/) {
                    $flush->();
                    print $out "$line\n";
                } elsif ($line =~ /\tgene\t/) {
                    $flush->();
                    push @cur_block, $line;
                    ($cur_id)        = $line =~ /\bID=([^;]+)/;
                    $cur_is_mender   = ($line =~ /merge_source=Mender/) ? 1 : 0;
                    $cur_is_excluded = (exists $exclude{$cur_id // ""}) ? 1 : 0;
                } else {
                    push @cur_block, $line;
                }
            }
            $flush->();
            close $in; close $out;
        }

        printf "  FAIL merges replaced with source genes:   %d\n", $n_fail_excl;
        printf "  REVIEW merges replaced with source genes: %d\n", $n_review_excl
            if $tv_final_incl_review =~ /^no$/i;
        printf "  REVIEW merges retained (PASS-equivalent): %d\n",
            (scalar(grep { $transl_attrs{$_}{result} eq "REVIEW" } keys %transl_attrs) - $n_review_excl)
            if %transl_attrs;
        print  "  Final validated GFF: $validated_gff\n";
    }
}

# ---------------------------------------------------------------------------
# RUN REPORT
# ---------------------------------------------------------------------------
# Written alongside the output GFF. Captures all parameters, step status,
# result counts, and output file inventory so the run is fully reproducible
# and auditable without re-reading the terminal log.
# ---------------------------------------------------------------------------
{
    my $end_stamp  = strftime("%Y-%m-%d %H:%M:%S", localtime);
    my $elapsed    = time() - $run_start;
    my $elapsed_str = sprintf "%dh %02dm %02ds",
                      int($elapsed/3600), int(($elapsed%3600)/60), $elapsed%60;

    my $tv_out = $tv_out_prefix || "$workdir/transl";
    (my $validated_gff = $output_gff) =~ s/(\.[^.]+)$/_validated$1/;
    (my $report_file   = $output_gff) =~ s/(\.[^.]+)$/_run_report.txt/;

    my @rpt;  # lines to write

    my $r = sub { push @rpt, @_ };  # append lines
    my $h = sub { push @rpt, "", "=" x 60, "  $_[0]", "=" x 60 };

    $h->("MENDER RUN REPORT");
    $r->("Start:       $start_stamp");
    $r->("End:         $end_stamp");
    $r->("Elapsed:     $elapsed_str");
    $r->("Config:      $config_file");
    $r->("Command:     $0 " . join(" ", @argv_orig));
    $r->("Dry run:     " . ($dry_run ? "YES" : "no"));

    $h->("INPUT");
    $r->("gff:             $gff");
    $r->("genome_fa:       $genome_fa");
    $r->("proteome_fa:     $proteome_fa");
    $r->("subject_fa:      $subject_fa");
    $r->("isoseq_gff:      " . ($isoseq_gff || "(not provided)"));

    $h->("OUTPUT");
    $r->("workdir:         $workdir");
    $r->("output_gff:      $output_gff");
    $r->("removed_gff:     $removed_gff");
    $r->("validated_gff:   $validated_gff  (" .
         ($tv_final_incl_review =~ /^no$/i
             ? "FAIL + REVIEW replaced with source genes"
             : "FAIL replaced with source genes; REVIEW retained") .
         ")");

    $h->("PARAMETERS");
    $r->("");
    $r->("[ids]");
    $r->("  gene_template:              $gene_template");
    $r->("  trans_template:             $trans_template");
    $r->("");
    $r->("[diamond]");
    $r->("  threads:                    $threads");
    $r->("  evalue:                     $evalue");
    $r->("");
    $r->("[merge_filters]");
    $r->("  skip_flags:                 $skip_flags");
    $r->("  flags:                      " . cfg("merge_filters", "flags", "all"));
    $r->("  fix_partial:                $fix_partial");
    $r->("  min_tiling:                 " . cfg("merge_filters", "min_tiling", "1"));
    $r->("  min_cov:                    " . cfg("merge_filters", "min_cov", "0"));
    $r->("  require_isoseq:             " . (cfg("merge_filters", "require_isoseq") || "(not set)"));
    $r->("  isoseq_min_spanning:        " . cfg("merge_filters", "isoseq_min_spanning", "0"));
    $r->("  large_span_warn:            ${large_span_warn} bp");
    $r->("  large_span_hard:            ${large_span_hard} bp");
    $r->("  asym_trim:                  " . cfg("merge_filters", "asym_trim", "yes"));
    $r->("");
    $r->("[validation]");
    $r->("  run_translation_validation: " . ($run_transl_val =~ /^no$/i ? "no" : "yes"));
    $r->("  min_junction_score:         $tv_min_junction");
    $r->("  min_merged_cov:             $tv_min_mrgd_cov");
    $r->("  min_ref_cov:                $tv_min_ref_cov");
    $r->("  aligner:                    $tv_aligner");
    $r->("  no_msa:                     " . ($tv_no_msa =~ /^yes$/i ? "yes" : "no"));
    $r->("  keep_msa:                   " . ($tv_keep_msa =~ /^yes$/i ? "yes" : "no"));
    $r->("  swissprot_fa:               " . ($tv_swissprot_fa || "(not set — using ref_fa)"));
    $r->("  transl_out_prefix:          " . ($tv_out_prefix || "(default: $workdir/transl)"));
    $r->("  final_gff_include_review:   " . ($tv_final_incl_review =~ /^no$/i ? "no" : "yes"));
    $r->("  run_agat:                   " . ($run_agat =~ /^no$/i ? "no" : "yes"));

    $h->("STEPS");
    my %step_names = (
        1 => "prepare",  2 => "diamond", 3 => "bedtools", 4 => "find",
        5 => "isoseq",   6 => "merge",   7 => "gt",       8 => "transl",
        9 => "agat",
    );
    for my $s (1..9) {
        my $name   = $step_names{$s};
        my $status = $run_steps{$s} ? "RUN" : "SKIPPED";
        my $reason = $step_skip{$s} ? "  ($step_skip{$s})" : "";
        $r->(sprintf "  %d  %-10s %s%s", $s, $name, $status, $reason);
    }

    # Results from TSV report
    my $report_tsv = "${tv_out}_report.tsv";
    if (-f $report_tsv) {
        $h->("RESULTS");
        my (%counts, $total);
        open my $fh, "<", $report_tsv or die;
        <$fh>;  # header
        while (<$fh>) {
            chomp;
            my @f = split /\t/;
            next unless @f >= 18;
            $counts{$f[17]}++;
            $total++;
        }
        close $fh;
        $total //= 0;
        $r->(sprintf "Merges processed:  %d", $total);
        for my $res (qw(PASS REVIEW FAIL)) {
            my $n = $counts{$res} // 0;
            my $pct = $total ? sprintf("%.1f%%", 100 * $n / $total) : "n/a";
            $r->(sprintf "  %-8s  %4d  (%s)", $res, $n, $pct);
        }
    }

    # Output file inventory
    $h->("OUTPUT FILES");
    my @out_files = (
        [$output_gff,    "complete merged annotation (all results)"],
        [$validated_gff, "validated annotation (" .
                         ($tv_final_incl_review =~ /^no$/i
                             ? "FAIL + REVIEW replaced with source genes"
                             : "FAIL replaced with source genes; REVIEW retained") .
                         ")"],
        [$removed_gff,   "source genes removed during merge (audit trail)"],
        ["${tv_out}_report.tsv",          "per-merge validation results"],
        ["${tv_out}_pass.gff3",           "PASS merges only"],
        ["${tv_out}_review.gff3",         "REVIEW merges (borderline MSA score)"],
        ["${tv_out}_fail.gff3",           "FAIL merges only"],
        ["${tv_out}_pass_proteins.fa",    "translated proteins for PASS merges"],
        ["$workdir/merge_candidates.txt", "merge candidates from step 4"],
        ["$workdir/isoseq_validated.txt", "IsoSeq-validated merge table (step 5)"],
    );
    for my $entry (@out_files) {
        my ($file, $desc) = @$entry;
        if (-f $file) {
            my $size = -s $file;
            my $hsize = $size > 1_000_000 ? sprintf("%.0fM", $size/1_000_000)
                      : $size > 1_000     ? sprintf("%.0fK", $size/1_000)
                      :                     "${size}B";
            $r->(sprintf "  %-50s  %6s  %s", $file, $hsize, $desc);
        } else {
            $r->(sprintf "  %-50s  %6s  %s", $file, "--", $desc);
        }
    }

    $r->("", "=" x 60, "");

    unless ($dry_run) {
        open my $rfh, ">", $report_file or die "Cannot write run report: $!\n";
        print $rfh "$_\n" for @rpt;
        close $rfh;
        print "\nRun report written: $report_file\n";
    }
}

print "\n", "=" x 60, "\n";
print $dry_run ? "DRY RUN COMPLETE\n" : "MENDER PIPELINE COMPLETE\n";
print "Output GFF (all merges):     $output_gff\n" unless $dry_run;
{
    my $tv_out = $tv_out_prefix || "$workdir/transl";
    (my $validated_gff = $output_gff) =~ s/(\.[^.]+)$/_validated$1/;
    print "Output GFF (validated only): $validated_gff\n"
        if -f $validated_gff && !$dry_run;
}
print "=" x 60, "\n";
