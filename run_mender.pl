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
##   10 fasta      export protein/CDS/transcript FASTAs from final validated GFF
##
## If --steps is not specified, all applicable steps are run in order.
## Steps 3 and 5 are skipped automatically if isoseq_gff is not set.
## Step 8 is skipped automatically if run_translation_validation = no in config.
## Step 9 is skipped automatically if run_agat = no in config.

use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use Cwd qw(abs_path);

# numbers ascending (a named comparator: sort passes the two elements in @_)
sub by_number($$) { my ($left, $right) = @_; return $left <=> $right }

# merges whose translation result is REVIEW
sub count_review {
    my ($attrs) = @_;
    my $n = 0;
    for my $gid (keys %$attrs) { $n++ if $attrs->{$gid}{result} eq "REVIEW" }
    return $n;
}

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
my $workdir        = cfg("output", "workdir",      "mender_workdir");
my $db_dir         = cfg("output", "db_dir",       "mender_db");
my $results_dir    = cfg("output", "results_dir",  "results");
my $run_prefix     = cfg("output", "run_prefix",   "")
    or die "ERROR: run_prefix is required in the [output] config section.\n" .
           "  Example: run_prefix = 20260414\n";
die "ERROR: run_prefix must not contain '/' or spaces: '$run_prefix'\n"
    if $run_prefix =~ m{[/\s]};
my $run_dir        = "$results_dir/$run_prefix";
my $run_workdir    = "$run_dir/work";
my $output_gff     = "$run_dir/merges.gff";
my $removed_gff    = "$run_dir/removed.gff";

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
my $skip_flags         = cfg("merge_filters", "skip_flags",         "SKIPPED_GENE,LOW_COV,PIECES_OVERLAP");
my $flags              = cfg("merge_filters", "flags",              "all");
my $min_tiling         = cfg("merge_filters", "min_tiling",         "1");
my $min_cov            = cfg("merge_filters", "min_cov",            "0");
my $require_isoseq     = cfg("merge_filters", "require_isoseq",     "");
my $isoseq_min_span    = cfg("merge_filters", "isoseq_min_spanning","0");
my $spanning_rescue    = cfg("merge_filters", "spanning_rescue",    "no");
my $asym_trim          = cfg("merge_filters", "asym_trim",          "yes");
my $large_span_warn    = cfg("merge_filters", "large_span_warn",    "500000");
my $large_span_hard    = cfg("merge_filters", "large_span_hard",    "2000000");
my $max_dist           = cfg("merge_filters", "max_dist",           "4");
my $wiggle             = cfg("merge_filters", "wiggle",             "15");
my $asym_threshold     = cfg("merge_filters", "asym_threshold",     "6");
my $low_cov_thresh     = cfg("merge_filters", "low_cov_thresh",     "60");
my $frag_min_qcov      = cfg("merge_filters", "frag_min_qcov",      "60");
my $end_tol            = cfg("merge_filters", "end_tol",            "30");
my $overlap_min_frac   = cfg("merge_filters", "overlap_min_frac",   "0.5");
my $overlap_min_qcov   = cfg("merge_filters", "overlap_min_qcov",   "60");

# validation (gt + agat)
my $run_gt    = cfg("validation", "run_gt",   "yes");
my $run_agat  = cfg("validation", "run_agat", "yes");

# translation validation (step 8)
my $run_transl_val  = cfg("validation", "run_translation_validation", "yes");
my $tv_min_mrgd_cov = cfg("validation", "min_merged_cov",             "0.80");
my $tv_min_ref_cov  = cfg("validation", "min_ref_cov",                "0.50");
my $tv_final_incl_review = cfg("validation", "final_gff_include_review", "yes");

# miniprot rebuild of merged genes (step 6): the tiling reference protein aligned
# across each merged locus; supported merges get the miniprot model
my $run_miniprot      = cfg("miniprot", "run_miniprot", "yes");
my $mp_flank          = cfg("miniprot", "flank",        "5000");
my $mp_min_prot_cov   = cfg("miniprot", "min_prot_cov", "0.80");
my $mp_max_refs       = cfg("miniprot", "max_refs",     "10");

# config keys of the MSA junction scoring, removed 2026-09-28 (it did not separate
# false merges from real ones): accepted but ignored, with a warning
my @obsolete_msa_keys;
for my $key (qw(min_junction_score no_msa keep_msa aligner w_conservation w_continuity w_gap
                min_msa_refs max_msa_refs swissprot_fa swissprot_db)) {
    push @obsolete_msa_keys, $key if exists $cfg{validation}{$key};
}
for my $key (qw(mafft_bin kalign_bin)) {
    push @obsolete_msa_keys, $key if exists $cfg{paths}{$key};
}

# paths
my $scripts_dir  = cfg("paths", "scripts_dir",  "") || dirname(abs_path($0));
my $diamond_bin  = cfg("paths", "diamond_bin",  "") || "diamond";
my $bedtools_bin = cfg("paths", "bedtools_bin", "") || "bedtools";
my $gt_bin       = cfg("paths", "gt_bin",       "") || "gt";
my $agat_bin     = cfg("paths", "agat_bin",     "") || "agat_convert_sp_gxf2gxf.pl";
my $gffread_bin  = cfg("paths", "gffread_bin",  "") || "gffread";
my $miniprot_bin = cfg("paths", "miniprot_bin", "") || "miniprot";
my $samtools_bin = cfg("paths", "samtools_bin", "") || "samtools";
my $perl_bin     = cfg("paths", "perl_bin",     "") || "perl";

# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------
{
    my (@errors, @warnings);

    # yes/no keys
    my %yesno = (
        fix_partial                => $fix_partial,
        asym_trim                  => $asym_trim,
        run_gt                     => $run_gt,
        run_agat                   => $run_agat,
        run_translation_validation => $run_transl_val,
        run_miniprot               => $run_miniprot,
        final_gff_include_review   => $tv_final_incl_review,
    );
    for my $key (sort keys %yesno) {
        push @errors, "$key: must be 'yes' or 'no', got '$yesno{$key}'"
            unless $yesno{$key} =~ /^(yes|no)$/i;
    }

    # enum keys
    if ($require_isoseq ne "") {
        push @errors, "require_isoseq: must be FULL_SPAN, PARTIAL_SPAN, or NO_SPANNERS (or blank), got '$require_isoseq'"
            unless $require_isoseq =~ /^(FULL_SPAN|PARTIAL_SPAN|NO_SPANNERS)$/;
    }
    push @warnings, "obsolete MSA key(s) ignored (the MSA junction scoring was removed): "
        . join(", ", @obsolete_msa_keys) if @obsolete_msa_keys;

    # numeric keys
    my %numeric = (
        threads              => $threads,
        min_tiling           => $min_tiling,
        min_cov              => $min_cov,
        max_dist             => $max_dist,
        wiggle               => $wiggle,
        asym_threshold       => $asym_threshold,
        low_cov_thresh       => $low_cov_thresh,
        frag_min_qcov        => $frag_min_qcov,
        end_tol              => $end_tol,
        overlap_min_qcov     => $overlap_min_qcov,
        isoseq_min_spanning  => $isoseq_min_span,
        flank                => $mp_flank,
        max_refs             => $mp_max_refs,
        large_span_warn      => $large_span_warn,
        large_span_hard      => $large_span_hard,
    );
    for my $key (sort keys %numeric) {
        push @errors, "$key: must be numeric, got '$numeric{$key}'"
            unless $numeric{$key} =~ /^-?\d+(\.\d+)?$/;
    }

    # 0-1 range
    my %zeroone = (
        min_merged_cov     => $tv_min_mrgd_cov,
        min_ref_cov        => $tv_min_ref_cov,
        min_prot_cov       => $mp_min_prot_cov,
        overlap_min_frac   => $overlap_min_frac,
    );
    for my $key (sort keys %zeroone) {
        my $val = $zeroone{$key};
        if ($val !~ /^\d+(\.\d+)?$/) {
            push @errors, "$key: must be numeric (0-1), got '$val'";
        } elsif ($val < 0 || $val > 1) {
            push @errors, "$key: must be between 0 and 1, got '$val'";
        }
    }

    # flag name validation in skip_flags / flags (warn only — unknown flags are skipped silently
    # at runtime, so this catches typos before the run wastes time)
    my %known_flags;
    for my $known (qw(
        CLEAN STRONG SINGLE_HIT WEAK_END WEAK_INTERNAL LOW_COV
        SKIPPED_GENE OPPOSITE_STRAND_SKIP TRANSITIVE_JOIN MULTI_ISOFORM_JOIN LARGE_SPAN LARGE_SPAN_EXTREME
        JUNCTION_OVERHANG PIECES_OVERLAP
    )) { $known_flags{$known} = 1 }
    for my $item (["skip_flags", $skip_flags], ["flags", $flags]) {
        my ($key, $val) = @$item;
        next if $val eq "all" || $val eq "";
        for my $flag (split /,/, $val) {
            $flag =~ s/^\s+|\s+$//g;
            push @warnings, "$key: unknown flag '$flag'"
                unless $known_flags{$flag};
        }
    }

    if (@warnings) {
        print "WARNING: config issue in $config_file:\n";
        for my $warning (@warnings) { print "  $warning\n" }
        print "\n";
    }
    if (@errors) {
        print "ERROR: invalid config in $config_file:\n";
        for my $error (@errors) { print "  $error\n" }
        die "\nAborting — fix the config errors above before running.\n";
    }
}

# intermediate file paths
my $clean_fa         = "$workdir/protein.nostops.fa";
my $diamond_db       = "$db_dir/subject.dmnd";
my $diamond_out      = "$workdir/diamond.out";
my $genes_gff        = "$workdir/models.gene.gff";
my $isoseq_mrna_gff  = "$workdir/isoseq.mrna.gff";
my $overlaps         = "$workdir/all.overlaps";
my $merge_table      = "$run_workdir/merge_candidates.txt";
my $validated_table  = "$run_workdir/isoseq_validated.txt";
my $find_script      = "$scripts_dir/find_split_genes.pl";
my $isoseq_script    = "$scripts_dir/validate_merge_with_isoseq.pl";
my $merge_script     = "$scripts_dir/merge_split_genes.pl";
my $transl_script    = "$scripts_dir/validate_merge_translation.pl";
my $miniprot_script  = "$scripts_dir/rebuild_with_miniprot.pl";

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
    @steps_requested = sort by_number keys %run_steps;
    # steps not in --steps are "not requested"
    for my $s (1..10) {
        $step_skip{$s} = "not in --steps" unless $run_steps{$s};
    }
} else {
    for my $s (1 .. 10) { $run_steps{$s} = 1 }
    @steps_requested = (1..10);
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

# Reuse events collected here for the run report
my @reuse_log;

# shared_step — coordinate concurrent runs sharing the same output files.
#
# $files    : arrayref of output file paths; first file is used as lock key
# $desc     : human-readable label (shown on screen and in run report)
# $generate : coderef that creates the output files
#
# Behaviour:
#   - If lock file exists and is < 6 hours old: wait until it clears, then reuse
#   - If lock file exists and is >= 6 hours old: treat as stale, remove and regenerate
#   - If all output files already exist (and no lock): print reuse message and return
#   - Otherwise: create lock, run $generate, remove lock
#
sub shared_step {
    my ($files, $desc, $generate) = @_;
    my $lock = $files->[0] . ".lock";

    # Stale lock check (6 hours)
    if (-f $lock) {
        my $age = time() - (stat($lock))[9];
        if ($age > 21600) {
            print "WARNING: stale lock for '$desc' (${age}s old) — removing and regenerating\n";
            unlink $lock unless $dry_run;
        }
    }

    # Wait if another process is currently generating
    if (-f $lock) {
        my $ts = strftime("%Y-%m-%d %H:%M:%S", localtime((stat($lock))[9]));
        print "\nWAIT: '$desc' is being generated by another process (lock created $ts)...\n";
        my $waited = 0;
        while (-f $lock) {
            sleep 15;
            $waited += 15;
            print "  ... still waiting for '$desc' (${waited}s elapsed)\n"
                if $waited % 60 == 0;
        }
        print "  Lock cleared after ${waited}s\n";
    }

    # Reuse if all output files already exist
    my @missing;
    for my $file (@$files) { push @missing, $file unless -f $file }
    unless (@missing) {
        print "\n--- $desc ---\n";
        print "REUSE: output already exists — skipping regeneration\n";
        for my $f (@$files) {
            my $mts = strftime("%Y-%m-%d %H:%M:%S", localtime((stat($f))[9]));
            print "  $f  (modified $mts)\n";
            push @reuse_log, "REUSE  $desc: $f (modified $mts)";
        }
        return;
    }

    # Create lock, run, remove lock
    unless ($dry_run) {
        open my $lk, ">", $lock or die "Cannot create lock $lock: $!\n";
        print $lk "pid=$$\nstarted=" . time() . "\n";
        close $lk;
    }
    my $ok = eval { $generate->(); 1 };
    my $err = $@;
    unlink $lock unless $dry_run;
    die $err if $err;
}

# ---------------------------------------------------------------------------
# Print run summary
# ---------------------------------------------------------------------------
print "=" x 60, "\n";
print "MENDER PIPELINE\n";
print "=" x 60, "\n";
print "Config:          $config_file\n";
print "Date:            $date_tag\n";
print "Steps:           ", join(", ", sort by_number keys %run_steps), "\n";
print "Dry run:         ", ($dry_run ? "YES" : "no"), "\n";
print "-" x 60, "\n";
print "Input GFF:       $gff\n";
print "Proteome FA:     $proteome_fa\n";
print "Subject FA:      $subject_fa\n";
print "IsoSeq GFF:      ", ($isoseq_gff || "(not provided)"), "\n";
print "Run dir:         $run_dir\n";
print "Run workdir:     $run_workdir\n";
print "Shared workdir:  $workdir\n";
print "DB dir:          $db_dir\n";
print "Gene template:   $gene_template\n";
print "Trans template:  $trans_template\n";
print "Skip flags:      $skip_flags\n";
print "Spanning rescue: $spanning_rescue\n";
print "Fix partial:     $fix_partial\n";
print "Large span warn: $large_span_warn bp\n";
print "Large span hard: $large_span_hard bp\n";
print "Transl val:      ", ($run_transl_val =~ /^no$/i ? "no" : "yes"), "\n";
if ($run_transl_val !~ /^no$/i) {
    print "  min_merged_cov:$tv_min_mrgd_cov\n";
    print "  min_ref_cov:   $tv_min_ref_cov\n";
}
print "miniprot rebuild: ", ($run_miniprot =~ /^no$/i ? "no" : "yes (flank $mp_flank bp, min_prot_cov $mp_min_prot_cov, up to $mp_max_refs references)"), "\n";
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

    # gffread: needed for translation validation (step 8) or gene set export (step 10)
    push @checks, [$gffread_bin,  "gffread",  $run_transl_val !~ /^no$/i || $run_steps{10}];

    # miniprot + samtools: the rebuild of merged genes in step 6 (run_miniprot = yes)
    my $need_mp = $run_miniprot !~ /^no$/i && $run_steps{6};
    push @checks, [$miniprot_bin, "miniprot", $need_mp];
    push @checks, [$samtools_bin, "samtools", $need_mp];

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
        for my $tool (@missing) { print "  $tool\n" }
        print "Load the relevant modules or set bin paths in the [paths] config section.\n";
        exit 1;
    }
    print "  All required tools found.\n";

    # miniprot reads the merged locus through a samtools faidx index of the genome
    if ($need_mp && !-e "$genome_fa.fai") {
        print "ERROR: run_miniprot = yes needs a genome index $genome_fa.fai (samtools faidx $genome_fa)\n";
        exit 1;
    }
}

# create workdir, db dir, results dir, run dir, and run workdir
unless ($dry_run) {
    mkdir $workdir unless -d $workdir;
    die "ERROR: could not create workdir: $workdir\n" unless -d $workdir;
    mkdir $db_dir unless -d $db_dir;
    die "ERROR: could not create db_dir: $db_dir\n" unless -d $db_dir;
    mkdir $results_dir unless -d $results_dir;
    die "ERROR: could not create results_dir: $results_dir\n" unless -d $results_dir;
    mkdir $run_dir unless -d $run_dir;
    die "ERROR: could not create run_dir: $run_dir\n" unless -d $run_dir;
    mkdir $run_workdir unless -d $run_workdir;
    die "ERROR: could not create run_workdir: $run_workdir\n" unless -d $run_workdir;
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

    my @step1_files = ($clean_fa, $genes_gff);
    push @step1_files, $isoseq_mrna_gff if $isoseq_gff;

    shared_step(\@step1_files, "Step 1 prepare", sub {
        # clean protein fasta for diamond:
        #   1. read full FASTA records (multi-line safe via $/ trick)
        #   2. strip terminal stop (* or .) — biologically normal, just upsets diamond
        #   3. remove entire sequence if internal stops or dots remain
        run_cmd(
            "perl -e 'local \$/ = \">\"; while (my \$rec = <>) { chomp \$rec; next unless \$rec =~ /\\S/; " .
            "my (\$h, \$s) = split(/\\n/, \$rec, 2); \$s =~ s/\\n//g; " .
            "\$s =~ s/[*.]\$//; " .
            "print \">\$h\\n\$s\\n\" unless \$s =~ /[*.]/ }' $proteome_fa > $clean_fa",
            "Clean protein fasta: strip terminal stops, remove sequences with internal stops"
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
    });
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

    shared_step([$diamond_db], "Build diamond database", sub {
        run_cmd(
            "$diamond_bin makedb --in $subject_fa --db $diamond_db",
            "Build diamond database"
        );
    });

    shared_step([$diamond_out], "Run diamond blastp", sub {
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
    });
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

    shared_step([$overlaps], "Bedtools intersect IsoSeq overlaps", sub {
        run_cmd(
            "$bedtools_bin intersect " .
            "-a $isoseq_mrna_gff -b $genes_gff -wa -wb | " .
            "$perl_bin -p -e " .
            "'s/.+\\sID=([^;]+).+\\sID=([^;]+).*\$/\$1\\t\$2/' " .
            "> $overlaps",
            "Bedtools intersect IsoSeq reads vs gene models"
        );
    });
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
        ($max_dist         ? "--max_dist $max_dist "                      : "") .
        ($wiggle           ? "--wiggle $wiggle "                          : "") .
        ($asym_threshold   ? "--asym_threshold $asym_threshold "         : "") .
        ($low_cov_thresh   ? "--low_cov_thresh $low_cov_thresh "         : "") .
        "--frag_min_qcov $frag_min_qcov --end_tol $end_tol --max_evalue $evalue " .
        "--overlap_min_frac $overlap_min_frac --overlap_min_qcov $overlap_min_qcov " .
        "$diamond_out " .
        "$gff " .
        "$subject_fa " .
        "$clean_fa",
        "Find split gene candidates"
    );

    # move outputs to workdir
    run_cmd("mv split_genes_summary.txt split_genes_detail.txt $run_workdir/ && " .
            "mv split_genes_merge.txt $run_workdir/merge_candidates.txt",
            "Move find_split_genes outputs to run workdir");
    $merge_table = "$run_workdir/merge_candidates.txt";
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
    $merge_cmd .= " --spanning_rescue"          if $spanning_rescue eq "yes";
    $merge_cmd .= " --flags $flags"             if $flags && $flags ne "all";
    $merge_cmd .= " --min_tiling $min_tiling"   if $min_tiling > 1;
    $merge_cmd .= " --min_cov $min_cov"         if $min_cov > 0;
    $merge_cmd .= " --require_isoseq $require_isoseq"  if $require_isoseq;
    $merge_cmd .= " --isoseq_min_spanning $isoseq_min_span" if $isoseq_min_span > 0;
    $merge_cmd .= " --removed $removed_gff";
    # the genome lets the merge trim a fragment's last codon only when it is a stop
    if ($genome_fa && -e "$genome_fa.fai") {
        $merge_cmd .= " --genome_fa $genome_fa";
    } else {
        print "WARNING: genome_fa or its .fai index (samtools faidx) missing: fragment stop codons are trimmed unchecked\n";
    }
    $merge_cmd .= " --gene_template  \"$gene_template\"";
    $merge_cmd .= " --trans_template \"$trans_template\"";
    $merge_cmd .= " $validated_table $gff $output_gff";

    run_cmd($merge_cmd, "Merge split genes into annotation");

    # Rebuild each merged gene from protein evidence: the reference protein the
    # fragments tiled against, aligned across the merged locus with miniprot. A
    # supported merge gets the miniprot model (real splice sites, missing exons
    # recovered); an unsupported one keeps the placeholder and fails step 8
    # (NO_PROTEIN_SPAN). The placeholder models are kept in merges.placeholder.gff.
    if ($run_miniprot !~ /^no$/i) {
        check_file($miniprot_script, "rebuild_with_miniprot.pl");
        my $placeholder_gff = "$run_dir/merges.placeholder.gff";
        run_cmd("mv $output_gff $placeholder_gff", "Keep the placeholder merged models");
        run_cmd("$perl_bin $miniprot_script"
              . " --merged_gff $placeholder_gff --removed_gff $removed_gff"
              . " --merge_table $validated_table --ref_fa $subject_fa --genome_fa $genome_fa"
              . " --out $output_gff --report $run_dir/miniprot.tsv --scratch $run_workdir"
              . " --flank $mp_flank --min_prot_cov $mp_min_prot_cov --max_refs $mp_max_refs --threads $threads"
              . " --miniprot_bin $miniprot_bin --samtools_bin $samtools_bin",
              "Rebuild merged genes with miniprot");
    }
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
        my $mender_only_gff = "$run_workdir/mender_only.gff";
        print "\n--- GT GFF3 validation (Mender genes only) ---\n";
        unless ($dry_run) {
            open my $mo, ">", $mender_only_gff
                or die "Cannot write $mender_only_gff: $!\n";
            print $mo "##gff-version 3\n";
            open my $fh, "<", $output_gff
                or die "Cannot open $output_gff: $!\n";
            while (my $line = <$fh>) { print $mo $line if $line =~ /\tMender\t/ }
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

    my $tv_out = $run_dir;

    my $tv_cmd = "$perl_bin $transl_script"
        . " --merged_gff $output_gff"
        . " --genome_fa $genome_fa"
        . " --ref_fa $subject_fa"
        . " --diamond_db $diamond_db"
        . " --merge_table $validated_table"
        . " --out $tv_out"
        . " --scratch $run_workdir"
        . " --threads $threads"
        . " --min_merged_cov $tv_min_mrgd_cov"
        . " --min_ref_cov $tv_min_ref_cov"
        . " --gffread_bin $gffread_bin"
        . " --diamond_bin $diamond_bin";

    run_cmd($tv_cmd, "Translation validation of merged genes");

    unless ($dry_run) {
        print "\n--- Translation validation summary ---\n";
        my $pass_gff   = "$tv_out/pass.gff3";
        my $report_tsv = "$tv_out/report.tsv";
        if (-f $report_tsv) {
            my $pass   = `grep -c 'PASS'   $report_tsv 2>/dev/null` || 0;
            my $fail   = `grep -c 'FAIL'   $report_tsv 2>/dev/null` || 0;
            my $review = `grep -c 'REVIEW' $report_tsv 2>/dev/null` || 0;
            chomp ($pass, $fail, $review);
            # subtract header line matches if any
            print "  Report:   $report_tsv\n";
            print "  Pass GFF: $tv_out/pass.gff3\n";
            print "  Fail GFF: $tv_out/fail.gff3\n";
        }
    }
}

# ---------------------------------------------------------------------------
# STEP 9: AGAT GENE-MODEL CHECK (on pass GFF only)
# ---------------------------------------------------------------------------
if ($run_steps{9}) {
    print "\n", "=" x 60, "\n";
    print "STEP 9: AGAT GENE-MODEL CHECK\n";
    print "=" x 60, "\n";

    my $tv_out = $run_dir;
    my $pass_gff  = "$tv_out/pass.gff3";

    # fall back to full merged GFF if translation validation was not run
    my $target_gff = (-f $pass_gff) ? $pass_gff : $output_gff;
    if ($target_gff eq $output_gff) {
        print "INFO: pass GFF not found — running AGAT on full merged GFF\n";
        print "      (run step 8 first to validate on pass merges only)\n";
    } else {
        # warn if pass.gff3 is older than merges.gff — step 8 may be stale
        if (-f $output_gff) {
            my $pass_mtime   = (stat($pass_gff))[9];
            my $merges_mtime = (stat($output_gff))[9];
            if ($pass_mtime < $merges_mtime) {
                my $pass_ts   = strftime("%Y-%m-%d %H:%M:%S", localtime($pass_mtime));
                my $merges_ts = strftime("%Y-%m-%d %H:%M:%S", localtime($merges_mtime));
                print "WARNING: pass.gff3 ($pass_ts) is older than merges.gff ($merges_ts)\n";
                print "         Step 8 may need to be re-run before step 9.\n";
                print "         Re-run with: --steps 8,9\n";
            }
        }
        print "  Running AGAT on: $pass_gff\n";
    }

    check_file($target_gff, "GFF for AGAT validation");

    (my $agat_out = $target_gff) =~ s/(\.[^.]+)$/_agat$1/;
    print "  AGAT corrected GFF: $agat_out\n";

    run_cmd(
        "$agat_bin -g $target_gff -o $agat_out 2>&1 | " .
        "grep -i 'error\\|warn' | head -20",
        "Gene-model check (AGAT)"
    );
}

# ---------------------------------------------------------------------------
# FINAL VALIDATED GFF
# ---------------------------------------------------------------------------
# Builds a complete annotation GFF where FAIL merges (and optionally REVIEW)
# are replaced with their original source genes. Adds the transl_result
# attribute to all retained Mender gene features.
# Generated automatically whenever translation validation outputs are present.
# ---------------------------------------------------------------------------
{
    my $tv_out     = $run_dir;
    my $fail_gff   = "$tv_out/fail.gff3";
    my $report_tsv = "$tv_out/report.tsv";

    if (-f $fail_gff && !$dry_run) {
        print "\n", "=" x 60, "\n";
        print "FINAL VALIDATED GFF\n";
        print "=" x 60, "\n";

        # Load validation results from the TSV report, by column NAME (a fixed column
        # number silently read the wrong column when the report's columns changed)
        my %transl_attrs;  # new_gene_id -> hashref
        if (-f $report_tsv) {
            open my $fh, "<", $report_tsv or die "Cannot open $report_tsv: $!\n";
            chomp(my $hdr = <$fh> // "");
            my @cols = split /\t/, $hdr;
            my %ci;
            for my $i (0 .. $#cols) { $ci{ $cols[$i] } = $i }
            die "ERROR: $report_tsv lacks new_gene_id / overall_result columns\n"
                unless defined $ci{new_gene_id} && defined $ci{overall_result};
            while (my $line = <$fh>) {
                chomp $line;
                my @f = split /\t/, $line;
                $transl_attrs{ $f[$ci{new_gene_id}] } = { result => $f[$ci{overall_result}] };
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
            while (my $line = <$fh>) {
                next unless $line =~ /\tgene\t/;
                chomp $line;
                my $attr_str = (split /\t/, $line)[8] // "";
                my ($gid)   = $attr_str =~ /\bID=([^;]+)/;
                my ($mfrom) = $attr_str =~ /\bmerged_from=([^;]+)/;
                next unless $gid && $mfrom;
                $exclude{$gid} = [split /,/, $mfrom];
            }
            close $fh;
        };

        $parse_excl_gff->($fail_gff);
        $n_fail_excl = scalar keys %exclude;

        my $review_gff = "$tv_out/review.gff3";
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
        my $validated_gff = "$run_dir/validated.gff";
        {
            open my $in,  "<", $output_gff   or die "Cannot open $output_gff: $!\n";
            open my $out, ">", $validated_gff or die "Cannot write $validated_gff: $!\n";

            my ($cur_id, @cur_block, $cur_is_mender, $cur_is_excluded);

            my $flush = sub {
                return unless @cur_block;
                if ($cur_is_excluded) {
                    # Replace excluded Mender gene with its original source gene blocks
                    for my $src_id (@{ $exclude{$cur_id} }) {
                        for my $line (@{ $removed_blocks{$src_id} // [] }) { print $out "$line\n" }
                    }
                } else {
                    # Annotate retained Mender genes with transl validation flags
                    if ($cur_is_mender && exists $transl_attrs{$cur_id}) {
                        my $ta = $transl_attrs{$cur_id};
                        $cur_block[0] .= ";transl_result=$ta->{result}";
                    }
                    for my $line (@cur_block) { print $out "$line\n" }
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
            (count_review(\%transl_attrs) - $n_review_excl)
            if %transl_attrs;
        print  "  Final validated GFF: $validated_gff\n";
    }
}

# ---------------------------------------------------------------------------
# STEP 10: EXPORT THE GENE SET (datastore layout) FROM THE FINAL VALIDATED GFF
# ---------------------------------------------------------------------------
if ($run_steps{10}) {
    print "\n", "=" x 60, "\n";
    print "STEP 10: EXPORT THE GENE SET (GFF, GTF, FASTAs, id maps)\n";
    print "=" x 60, "\n";

    check_file($genome_fa, "genome_fa");

    # Use validated.gff if present (built above), otherwise fall back to merges.gff
    my $validated_gff = "$run_dir/validated.gff";
    my $fasta_src_gff = (-f $validated_gff) ? $validated_gff : $output_gff;
    if ($fasta_src_gff eq $output_gff) {
        print "INFO: validated.gff not found — exporting from full merged GFF\n";
        print "      (run step 8 first to get a validated.gff)\n";
    } else {
        print "  Exporting FASTAs from: $fasta_src_gff\n";
    }

    check_file($fasta_src_gff, "GFF for FASTA export");

    # the gene set, laid out for the genome datastore (genes.gff, genes.gtf, protein.aa.fa,
    # cds.nt.fa, transcript.nt.fa, protein2gene.txt, cds2gene.txt, transcript2gene.txt);
    # descriptions: the input proteome's for untouched genes, merge provenance for Mender genes
    my $geneset_dir = "$run_dir/geneset";
    run_cmd(
        "$perl_bin $scripts_dir/export_geneset.pl --gff $fasta_src_gff --genome_fa $genome_fa " .
        "--proteome_fa $proteome_fa --out_dir $geneset_dir --gffread_bin $gffread_bin",
        "Export the gene set (GFF, GTF, FASTAs, id-to-gene maps) for the datastore"
    );
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

    my $tv_out        = $run_dir;
    my $validated_gff = "$run_dir/validated.gff";
    my $report_file   = "$run_dir/run_report.txt";

    my @rpt;  # lines to write

    my $r = sub { push @rpt, @_ };  # append lines
    my $h = sub { my ($title) = @_; push @rpt, "", "=" x 60, "  $title", "=" x 60 };

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
    $r->("shared_workdir:  $workdir");
    $r->("db_dir:          $db_dir");
    $r->("run_dir:         $run_dir");
    $r->("run_workdir:     $run_workdir");
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
    $r->("  spanning_rescue:            $spanning_rescue");
    $r->("  flags:                      " . cfg("merge_filters", "flags", "all"));
    $r->("  fix_partial:                $fix_partial");
    $r->("  min_tiling:                 " . cfg("merge_filters", "min_tiling", "1"));
    $r->("  min_cov:                    " . cfg("merge_filters", "min_cov", "0"));
    $r->("  require_isoseq:             " . (cfg("merge_filters", "require_isoseq") || "(not set)"));
    $r->("  isoseq_min_spanning:        " . cfg("merge_filters", "isoseq_min_spanning", "0"));
    $r->("  large_span_warn:            ${large_span_warn} bp");
    $r->("  large_span_hard:            ${large_span_hard} bp");
    $r->("  asym_trim:                  " . cfg("merge_filters", "asym_trim", "yes"));
    $r->("  max_dist:                   ${max_dist}");
    $r->("  wiggle:                     ${wiggle} aa");
    $r->("  asym_threshold:             ${asym_threshold}");
    $r->("  low_cov_thresh:             ${low_cov_thresh}%");
    $r->("  frag_min_qcov:              ${frag_min_qcov}%");
    $r->("  end_tol:                    ${end_tol} aa");
    $r->("  overlap_min_frac:           ${overlap_min_frac}");
    $r->("  overlap_min_qcov:           ${overlap_min_qcov}%");
    $r->("");
    $r->("[validation]");
    $r->("  run_translation_validation: " . ($run_transl_val =~ /^no$/i ? "no" : "yes"));
    $r->("  min_merged_cov:             $tv_min_mrgd_cov");
    $r->("  min_ref_cov:                $tv_min_ref_cov");
    $r->("");
    $r->("[miniprot]");
    $r->("  run_miniprot:               " . ($run_miniprot =~ /^no$/i ? "no" : "yes"));
    $r->("  flank:                      $mp_flank bp");
    $r->("  min_prot_cov:               $mp_min_prot_cov");
    $r->("  max_refs:                   $mp_max_refs");
    $r->("  run_dir:                    $run_dir");
    $r->("  run_workdir:                $run_workdir");
    $r->("  final_gff_include_review:   " . ($tv_final_incl_review =~ /^no$/i ? "no" : "yes"));
    $r->("  run_agat:                   " . ($run_agat =~ /^no$/i ? "no" : "yes"));

    $h->("STEPS");
    my %step_names = (
        1 => "prepare",  2 => "diamond", 3 => "bedtools", 4 => "find",
        5 => "isoseq",   6 => "merge",   7 => "gt",       8 => "transl",
        9 => "agat",    10 => "fasta",
    );
    for my $s (1..10) {
        my $name   = $step_names{$s};
        my $status = $run_steps{$s} ? "RUN" : "SKIPPED";
        my $reason = $step_skip{$s} ? "  ($step_skip{$s})" : "";
        $r->(sprintf "  %d  %-10s %s%s", $s, $name, $status, $reason);
    }

    # Results from TSV report
    my $report_tsv = "$tv_out/report.tsv";
    if (-f $report_tsv) {
        $h->("RESULTS");
        my (%counts, $total);
        open my $fh, "<", $report_tsv or die;
        chomp(my $hdr = <$fh> // "");
        my @cols = split /\t/, $hdr;
        my $c_res;
        for my $i (0 .. $#cols) { $c_res //= $i if $cols[$i] eq "overall_result" }
        while (my $line = <$fh>) {
            chomp $line;
            my @f = split /\t/, $line;
            next unless defined $c_res && defined $f[$c_res];
            $counts{$f[$c_res]}++;
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
        ["$tv_out/report.tsv",          "per-merge validation results"],
        ["$tv_out/pass.gff3",           "PASS merges only"],
        ["$tv_out/review.gff3",         "REVIEW merges (merged protein covers < min_ref_cov of its reference)"],
        ["$tv_out/fail.gff3",           "FAIL merges only"],
        ["$tv_out/miniprot.tsv",        "miniprot rebuild: per merge, SUPPORTED / NOT_SUPPORTED and why"],
        ["$tv_out/merges.placeholder.gff", "merged models before the miniprot rebuild"],
        ["$tv_out/pass_proteins.fa",       "translated proteins for PASS merges (step 8)"],
        ["$tv_out/geneset/genes.gff",        "gene set for the datastore: annotation (step 10)"],
        ["$tv_out/geneset/genes.gtf",        "gene set: annotation as GTF (step 10)"],
        ["$tv_out/geneset/protein.aa.fa",    "gene set: proteins (step 10)"],
        ["$tv_out/geneset/cds.nt.fa",        "gene set: CDS (step 10)"],
        ["$tv_out/geneset/transcript.nt.fa", "gene set: transcripts, UTRs included (step 10)"],
        ["$tv_out/geneset/protein2gene.txt", "gene set: protein id -> gene id (also cds2gene.txt, transcript2gene.txt)"],
        ["$run_workdir/merge_candidates.txt", "merge candidates from step 4"],
        ["$run_workdir/isoseq_validated.txt", "IsoSeq-validated merge table (step 5)"],
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

    if (@reuse_log) {
        $r->("", "Shared files reused from previous run:");
        for my $reused (@reuse_log) { $r->("  $reused") }
    }

    $r->("", "=" x 60, "");

    unless ($dry_run) {
        open my $rfh, ">", $report_file or die "Cannot write run report: $!\n";
        for my $line (@rpt) { print $rfh "$line\n" }
        close $rfh;
        print "\nRun report written: $report_file\n";
    }
}

print "\n", "=" x 60, "\n";
print $dry_run ? "DRY RUN COMPLETE\n" : "MENDER PIPELINE COMPLETE\n";
print "Output GFF (all merges):     $output_gff\n" unless $dry_run;
{
    my $validated_gff = "$run_dir/validated.gff";
    print "Output GFF (validated only): $validated_gff\n"
        if -f $validated_gff && !$dry_run;
}
print "=" x 60, "\n";
