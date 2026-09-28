#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

## rebuild_with_miniprot.pl
##
## Rebuilds each Mender-merged gene from protein evidence, and says whether that
## evidence supports the merge at all.
##
## merge_split_genes.pl joins the source fragments' own exons, so the joined
## "intron" runs from one fragment's last codon to the next fragment's start
## codon (no real splice sites) and any exon the original annotation missed at
## the break is still missing. Here the reference protein the fragments tiled
## against (best_hit in the merge table) is aligned to the merged locus with
## miniprot (spliced protein-to-genome alignment). A true split gives ONE
## alignment covering most of the reference protein and every source fragment,
## with real splice sites, no frameshift and no in-frame stop: that alignment
## becomes the merged gene's model. Two different genes joined by mistake give an
## alignment that covers only one of them.
##
## Per merged gene (Mender gene lines in the merged GFF):
##   region  = merged gene span +/- --flank bp (clipped to the sequence)
##   protein = best_hit of its merge_id in the merge table
##   miniprot --gff on (region, protein); the Rank=1 alignment is judged:
##     SUPPORTED     protein coverage >= --min_prot_cov, no frameshift (fs:i:0),
##                   no in-frame stop (st:i:0), and its CDS overlap every source gene
##     NOT_SUPPORTED otherwise (the reason is recorded)
##   SUPPORTED:     the gene's transcripts are replaced by the miniprot model: one
##                  mRNA, exon = CDS (no UTRs), phases from miniprot; gene span = model
##   NOT_SUPPORTED: the placeholder model is kept; validate_merge_translation.pl
##                  fails it (NO_PROTEIN_SPAN)
## Attributes on the gene (and mRNA): miniprot_support=SUPPORTED|NOT_SUPPORTED:<reason>,
## miniprot_ref, miniprot_prot_cov, miniprot_identity, miniprot_positive.
##
## Usage:
##   perl rebuild_with_miniprot.pl --merged_gff merges.gff --removed_gff removed.gff \
##        --merge_table merge_candidates.txt \
##        --ref_fa ref.fa --genome_fa genome.fa --out merges.miniprot.gff --report miniprot.tsv
##        [--flank 5000] [--min_prot_cov 0.80] [--miniprot_bin miniprot] [--samtools_bin samtools]
##        [--threads 1] [--scratch DIR]

my ($merged_gff, $merge_table, $ref_fa, $genome_fa, $out_gff, $report, $removed_gff);
my $flank        = 5000;
my $min_prot_cov = 0.80;
my $miniprot_bin = "miniprot";
my $samtools_bin = "samtools";
my $threads      = 1;
my $scratch      = ".";
GetOptions(
    "merged_gff=s"   => \$merged_gff,
    "removed_gff=s"  => \$removed_gff,
    "merge_table=s"  => \$merge_table,
    "ref_fa=s"       => \$ref_fa,
    "genome_fa=s"    => \$genome_fa,
    "out=s"          => \$out_gff,
    "report=s"       => \$report,
    "flank=i"        => \$flank,
    "min_prot_cov=f" => \$min_prot_cov,
    "miniprot_bin=s" => \$miniprot_bin,
    "samtools_bin=s" => \$samtools_bin,
    "threads=i"      => \$threads,
    "scratch=s"      => \$scratch,
) or die "Error in options\n";
for my $required ([merged_gff => $merged_gff], [removed_gff => $removed_gff], [merge_table => $merge_table], [ref_fa => $ref_fa],
                  [genome_fa => $genome_fa], [out => $out_gff], [report => $report]) {
    die "--$required->[0] is required\n" unless defined $required->[1];
}
die "genome index $genome_fa.fai missing (samtools faidx $genome_fa)\n" unless -e "$genome_fa.fai";

# sequence lengths, for clipping the flank
my %seq_len;
open my $fai, "<", "$genome_fa.fai" or die "cant open $genome_fa.fai: $!\n";
while (my $line = <$fai>) { my ($name, $len) = split /\t/, $line; $seq_len{$name} = $len }
close $fai;

# merge_id -> best_hit (the reference protein the fragments tiled against)
my %best_hit;
open my $mt, "<", $merge_table or die "cant open $merge_table: $!\n";
my @hdr = split /\t/, scalar <$mt>;
chomp @hdr;
my %col = map { $hdr[$_] => $_ } 0 .. $#hdr;
die "merge table needs merge_id and best_hit columns\n" unless defined $col{merge_id} && defined $col{best_hit};
while (my $line = <$mt>) {
    chomp $line;
    my @f = split /\t/, $line;
    $best_hit{$f[$col{merge_id}]} = $f[$col{best_hit}];
}
close $mt;

# reference proteins needed
my %wanted = map { $_ => 1 } grep { defined && $_ ne "" && $_ ne "?" } values %best_hit;
my (%ref_seq, $id);
open my $rf, "<", $ref_fa or die "cant open $ref_fa: $!\n";
while (my $line = <$rf>) {
    chomp $line;
    if ($line =~ /^>(\S+)/) { $id = $wanted{$1} ? $1 : undef; $ref_seq{$id} = "" if defined $id; next }
    $ref_seq{$id} .= $line if defined $id;
}
close $rf;

# read the merged GFF: Mender gene blocks (gene line, then its mRNA and children) and everything else
my (@order, %block, %gene_info);
my $current;
open my $in, "<", $merged_gff or die "cant open $merged_gff: $!\n";
while (my $line = <$in>) {
    chomp $line;
    my @f = split /\t/, $line;
    if (@f >= 9 && $f[2] eq "gene") {
        my ($gid) = $f[8] =~ /ID=([^;]+)/;
        if ($f[1] eq "Mender") {
            my ($mid) = $f[8] =~ /merge_id=([^;]+)/;
            my ($from) = $f[8] =~ /merged_from=([^;]+)/;
            $current = $gid;
            $gene_info{$gid} = { line => $line, chr => $f[0], start => $f[3], end => $f[4], strand => $f[6],
                                 merge_id => $mid, sources => [ split /,/, $from // "" ] };
            push @order, ["mender", $gid];
            next;
        }
        $current = undef;
    } elsif (@f >= 9 && defined $current) {
        push @{ $block{$current} }, $line;
        next;
    } elsif ($line =~ /^#/ || @f < 9) {
        $current = undef if $line !~ /^#/;
    }
    push @order, ["line", $line];
}
close $in;

# source gene spans: the source genes are gone from the merged GFF; merge_split_genes.pl
# wrote them (with their merge group) to its removed-genes GFF
my %source_span;   # source gene id -> [start, end]
my (%tx_gene, %tx_cds);   # source transcript -> gene; -> [[start, end, phase]]
open my $rm, "<", $removed_gff or die "cant open $removed_gff: $!\n";
while (my $line = <$rm>) {
    chomp $line;
    my @f = split /\t/, $line;
    next unless @f >= 9;
    my ($fid) = $f[8] =~ /ID=([^;]+)/;
    my ($parent) = $f[8] =~ /Parent=([^;]+)/;
    if ($f[2] eq "gene") {
        $source_span{$fid} = [$f[3], $f[4]] if defined $fid;
    } elsif ($f[2] eq "mRNA" && defined $fid && defined $parent) {
        $tx_gene{$fid} = $parent;
    } elsif ($f[2] eq "CDS" && defined $parent) {
        push @{ $tx_cds{$parent} }, [$f[3], $f[4], ($f[7] =~ /^[012]$/ ? $f[7] : 0)];
    }
}
close $rm;
# per source gene, the CDS of its longest-CDS transcript: used to complete a rebuilt
# model's 5' end with the fragment's own start codon
my %source_cds;
for my $tx (sort keys %tx_cds) {
    my $gene = $tx_gene{$tx} // next;
    my $len = 0; $len += $_->[1] - $_->[0] + 1 for @{ $tx_cds{$tx} };
    my $cur = $source_cds{$gene};
    $source_cds{$gene} = { len => $len, cds => $tx_cds{$tx} } if !$cur || $len > $cur->{len};
}

# genome bases chr:start-end (1-based), on the given strand (reverse complement on -)
sub genome_bases {
    my ($chr, $start, $end, $strand) = @_;
    my @lines = `$samtools_bin faidx $genome_fa '$chr:$start-$end'`;
    shift @lines;
    chomp @lines;
    my $seq = uc join("", @lines);
    if ($strand eq "-") { $seq = reverse $seq; $seq =~ tr/ACGT/TGCA/ }
    return $seq;
}

# GFF3 phases in transcription order from the first CDS's phase (list of [start, end, phase])
sub set_phases {
    my ($tx, $first_phase) = @_;
    my $coding = -$first_phase;
    for my $c (@$tx) {
        $c->[2] = (3 - $coding % 3) % 3;
        $coding += $c->[1] - $c->[0] + 1;
    }
}

open my $rep, ">", $report or die "cant write $report: $!\n";
print $rep join("\t", qw(new_gene_id merge_id ref region prot_len prot_cov identity positive frameshifts stops
                         n_cds sources_covered result reason)), "\n";

my %new_block;   # gid -> [lines] for SUPPORTED genes (replacing the placeholder transcripts)
my %gene_attr;   # gid -> attribute string to append
my ($n_sup, $n_not) = (0, 0);
mkdir $scratch unless -d $scratch;

for my $gid (map { $_->[1] } grep { $_->[0] eq "mender" } @order) {
    my $g = $gene_info{$gid};
    my $ref = $best_hit{ $g->{merge_id} // "" };
    my @report_row = ($gid, $g->{merge_id} // ".", $ref // ".");
    my ($result, $reason);
    unless (defined $ref && defined $ref_seq{$ref} && $ref_seq{$ref} ne "") {
        ($result, $reason) = ("NOT_SUPPORTED", "no_reference_protein");
        print $rep join("\t", @report_row, (".") x 9, $result, $reason), "\n";
        $gene_attr{$gid} = ";miniprot_support=$result:$reason";
        $n_not++;
        next;
    }
    my $r_start = $g->{start} - $flank;  $r_start = 1 if $r_start < 1;
    my $r_end   = $g->{end} + $flank;    $r_end = $seq_len{$g->{chr}} if $r_end > ($seq_len{$g->{chr}} // $r_end);
    my $region  = "$g->{chr}:$r_start-$r_end";
    my $region_fa = "$scratch/mp_region.$$.fa";
    my $prot_fa   = "$scratch/mp_prot.$$.fa";
    system("$samtools_bin faidx $genome_fa '$region' > $region_fa") == 0 or die "samtools faidx $region failed\n";
    open my $pf, ">", $prot_fa or die; print $pf ">$ref\n$ref_seq{$ref}\n"; close $pf;
    my @out = `$miniprot_bin --gff -t$threads $region_fa $prot_fa 2>/dev/null`;
    unlink $region_fa, $prot_fa;

    # the Rank=1 alignment: its PAF line (frameshifts, stops), mRNA and CDS lines
    my ($paf, $mrna, @cds);
    my ($mp_id, $mp_stop);
    for my $line (@out) {
        chomp $line;
        if ($line =~ /^##PAF\t/ && !defined $paf) { $paf = $line; next }
        my @f = split /\t/, $line;
        next unless @f >= 9;
        if ($f[2] eq "mRNA" && !defined $mrna && $f[8] =~ /Rank=1\b/) { $mrna = \@f; ($mp_id) = $f[8] =~ /ID=([^;]+)/; next }
        push @cds, \@f if $f[2] eq "CDS" && defined $mp_id && $f[8] =~ /Parent=\Q$mp_id\E;/;
        $mp_stop = 1 if $f[2] eq "stop_codon" && defined $mp_id && $f[8] =~ /Parent=\Q$mp_id\E;/;
    }
    my $prot_len = length $ref_seq{$ref};
    unless ($mrna && @cds) {
        ($result, $reason) = ("NOT_SUPPORTED", "no_alignment");
        print $rep join("\t", @report_row, $region, $prot_len, (".") x 7, $result, $reason), "\n";
        $gene_attr{$gid} = ";miniprot_support=$result:$reason;miniprot_ref=$ref";
        $n_not++;
        next;
    }
    my ($t_start, $t_end) = $mrna->[8] =~ /Target=\S+ (\d+) (\d+)/;
    my ($identity) = $mrna->[8] =~ /Identity=([\d.]+)/;
    my ($positive) = $mrna->[8] =~ /Positive=([\d.]+)/;
    my ($fs) = ($paf // "") =~ /\tfs:i:(\d+)/;
    my ($st) = ($paf // "") =~ /\tst:i:(\d+)/;
    $fs //= 0; $st //= 0;
    my $prot_cov = ($t_end - $t_start + 1) / $prot_len;

    # to genome coordinates
    my $shift = $r_start - 1;
    my @gcds = map { [ $_->[3] + $shift, $_->[4] + $shift, $_->[6], $_->[7] ] } @cds;

    # every source gene must be overlapped by the model's CDS
    my @sources = @{ $g->{sources} };
    my $covered = 0;
    for my $sg (@sources) {
        my $span = $source_span{$sg} or next;
        $covered++ if grep { $_->[0] <= $span->[1] && $_->[1] >= $span->[0] } @gcds;
    }
    my $n_known = grep { $source_span{$_} } @sources;

    my @why;
    push @why, sprintf("protein_cov_%.2f", $prot_cov) if $prot_cov < $min_prot_cov;
    push @why, "frameshift" if $fs > 0;
    push @why, "inframe_stop" if $st > 0;
    push @why, "covers_${covered}_of_" . scalar(@sources) . "_sources" if $n_known == @sources && $covered < @sources;
    push @why, "source_spans_unknown" if $n_known < @sources;
    ($result, $reason) = @why ? ("NOT_SUPPORTED", join(",", @why)) : ("SUPPORTED", ".");

    print $rep join("\t", @report_row, $region, $prot_len, sprintf("%.3f", $prot_cov), $identity // ".", $positive // ".",
                    $fs, $st, scalar(@gcds), "$covered/" . scalar(@sources), $result, $reason), "\n";
    my $attrs = sprintf(";miniprot_support=%s%s;miniprot_ref=%s;miniprot_prot_cov=%.3f;miniprot_identity=%s;miniprot_positive=%s",
                        $result, ($reason ne "." ? ":$reason" : ""), $ref, $prot_cov, $identity // ".", $positive // ".");
    $gene_attr{$gid} = $attrs;
    if ($result eq "NOT_SUPPORTED") { $n_not++; next }
    $n_sup++;

    # the rebuilt transcript, in transcription order: [start, end, phase]
    my $strand = $g->{strand};
    my @tx = map { [$_->[0], $_->[1], ($_->[3] =~ /^[012]$/ ? $_->[3] : 0)] }
             ($strand eq "-" ? sort { $b->[1] <=> $a->[1] } @gcds : sort { $a->[0] <=> $b->[0] } @gcds);

    # Start codon. miniprot's model begins where the reference alignment begins, which is
    # not a start codon when the N-terminus is poorly conserved. The 5'-most source fragment
    # was called with its own start codon: when miniprot's first coding base lies inside one
    # of that fragment's CDS segments, in the fragment's reading frame, the fragment's CDS
    # upstream of that point (from its ATG) completes the model.
    my $start_note = "incomplete";
    my $first_codon = $tx[0][2] == 0 && $tx[0][1] - $tx[0][0] >= 2
        ? genome_bases($g->{chr}, ($strand eq "-" ? ($tx[0][1] - 2, $tx[0][1]) : ($tx[0][0], $tx[0][0] + 2)), $strand) : "";
    if ($first_codon eq "ATG") {
        $start_note = "present";
    } else {
        my ($five) = sort { $strand eq "-" ? $source_span{$b}[1] <=> $source_span{$a}[1] : $source_span{$a}[0] <=> $source_span{$b}[0] }
                     grep { $source_span{$_} && $source_cds{$_} } @sources;
        my @fcds = $five ? map { [@$_] } @{ $source_cds{$five}{cds} } : ();
        @fcds = $strand eq "-" ? sort { $b->[1] <=> $a->[1] } @fcds : sort { $a->[0] <=> $b->[0] } @fcds;
        if (@fcds && $fcds[0][2] == 0) {
            my $f_codon = genome_bases($g->{chr}, ($strand eq "-" ? ($fcds[0][1] - 2, $fcds[0][1]) : ($fcds[0][0], $fcds[0][0] + 2)), $strand);
            # miniprot's first coding base (after its phase), and the fragment segment holding it
            my $mp_first = $strand eq "-" ? $tx[0][1] - $tx[0][2] : $tx[0][0] + $tx[0][2];
            my ($k) = grep { $fcds[$_][0] <= $mp_first && $mp_first <= $fcds[$_][1] } 0 .. $#fcds;
            if ($f_codon eq "ATG" && defined $k) {
                my $offset = 0;   # fragment coding bases before miniprot's first coding base
                $offset += $fcds[$_][1] - $fcds[$_][0] + 1 for 0 .. $k - 1;
                $offset += $strand eq "-" ? $fcds[$k][1] - $mp_first : $mp_first - $fcds[$k][0];
                if ($offset > 0 && $offset % 3 == 0) {
                    my @head = @fcds[0 .. $k - 1];
                    my $joined = $strand eq "-" ? [$tx[0][0], $fcds[$k][1], 0] : [$fcds[$k][0], $tx[0][1], 0];
                    @tx = (@head, $joined, @tx[1 .. $#tx]);
                    $start_note = "from_fragment:$five";
                }
            }
        }
    }
    set_phases(\@tx, $start_note eq "incomplete" ? $tx[0][2] : 0);
    my $end_note = $mp_stop ? "present" : "absent";
    $attrs .= ";start_codon=$start_note;stop_codon=$end_note";
    $gene_attr{$gid} = $attrs;
    my @sorted = sort { $a->[0] <=> $b->[0] } @tx;
    my ($m_start, $m_end) = ($sorted[0][0], $sorted[-1][1]);
    my ($old_tid) = map { /\tmRNA\t.*ID=([^;]+)/ ? $1 : () } @{ $block{$gid} // [] };
    my $tid = $old_tid // "$gid.t1";
    my @lines = (join("\t", $g->{chr}, "Mender", "mRNA", $m_start, $m_end, ".", $g->{strand}, ".",
                      "ID=$tid;Parent=$gid;Name=$tid;merge_source=Mender;merge_id=$g->{merge_id};model_source=miniprot$attrs"));
    my $n = 0;
    for my $c (@tx) {
        $n++;
        push @lines, join("\t", $g->{chr}, "Mender", "exon", $c->[0], $c->[1], ".", $g->{strand}, ".", "ID=$tid.exon.$n;Parent=$tid");
    }
    $n = 0;
    for my $c (@tx) {
        $n++;
        push @lines, join("\t", $g->{chr}, "Mender", "CDS", $c->[0], $c->[1], ".", $g->{strand}, $c->[2], "ID=$tid.CDS.$n;Parent=$tid");
    }
    $new_block{$gid} = { lines => \@lines, start => $m_start, end => $m_end };
}
close $rep;

# write: SUPPORTED genes get the rebuilt model (gene span = model span); others keep theirs
open my $out, ">", $out_gff or die "cant write $out_gff: $!\n";
for my $item (@order) {
    if ($item->[0] eq "line") { print $out "$item->[1]\n"; next }
    my $gid = $item->[1];
    my @f = split /\t/, $gene_info{$gid}{line};
    if (my $nb = $new_block{$gid}) {
        @f[3, 4] = ($nb->{start}, $nb->{end});
        $f[8] .= ($gene_attr{$gid} // "") . ";model_source=miniprot";
        print $out join("\t", @f), "\n";
        print $out "$_\n" for @{ $nb->{lines} };
    } else {
        $f[8] .= $gene_attr{$gid} // "";
        print $out join("\t", @f), "\n";
        print $out "$_\n" for @{ $block{$gid} // [] };
    }
}
close $out;
printf STDERR "miniprot rebuild: %d SUPPORTED (rebuilt), %d NOT_SUPPORTED (placeholder kept)\n", $n_sup, $n_not;
