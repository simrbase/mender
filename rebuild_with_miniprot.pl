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
## the break is still missing. Here a reference protein the fragments tiled
## against is aligned to the merged locus with miniprot (spliced
## protein-to-genome alignment). A true split gives ONE
## alignment covering most of the reference protein and every source fragment,
## with real splice sites, no frameshift and no in-frame stop: that alignment
## becomes the merged gene's model. Two different genes joined by mistake give an
## alignment that covers only one of them.
##
## Per merged gene (Mender gene lines in the merged GFF):
##   region  = merged gene span +/- --flank bp (clipped to the sequence)
##   proteins = rebuild_candidates of its merge_id in the merge table (ranked by
##              find_split_genes.pl: most chain genes, then identity-weighted
##              coverage), or tiling_ref when the table has no rebuild_candidates
##              column
##   miniprot --gff on (region, protein), one candidate at a time, best-ranked first,
##   at most --max_refs; the Rank=1 alignment is judged, and the first SUPPORTED
##   candidate is used. A distant reference aligned to the genome makes spurious
##   frameshifts; a closer one of the same family often aligns cleanly:
##     SUPPORTED     protein coverage >= --min_prot_cov, no frameshift (fs:i:0),
##                   no in-frame stop (st:i:0), and its CDS overlap every source gene
##     NOT_SUPPORTED otherwise (the reason is recorded)
##   SUPPORTED:     the gene's transcripts are replaced by the miniprot model: one
##                  mRNA, phases from miniprot; gene span = model. miniprot gives the CDS
##                  only; UTRs come back from the source genes (see UTRs below)
##   NOT_SUPPORTED: no candidate supports it; the placeholder model is kept and
##                  validate_merge_translation.pl fails it (NO_PROTEIN_SPAN). The
##                  best-ranked candidate's alignment is what is reported.
## Every candidate must itself tile the fragments (find_split_genes.pl lists only
## references hitting two or more of them), and SUPPORTED always needs one
## alignment covering every source gene, so trying several does not let an
## unrelated protein through.
## UTRs: a protein has none, so a rebuilt model is CDS-only unless its ends are the source
## genes' ends. When the rebuilt start codon is at the same genome position as the 5' source
## gene's CDS start, that gene's 5' UTR (its exons upstream of the start, and the untranslated
## part of the first coding exon) is put back; the same for the 3' UTR when the rebuilt CDS
## ends where the 3' source gene's CDS ends (miniprot's CDS includes the stop codon, as the
## source's does). A UTR is never moved to a start or stop that differs: it would sit inside
## the new CDS or away from it. Recorded as utr5= / utr3= from_fragment:<transcript> | none,
## with five_prime_UTR / three_prime_UTR features.
## Isoforms: a rebuilt gene has ONE transcript. The reference protein defines one coding
## model; alternative splicing is evidence it cannot give. Source transcripts contribute only
## the ends (start / stop codon completion, UTRs), each from the first transcript of the end
## fragment, longest CDS first, that fits there; the notes name it (from_fragment:<transcript>).
## When a source gene had several transcripts (the placeholder carried their N x M
## cross-product), the gene records placeholder_transcripts=, source_isoforms=<gene>:<n>, and
## isoforms_not_used= (the multi-isoform source transcripts nothing was taken from), so no
## isoform is dropped silently.
## Attributes on the gene (and mRNA): miniprot_support=SUPPORTED|NOT_SUPPORTED:<reason>,
## miniprot_ref, miniprot_ref_rank (1 = best-ranked candidate), miniprot_refs_tried,
## miniprot_prot_cov, miniprot_identity, miniprot_positive.
## Report (--report, one row per merged gene): rebuild_ref is the reference the model
## was rebuilt from (SUPPORTED), or the best-ranked candidate whose alignment is
## reported (NOT_SUPPORTED); rebuild_ref_desc is its description from --ref_fa.
##
## Usage:
##   perl rebuild_with_miniprot.pl --merged_gff merges.gff --removed_gff removed.gff \
##        --merge_table merge_candidates.txt \
##        --ref_fa ref.fa --genome_fa genome.fa --out merges.miniprot.gff --report miniprot.tsv
##        [--flank 5000] [--min_prot_cov 0.80] [--miniprot_bin miniprot] [--samtools_bin samtools]
##        [--threads 1] [--scratch DIR] [--max_refs 10]

my ($merged_gff, $merge_table, $ref_fa, $genome_fa, $out_gff, $report, $removed_gff);
my $flank        = 5000;
my $min_prot_cov = 0.80;
my $miniprot_bin = "miniprot";
my $samtools_bin = "samtools";
my $threads      = 1;
my $scratch      = ".";
my $max_refs     = 10;   # = the candidates find_split_genes.pl writes; often several isoforms of one gene
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
    "max_refs=i"     => \$max_refs,
) or die "Error in options\n";
for my $required ([merged_gff => $merged_gff], [removed_gff => $removed_gff], [merge_table => $merge_table], [ref_fa => $ref_fa],
                  [genome_fa => $genome_fa], [out => $out_gff], [report => $report]) {
    die "--$required->[0] is required\n" unless defined $required->[1];
}
die "--max_refs must be at least 1\n" if $max_refs < 1;
die "genome index $genome_fa.fai missing (samtools faidx $genome_fa)\n" unless -e "$genome_fa.fai";

# sequence lengths, for clipping the flank
my %seq_len;
open my $fai, "<", "$genome_fa.fai" or die "cant open $genome_fa.fai: $!\n";
while (my $line = <$fai>) { my ($name, $len) = split /\t/, $line; $seq_len{$name} = $len }
close $fai;

# merge_id -> reference proteins to try, best-ranked first: rebuild_candidates, then
# tiling_ref (already in the list when find_split_genes.pl wrote it)
my %candidates;
open my $mt, "<", $merge_table or die "cant open $merge_table: $!\n";
my $hdr_line = <$mt>;
chomp $hdr_line;
my @hdr = split /\t/, $hdr_line;
my %col;
for my $i (0 .. $#hdr) { $col{ $hdr[$i] } = $i }
die "merge table needs merge_id and rebuild_candidates or tiling_ref columns\n"
    unless defined $col{merge_id} && (defined $col{rebuild_candidates} || defined $col{tiling_ref});
while (my $line = <$mt>) {
    chomp $line;
    my @f = split /\t/, $line;
    my $listed = defined $col{rebuild_candidates} ? ($f[ $col{rebuild_candidates} ] // "") : "";
    my $tiling = defined $col{tiling_ref} ? ($f[ $col{tiling_ref} ] // "") : "";
    my @refs;
    my %seen;
    for my $ref (split(/,/, $listed), $tiling) {
        next if $ref eq "" || $ref eq "?" || $seen{$ref}++;
        push @refs, $ref;
    }
    splice @refs, $max_refs if @refs > $max_refs;
    $candidates{ $f[ $col{merge_id} ] } = \@refs;
}
close $mt;

# reference proteins needed
my %wanted;
for my $refs (values %candidates) {
    for my $ref (@$refs) { $wanted{$ref} = 1 }
}
# descriptions as find_split_genes.pl reads them: gene_symbol: and description: fields
my (%ref_seq, %ref_desc, $id);
open my $rf, "<", $ref_fa or die "cant open $ref_fa: $!\n";
while (my $line = <$rf>) {
    chomp $line;
    if ($line =~ /^>(\S+)/) {
        $id = $wanted{$1} ? $1 : undef;
        next unless defined $id;
        $ref_seq{$id} = "";
        my ($symbol)  = $line =~ /gene_symbol:(\S+)/;
        my ($text)    = $line =~ /description:(.+?)(?:\s+\[|$)/;
        my @parts;
        for my $part ($symbol, $text) { push @parts, $part if defined $part }
        $ref_desc{$id} = @parts ? join(" | ", @parts) : ".";
        next;
    }
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
my (%tx_gene, %tx_cds, %tx_exons);   # source transcript -> gene; -> [[start, end, phase]]; -> [[start, end]]
my %transcripts_of;                  # source gene -> [its transcripts]
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
        push @{ $transcripts_of{$parent} }, $fid;
    } elsif ($f[2] eq "CDS" && defined $parent) {
        push @{ $tx_cds{$parent} }, [$f[3], $f[4], ($f[7] =~ /^[012]$/ ? $f[7] : 0)];
    } elsif ($f[2] eq "exon" && defined $parent) {
        push @{ $tx_exons{$parent} }, [$f[3], $f[4]];
    }
}
close $rm;
# per source gene, its transcripts (CDS and exons), longest CDS first: the ends of a
# rebuilt model are completed (start / stop codon) and given UTRs from the first transcript
# of the end fragment that fits there -- not blindly the longest, which in a multi-isoform
# gene can be an alternative ORF that has nothing to do with the rebuilt model
my %source_txs;
sub by_cds_length_desc($$) { my ($left, $right) = @_; return $right->{len} <=> $left->{len} || $left->{tx} cmp $right->{tx} }
for my $tx (sort keys %tx_cds) {
    my $gene = $tx_gene{$tx} // next;
    my $len = 0;
    for my $seg (@{ $tx_cds{$tx} }) { $len += $seg->[1] - $seg->[0] + 1 }
    push @{ $source_txs{$gene} }, { len => $len, cds => $tx_cds{$tx}, exons => $tx_exons{$tx} // [], tx => $tx };
}
for my $gene (keys %source_txs) { @{ $source_txs{$gene} } = sort by_cds_length_desc @{ $source_txs{$gene} } }

# segment comparators for sort ([start, end, ...]); sort passes the two elements in @_
sub by_start_ascending($$) { my ($left, $right) = @_; return $left->[0] <=> $right->[0] }
sub by_end_descending($$)  { my ($left, $right) = @_; return $right->[1] <=> $left->[1] }

# copies of the segments, in transcription order (5' to 3' on the strand)
sub transcription_order {
    my ($strand, @segments) = @_;
    my @copies;
    for my $seg (@segments) { push @copies, [@$seg] }
    my @ordered = $strand eq "-" ? sort(by_end_descending @copies) : sort(by_start_ascending @copies);
    return @ordered;
}

# index of the segment holding genome position $pos, or undef
sub segment_holding {
    my ($pos, $segments) = @_;
    for my $i (0 .. $#$segments) {
        return $i if $segments->[$i][0] <= $pos && $pos <= $segments->[$i][1];
    }
    return undef;
}

# the source gene at the 5' (or 3') end of the merged gene, among those with a known span and CDS
sub end_source {
    my ($which, $strand, @sources) = @_;
    my $want_low = ($which eq "5prime") == ($strand ne "-");   # lowest start: 5' on +, 3' on -
    my ($best, $best_pos);
    for my $sg (@sources) {
        next unless $source_span{$sg} && $source_txs{$sg};
        my $pos = $want_low ? $source_span{$sg}[0] : $source_span{$sg}[1];
        if (!defined $best || ($want_low ? $pos < $best_pos : $pos > $best_pos)) { ($best, $best_pos) = ($sg, $pos) }
    }
    return $best;
}

# A source transcript's UTR at one end of its CDS: $boundary is its first coding base
# ("5prime") or last ("3prime"). Returns the genome coordinate the terminal coding exon
# extends to (the UTR part of the exon holding the boundary; undef if it has none) and the
# whole UTR exons beyond it, in transcription order.
sub source_utr {
    my ($which, $strand, $boundary, $exons) = @_;
    my $below = ($which eq "5prime") == ($strand ne "-");   # the UTR lies at lower coordinates
    my ($extend, @whole);
    for my $exon (@$exons) {
        if ($below) {
            if    ($exon->[1] < $boundary)                            { push @whole, [@$exon] }
            elsif ($exon->[0] < $boundary && $boundary <= $exon->[1]) { $extend = $exon->[0] }
        } else {
            if    ($exon->[0] > $boundary)                            { push @whole, [@$exon] }
            elsif ($exon->[0] <= $boundary && $boundary < $exon->[1]) { $extend = $exon->[1] }
        }
    }
    return ($extend, transcription_order($strand, @whole));
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

# the first codon of a segment in transcription order, and the last
sub first_codon_of { my ($chr, $seg, $strand) = @_; return genome_bases($chr, ($strand eq "-" ? ($seg->[1] - 2, $seg->[1]) : ($seg->[0], $seg->[0] + 2)), $strand) }
sub last_codon_of  { my ($chr, $seg, $strand) = @_; return genome_bases($chr, ($strand eq "-" ? ($seg->[0], $seg->[0] + 2) : ($seg->[1] - 2, $seg->[1])), $strand) }

# GFF3 phases in transcription order from the first CDS's phase (list of [start, end, phase])
sub set_phases {
    my ($tx, $first_phase) = @_;
    my $coding = -$first_phase;
    for my $c (@$tx) {
        $c->[2] = (3 - $coding % 3) % 3;
        $coding += $c->[1] - $c->[0] + 1;
    }
}

# Align one reference protein to the region and judge the Rank=1 alignment.
# Returns { ref, result, reason, prot_len, prot_cov, identity, positive, fs, st, gcds, covered, stop }.
sub try_reference {
    my ($g, $ref, $region_fa, $r_start) = @_;
    my %try = (ref => $ref, prot_len => length $ref_seq{$ref});
    my $prot_fa = "$scratch/mp_prot.$$.fa";
    open my $pf, ">", $prot_fa or die "cant write $prot_fa: $!\n";
    print $pf ">$ref\n$ref_seq{$ref}\n";
    close $pf;
    my @out = `$miniprot_bin --gff -t$threads $region_fa $prot_fa 2>/dev/null`;
    unlink $prot_fa;

    # the Rank=1 alignment: its PAF line (frameshifts, stops), mRNA and CDS lines
    my ($paf, $mrna, @cds, $mp_id);
    for my $line (@out) {
        chomp $line;
        if ($line =~ /^##PAF\t/ && !defined $paf) { $paf = $line; next }
        my @f = split /\t/, $line;
        next unless @f >= 9;
        if ($f[2] eq "mRNA" && !defined $mrna && $f[8] =~ /Rank=1\b/) { $mrna = \@f; ($mp_id) = $f[8] =~ /ID=([^;]+)/; next }
        next unless defined $mp_id && $f[8] =~ /Parent=\Q$mp_id\E(?:;|$)/;
        push @cds, \@f if $f[2] eq "CDS";
        $try{stop} = 1 if $f[2] eq "stop_codon";
    }
    unless ($mrna && @cds) {
        @try{qw(result reason)} = ("NOT_SUPPORTED", "no_alignment");
        return \%try;
    }
    my ($t_start, $t_end) = $mrna->[8] =~ /Target=\S+ (\d+) (\d+)/;
    ($try{identity}) = $mrna->[8] =~ /Identity=([\d.]+)/;
    ($try{positive}) = $mrna->[8] =~ /Positive=([\d.]+)/;
    ($try{fs}) = ($paf // "") =~ /\tfs:i:(\d+)/;
    ($try{st}) = ($paf // "") =~ /\tst:i:(\d+)/;
    $try{fs} //= 0;
    $try{st} //= 0;
    $try{prot_cov} = ($t_end - $t_start + 1) / $try{prot_len};

    # to genome coordinates: [start, end, strand, phase]
    my $shift = $r_start - 1;
    my @gcds;
    for my $c (@cds) { push @gcds, [ $c->[3] + $shift, $c->[4] + $shift, $c->[6], $c->[7] ] }
    $try{gcds} = \@gcds;

    # every source gene must be overlapped by the model's CDS
    my @sources = @{ $g->{sources} };
    my ($covered, $n_known) = (0, 0);
    for my $sg (@sources) {
        my $span = $source_span{$sg} or next;
        $n_known++;
        for my $c (@gcds) {
            if ($c->[0] <= $span->[1] && $c->[1] >= $span->[0]) { $covered++; last }
        }
    }
    $try{covered} = $covered;

    my @why;
    push @why, sprintf("protein_cov_%.2f", $try{prot_cov}) if $try{prot_cov} < $min_prot_cov;
    push @why, "frameshift" if $try{fs} > 0;
    push @why, "inframe_stop" if $try{st} > 0;
    push @why, "covers_${covered}_of_" . scalar(@sources) . "_sources" if $n_known == @sources && $covered < @sources;
    push @why, "source_spans_unknown" if $n_known < @sources;
    @try{qw(result reason)} = @why ? ("NOT_SUPPORTED", join(",", @why)) : ("SUPPORTED", ".");
    return \%try;
}

open my $rep, ">", $report or die "cant write $report: $!\n";
print $rep join("\t", qw(new_gene_id merge_id rebuild_ref region prot_len prot_cov identity positive frameshifts stops
                         n_cds sources_covered result reason ref_rank refs_tried rebuild_ref_desc)), "\n";

my %new_block;   # gid -> [lines] for SUPPORTED genes (replacing the placeholder transcripts)
my %gene_attr;   # gid -> attribute string to append
my ($n_sup, $n_not, $n_sup_later, $n_utr5, $n_utr3, $n_iso_genes, $n_iso_dropped) = (0, 0, 0, 0, 0, 0, 0);
mkdir $scratch unless -d $scratch;

for my $item (@order) {
    next unless $item->[0] eq "mender";
    my $gid = $item->[1];
    my $g = $gene_info{$gid};
    my @refs;
    for my $ref (@{ $candidates{ $g->{merge_id} // "" } // [] }) {
        push @refs, $ref if defined $ref_seq{$ref} && $ref_seq{$ref} ne "";
    }
    unless (@refs) {
        my $listed = $candidates{ $g->{merge_id} // "" } // [];
        print $rep join("\t", $gid, $g->{merge_id} // ".", $listed->[0] // ".", (".") x 9,
                        "NOT_SUPPORTED", "no_reference_protein", ".", 0, "."), "\n";
        $gene_attr{$gid} = ";miniprot_support=NOT_SUPPORTED:no_reference_protein";
        $n_not++;
        next;
    }
    my $r_start = $g->{start} - $flank;  $r_start = 1 if $r_start < 1;
    my $r_end   = $g->{end} + $flank;    $r_end = $seq_len{$g->{chr}} if $r_end > ($seq_len{$g->{chr}} // $r_end);
    my $region  = "$g->{chr}:$r_start-$r_end";
    my $region_fa = "$scratch/mp_region.$$.fa";
    system("$samtools_bin faidx $genome_fa '$region' > $region_fa") == 0 or die "samtools faidx $region failed\n";

    # candidates in rank order; the first SUPPORTED one is used, else the best-ranked is reported
    my ($chosen, $first, $rank, $tried) = (undef, undef, 0, 0);
    for my $i (0 .. $#refs) {
        my $try = try_reference($g, $refs[$i], $region_fa, $r_start);
        $tried++;
        $first //= $try;
        if ($try->{result} eq "SUPPORTED") { ($chosen, $rank) = ($try, $i + 1); last }
    }
    unlink $region_fa;
    my $use = $chosen // $first;
    $rank = 1 unless $chosen;

    my ($result, $reason) = @{$use}{qw(result reason)};
    my @row = ($gid, $g->{merge_id} // ".", $use->{ref}, $region, $use->{prot_len});
    if ($use->{gcds}) {
        push @row, sprintf("%.3f", $use->{prot_cov}), $use->{identity} // ".", $use->{positive} // ".",
                   $use->{fs}, $use->{st}, scalar(@{ $use->{gcds} }), "$use->{covered}/" . scalar(@{ $g->{sources} });
    } else {
        push @row, (".") x 7;
    }
    print $rep join("\t", @row, $result, $reason, $rank, $tried, $ref_desc{ $use->{ref} } // "."), "\n";

    my $attrs = ";miniprot_support=$result" . ($reason ne "." ? ":$reason" : "")
              . ";miniprot_ref=$use->{ref};miniprot_ref_rank=$rank;miniprot_refs_tried=$tried";
    $attrs .= sprintf(";miniprot_prot_cov=%.3f;miniprot_identity=%s;miniprot_positive=%s",
                      $use->{prot_cov}, $use->{identity} // ".", $use->{positive} // ".") if $use->{gcds};
    $gene_attr{$gid} = $attrs;
    unless ($chosen) { $n_not++; next }
    $n_sup++;
    $n_sup_later++ if $rank > 1;

    # the rebuilt transcript, in transcription order: [start, end, phase]
    my $strand = $g->{strand};
    my @tx;
    for my $c (transcription_order($strand, @{ $use->{gcds} })) {
        push @tx, [$c->[0], $c->[1], ($c->[3] =~ /^[012]$/ ? $c->[3] : 0)];
    }
    my @sources = @{ $g->{sources} };

    # Start codon. miniprot's model begins where the reference alignment begins, which is
    # not a start codon when the N-terminus is poorly conserved. The 5'-most source fragment
    # was called with its own start codon: when miniprot's first coding base lies inside one
    # of that fragment's CDS segments, in the fragment's reading frame, the fragment's CDS
    # upstream of that point (from its ATG) completes the model.
    my $start_note = "incomplete";
    my $first_codon = $tx[0][2] == 0 && $tx[0][1] - $tx[0][0] >= 2 ? first_codon_of($g->{chr}, $tx[0], $strand) : "";
    my $five = end_source("5prime", $strand, @sources);
    if ($first_codon eq "ATG") {
        $start_note = "present";
    } elsif ($five) {
        # the fragment's transcripts, longest CDS first; the first that fits completes it
        for my $cand (@{ $source_txs{$five} }) {
            my @fcds = transcription_order($strand, @{ $cand->{cds} });
            next unless $fcds[0][2] == 0 && first_codon_of($g->{chr}, $fcds[0], $strand) eq "ATG";
            # miniprot's first coding base (after its phase), and the fragment segment holding it
            my $mp_first = $strand eq "-" ? $tx[0][1] - $tx[0][2] : $tx[0][0] + $tx[0][2];
            my $k = segment_holding($mp_first, \@fcds);
            next unless defined $k;
            my $offset = 0;   # fragment coding bases before miniprot's first coding base
            for my $i (0 .. $k - 1) { $offset += $fcds[$i][1] - $fcds[$i][0] + 1 }
            $offset += $strand eq "-" ? $fcds[$k][1] - $mp_first : $mp_first - $fcds[$k][0];
            next unless $offset > 0 && $offset % 3 == 0;
            my @head = @fcds[0 .. $k - 1];
            my $joined = $strand eq "-" ? [$tx[0][0], $fcds[$k][1], 0] : [$fcds[$k][0], $tx[0][1], 0];
            @tx = (@head, $joined, @tx[1 .. $#tx]);
            $start_note = "from_fragment:$cand->{tx}";
            last;
        }
    }
    set_phases(\@tx, $start_note eq "incomplete" ? $tx[0][2] : 0);

    # Stop codon, the same way at the 3' end: when miniprot's model ends without a stop codon
    # (the reference's C-terminus is not conserved), the 3'-most source fragment was called
    # with its own stop. When miniprot's last coding base lies inside one of that fragment's
    # CDS segments, at the end of a codon in the fragment's reading frame, and the fragment
    # ends in a stop codon, the fragment's CDS from there on completes the model. Otherwise
    # the model is kept as it is, marked stop_codon=absent -- never extended into genome
    # sequence no annotation covered, where the first stop found is usually the wrong one.
    my $end_note = $use->{stop} ? "present" : "absent";
    my $three = end_source("3prime", $strand, @sources);
    if (!$use->{stop} && $three) {
        my $model_coding = -$tx[0][2];
        for my $seg (@tx) { $model_coding += $seg->[1] - $seg->[0] + 1 }
        # the fragment's transcripts, longest CDS first; the first that fits completes it
        for my $cand ($model_coding % 3 == 0 ? @{ $source_txs{$three} } : ()) {
            my @fcds = transcription_order($strand, @{ $cand->{cds} });
            next unless last_codon_of($g->{chr}, $fcds[-1], $strand) =~ /^(?:TAA|TAG|TGA)$/;
            my $mp_last = $strand eq "-" ? $tx[-1][0] : $tx[-1][1];
            my $k = segment_holding($mp_last, \@fcds);
            next unless defined $k;
            my $offset = -$fcds[0][2];   # fragment coding bases up to and including miniprot's last base
            for my $i (0 .. $k - 1) { $offset += $fcds[$i][1] - $fcds[$i][0] + 1 }
            $offset += $strand eq "-" ? $fcds[$k][1] - $mp_last + 1 : $mp_last - $fcds[$k][0] + 1;
            my $at_segment_end = $strand eq "-" ? $mp_last == $fcds[$k][0] : $mp_last == $fcds[$k][1];
            next unless $offset > 0 && $offset % 3 == 0 && !($at_segment_end && $k == $#fcds);
            if ($strand eq "-") { $tx[-1][0] = $fcds[$k][0] } else { $tx[-1][1] = $fcds[$k][1] }
            push @tx, @fcds[$k + 1 .. $#fcds];
            $end_note = "from_fragment:$cand->{tx}";
            set_phases(\@tx, $tx[0][2]);
            last;
        }
    }
    # UTRs from the source genes, only where the rebuilt CDS starts (ends) exactly where the
    # 5' (3') source gene's CDS does
    my @exons;
    for my $c (@tx) { push @exons, [$c->[0], $c->[1]] }
    my %utr_note = (5 => "none", 3 => "none");
    for my $which ("5prime", "3prime") {
        my $source = ($which eq "5prime" ? $five : $three) // next;
        my $is5   = $which eq "5prime";
        my $low   = ($is5 == ($strand ne "-"));   # this end is the lower genome coordinate
        my $model_end = $is5 ? ($low ? $tx[0][0] : $tx[0][1]) : ($low ? $tx[-1][0] : $tx[-1][1]);
        # the fragment's transcripts, longest CDS first: the first one whose CDS ends here and
        # that has a UTR there
        for my $cand (@{ $source_txs{$source} }) {
            my @scds = transcription_order($strand, @{ $cand->{cds} });
            my $source_end = $is5 ? ($low ? $scds[0][0] : $scds[0][1]) : ($low ? $scds[-1][0] : $scds[-1][1]);
            next unless $model_end == $source_end;
            my ($extend, @whole) = source_utr($which, $strand, $source_end, $cand->{exons});
            next unless defined $extend || @whole;
            my $terminal = $is5 ? $exons[0] : $exons[-1];
            $terminal->[$low ? 0 : 1] = $extend if defined $extend;
            if ($is5) { unshift @exons, @whole } else { push @exons, @whole }
            $utr_note{ $is5 ? 5 : 3 } = "from_fragment:$cand->{tx}";
            last;
        }
    }
    $n_utr5++ if $utr_note{5} ne "none";
    $n_utr3++ if $utr_note{3} ne "none";
    # UTR features: the parts of exons outside the CDS span
    my @cds_sorted = sort by_start_ascending @tx;
    my ($cds_lo, $cds_hi) = ($cds_sorted[0][0], $cds_sorted[-1][1]);
    my @utrs;   # [type, start, end]
    my ($low_type, $high_type) = $strand eq "-" ? ("three_prime_UTR", "five_prime_UTR") : ("five_prime_UTR", "three_prime_UTR");
    for my $exon (sort by_start_ascending @exons) {
        push @utrs, [$low_type,  $exon->[0], ($exon->[1] < $cds_lo ? $exon->[1] : $cds_lo - 1)] if $exon->[0] < $cds_lo;
        push @utrs, [$high_type, ($exon->[0] > $cds_hi ? $exon->[0] : $cds_hi + 1), $exon->[1]] if $exon->[1] > $cds_hi;
    }

    $attrs .= ";start_codon=$start_note;stop_codon=$end_note;utr5=$utr_note{5};utr3=$utr_note{3}";

    # isoforms: one rebuilt transcript replaces the placeholder's cross-product; record the
    # source transcripts it does not draw on
    my $n_placeholder = 0;
    for my $line (@{ $block{$gid} // [] }) { $n_placeholder++ if $line =~ /\tmRNA\t/ }
    my %used;
    for my $note ($start_note, $end_note, $utr_note{5}, $utr_note{3}) { $used{$1} = 1 if $note =~ /^from_fragment:(.+)$/ }
    my (@multi, @not_used);
    for my $source (@sources) {
        my @source_tx = sort @{ $transcripts_of{$source} // [] };
        next unless @source_tx > 1;
        push @multi, "$source:" . scalar(@source_tx);
        for my $source_tx (@source_tx) { push @not_used, $source_tx unless $used{$source_tx} }
    }
    if (@multi) {
        $attrs .= ";placeholder_transcripts=$n_placeholder;source_isoforms=" . join(",", @multi)
                . ";isoforms_not_used=" . (@not_used ? join(",", @not_used) : "none");
        $n_iso_genes++;
        $n_iso_dropped += @not_used;
    }
    $gene_attr{$gid} = $attrs;
    my @sorted = sort by_start_ascending @exons;
    my ($m_start, $m_end) = ($sorted[0][0], $sorted[-1][1]);
    my $old_tid;
    for my $line (@{ $block{$gid} // [] }) {
        if ($line =~ /\tmRNA\t.*ID=([^;]+)/) { $old_tid = $1; last }
    }
    my $tid = $old_tid // "$gid.t1";
    my @lines = (join("\t", $g->{chr}, "Mender", "mRNA", $m_start, $m_end, ".", $strand, ".",
                      "ID=$tid;Parent=$gid;Name=$tid;merge_source=Mender;merge_id=$g->{merge_id};model_source=miniprot$attrs"));
    my $n = 0;
    for my $exon (@exons) {
        $n++;
        push @lines, join("\t", $g->{chr}, "Mender", "exon", $exon->[0], $exon->[1], ".", $strand, ".", "ID=$tid.exon.$n;Parent=$tid");
    }
    $n = 0;
    for my $c (@tx) {
        $n++;
        push @lines, join("\t", $g->{chr}, "Mender", "CDS", $c->[0], $c->[1], ".", $strand, $c->[2], "ID=$tid.CDS.$n;Parent=$tid");
    }
    my %n_type;
    for my $utr (@utrs) {
        my ($type, $u_start, $u_end) = @$utr;
        $n_type{$type}++;
        push @lines, join("\t", $g->{chr}, "Mender", $type, $u_start, $u_end, ".", $strand, ".", "ID=$tid.$type.$n_type{$type};Parent=$tid");
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
    my $nb = $new_block{$gid};
    if ($nb) {
        @f[3, 4] = ($nb->{start}, $nb->{end});
        $f[8] .= ($gene_attr{$gid} // "") . ";model_source=miniprot";
    } else {
        $f[8] .= $gene_attr{$gid} // "";
    }
    print $out join("\t", @f), "\n";
    for my $line ($nb ? @{ $nb->{lines} } : @{ $block{$gid} // [] }) { print $out "$line\n" }
}
close $out;
printf STDERR "miniprot rebuild: %d SUPPORTED (rebuilt; %d by a lower-ranked reference), %d NOT_SUPPORTED (placeholder kept)\n",
    $n_sup, $n_sup_later, $n_not;
printf STDERR "  UTRs put back from the source genes: 5' on %d, 3' on %d of the rebuilt models\n", $n_utr5, $n_utr3;
printf STDERR "  multi-isoform sources: %d rebuilt genes, %d source transcripts not used (isoforms_not_used=)\n",
    $n_iso_genes, $n_iso_dropped if $n_iso_genes;
