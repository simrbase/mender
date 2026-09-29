#!/usr/bin/perl
use strict;
use warnings;
use File::Temp qw(tempdir);
use FindBin;

# tests/test_miniprot_rebuild.pl -- merge_split_genes.pl + rebuild_with_miniprot.pl on a made-up genome:
#   M6  a gene split in two by the annotation, fragment A also 20 codons short at the break:
#       SUPPORTED; the rebuilt model has exactly the gene's four exons (real splice sites) and
#       translates to the reference protein
#   M7  two different adjacent genes merged by mistake: NOT_SUPPORTED (the reference covers one)
#   M8  the reference's N-terminus differs from the gene's: miniprot starts downstream of the start
#       codon; the 5' fragment's own ATG completes the model (start_codon=from_fragment)
#   M9  the reference's C-terminus differs: miniprot ends before the stop codon; the 3' fragment's own
#       stop completes the model (stop_codon=from_fragment)
#   M10 ref_candidates lists an unrelated protein first, then the gene's: the first is tried and
#       fails, the second SUPPORTS it (miniprot_ref_rank=2, refs_tried=2)
#   M11 the false merge given two candidates is still NOT_SUPPORTED, and the best-ranked
#       candidate's result is the one reported
#   M12 --max_refs 1: only the first candidate is tried, so M10's merge is NOT_SUPPORTED
# Needs miniprot and samtools in PATH (skipped otherwise). Exit 0 = all pass.

for my $tool (qw(miniprot samtools)) {
  if (system("which $tool > /dev/null 2>&1") != 0) { print "SKIP: $tool not in PATH\n"; exit 0 }
}
my $dir = tempdir(CLEANUP => 1);
my ($passed, $failed) = (0, 0);
# a regex match passed straight to check() is wrapped in scalar(): in list context a FAILED match
# returns an empty list, which shifts the arguments and makes the description the (true) result
sub check {
  my ($ok, $what, $detail) = @_;
  if ($ok) { $passed++ } else { $failed++; print "FAIL: $what\n" . (defined $detail ? "      got: $detail\n" : '') }
}
sub write_file { my ($p, $t) = @_; open my $fh, '>', $p or die "$p: $!"; print $fh $t; close $fh }

srand(11);
my @AA = split //, 'ACDEFGHIKLMNPQRSTVWY';
my %BACK = (A=>'GCT',C=>'TGT',D=>'GAT',E=>'GAA',F=>'TTT',G=>'GGT',H=>'CAT',I=>'ATT',K=>'AAA',L=>'CTG',
            M=>'ATG',N=>'AAT',P=>'CCT',Q=>'CAA',R=>'CGT',S=>'TCT',T=>'ACT',V=>'GTT',W=>'TGG',Y=>'TAT');
my %CODON = reverse %BACK;
@CODON{qw(TAA TAG TGA)} = ('*') x 3;
sub random_residues { my ($n) = @_; my $p = ''; for my $i (1 .. $n) { $p .= $AA[int rand @AA] } return $p }
sub random_protein { my ($n) = @_; return 'M' . random_residues($n - 1) }
sub random_dna { my ($n) = @_; my @nt = qw(A C G T); my $s = ''; for my $i (1 .. $n) { $s .= $nt[int rand 4] } return $s }
sub intron { my ($n) = @_; my $mid = random_dna($n - 4); $mid =~ s/AG$/AC/; return "GT$mid" . "AG" }
sub encode { my ($p) = @_; my $s = ''; for my $aa (split //, $p) { $s .= $BACK{$aa} } return $s }
# named comparator: sort passes the two elements in @_ (the ($$) prototype)
sub by_start($$) { my ($left, $right) = @_; return $left->[1] <=> $right->[1] }
sub translate { my ($s) = @_; my $p = ''; for (my $i = 0; $i + 3 <= length $s; $i += 3) { $p .= $CODON{substr $s, $i, 3} // 'X' } return $p }

# ---- a gene of four exons; the annotation split it after exon 2 (fragment A: exons 1-2 + a stop
# read into intron 2; fragment B: from an internal ATG at the start of exon 3)
my $protein = random_protein(240);
substr($protein, 120, 1) = 'M';                       # exon 3 begins at codon 121: an ATG B can start at
my $gene_cds = encode($protein) . 'TAA';
my @cut = (0, 150, 360, 480, length $gene_cds);           # exon boundaries in the CDS (codon-aligned at 360)
my @exon_seq;
for my $i (0 .. 3) { push @exon_seq, substr($gene_cds, $cut[$i], $cut[$i + 1] - $cut[$i]) }
my @intron_len = (220, 400, 260);

# chrS: the split gene;  chrW: two unrelated genes side by side (M7)
my ($seq, @exons) = (random_dna(3000));
for my $i (0 .. 3) {
  push @exons, [length($seq) + 1, length($seq) + length $exon_seq[$i]];
  $seq .= $exon_seq[$i];
  $seq .= intron($intron_len[$i]) if $i < 3;
}
$seq .= random_dna(3000);
my %chr = (chrS => $seq);
# the genome sequence of CDS segments [chr, start, end, phase], joined in order
sub cds_sequence { my @segs = @_; my $nt = ''; for my $seg (@segs) { $nt .= substr($chr{ $seg->[0] }, $seg->[1] - 1, $seg->[2] - $seg->[1] + 1) } return $nt }

# M8: the same genome, but the reference protein's first 30 aa differ from the gene's
my $ref_ntdiff = random_protein(30) . substr($protein, 30);
my $ref_ctdiff = substr($protein, 0, -30) . random_residues(30);   # M9

# M7: gene X (one exon) then gene Y (one exon), unrelated proteins
my $px = random_protein(200); my $py = random_protein(180);
my $w = random_dna(3000);
my $x_start = length($w) + 1; $w .= encode($px) . 'TAA'; my $x_end = length $w;
$w .= random_dna(1500);
my $y_start = length($w) + 1; $w .= encode($py) . 'TAG'; my $y_end = length $w;
$w .= random_dna(3000);
$chr{chrW} = $w;

my $fa = '';
for my $c (sort keys %chr) { $fa .= ">$c\n" . join("\n", $chr{$c} =~ /(.{1,60})/g) . "\n" }
write_file("$dir/genome.fa", $fa);
system("samtools faidx $dir/genome.fa") == 0 or die "samtools faidx failed\n";

# annotation: A = exons 1-2 with a stop codon read on into intron 2 (the annotation's error);
# B = exons 3-4. Coordinates of the stop: first 3 bases of intron 2 are GT? -> use the TAA we place.
# To give A a real stop, A's last CDS ends at exon 2's end and the next 3 bases (intron start "GT.")
# are not a stop; so A is given the in-frame stop by extending its CDS 3 bases into the intron only
# when those bases are a stop. Simplest honest model: A ends at exon 2 without a stop (partial model).
my @gff = ('##gff-version 3');
my $g = sub { my ($c, $id, $strand, @ex) = @_;
  my ($s, $e) = ($ex[0][0], $ex[-1][1]);
  push @gff, join("\t", $c, 'maker', 'gene', $s, $e, '.', $strand, '.', "ID=$id"),
             join("\t", $c, 'maker', 'mRNA', $s, $e, '.', $strand, '.', "ID=$id-RA;Parent=$id");
  my $coding = 0; my $n = 0;
  for my $x (@ex) { $n++;
    push @gff, join("\t", $c, 'maker', 'exon', @$x, '.', $strand, '.', "ID=$id-RA:e$n;Parent=$id-RA");
    push @gff, join("\t", $c, 'maker', 'CDS', @$x, '.', $strand, (3 - $coding % 3) % 3, "ID=$id-RA:c$n;Parent=$id-RA");
    $coding += $x->[1] - $x->[0] + 1 }
};
# fragment A is also wrong at the break: its exon 2 ends 60 nt early (20 codons lost), as a
# predictor's error there would be; only a rebuild from the protein recovers them
$g->('chrS', 'A', '+', $exons[0], [$exons[1][0], $exons[1][1] - 60]);
$g->('chrS', 'B', '+', @exons[2, 3]);
$g->('chrW', 'X', '+', [$x_start, $x_end]);
$g->('chrW', 'Y', '+', [$y_start, $y_end]);
write_file("$dir/in.gff", join("\n", @gff) . "\n");
write_file("$dir/ref.fa", ">REFP\n$protein\n>REFND\n$ref_ntdiff\n>REFCD\n$ref_ctdiff\n>REFX\n$px\n");
write_file("$dir/merges.txt", join("\t", qw(merge_id genes_in_order flag max_tiling_hits best_combined_cov_pct best_hit)) . "\n"
  . "merge_1\tA,B\tCLEAN\t3\t95\tREFP\nmerge_2\tX,Y\tCLEAN\t3\t95\tREFX\n");
write_file("$dir/merges_nd.txt", join("\t", qw(merge_id genes_in_order flag max_tiling_hits best_combined_cov_pct best_hit)) . "\n"
  . "merge_1\tA,B\tCLEAN\t3\t95\tREFND\n");
write_file("$dir/merges_cd.txt", join("\t", qw(merge_id genes_in_order flag max_tiling_hits best_combined_cov_pct best_hit)) . "\n"
  . "merge_1\tA,B\tCLEAN\t3\t95\tREFCD\n");
write_file("$dir/merges_multi.txt", join("\t", qw(merge_id genes_in_order flag max_tiling_hits best_combined_cov_pct best_hit ref_candidates)) . "\n"
  . "merge_1\tA,B\tCLEAN\t3\t95\tREFX\tREFX,REFP\nmerge_2\tX,Y\tCLEAN\t3\t95\tREFX\tREFX,REFP\n");

my $bin = "$FindBin::Bin/..";
system("perl $bin/merge_split_genes.pl --genome_fa $dir/genome.fa --removed $dir/removed.gff $dir/merges.txt $dir/in.gff $dir/merged.gff > $dir/merge.log 2>&1") == 0
  or die "merge failed: " . `tail -3 $dir/merge.log`;
my $rebuild = sub { my ($table, $tag, $extra) = @_;
  my $st = system("perl $bin/rebuild_with_miniprot.pl --merged_gff $dir/merged.gff --removed_gff $dir/removed.gff "
                . "--merge_table $table --ref_fa $dir/ref.fa --genome_fa $dir/genome.fa --out $dir/rebuilt_$tag.gff "
                . "--report $dir/mp_$tag.tsv --scratch $dir --flank 2000 " . ($extra // '') . " > $dir/rebuild_$tag.log 2>&1");
  check($st == 0, "rebuild_with_miniprot.pl runs ($tag)", `tail -3 $dir/rebuild_$tag.log`);
  my (%gene, %cds, %tx_merge);   # %cds: merge_id -> CDS segments
  open my $fh, '<', "$dir/rebuilt_$tag.gff" or return ({}, {});
  while (my $line = <$fh>) { chomp $line; my @f = split /\t/, $line; next unless @f >= 9 && $f[1] eq 'Mender';
    my ($m) = $f[8] =~ /merge_id=([^;]+)/;
    my ($id) = $f[8] =~ /ID=([^;]+)/;
    $gene{$m} = $f[8] if $f[2] eq 'gene';
    $tx_merge{$id} = $m if $f[2] eq 'mRNA';
    if ($f[2] eq 'CDS') { my ($p) = $f[8] =~ /Parent=([^;]+)/; push @{$cds{ $tx_merge{$p} // '' }}, [@f[0, 3, 4, 7]] } }
  close $fh;
  return (\%gene, \%cds);
};
my ($gene, $cds) = $rebuild->("$dir/merges.txt", 'main');

# M6: SUPPORTED, exons = the gene's, translation = the protein
my @segs = sort by_start @{ $cds->{merge_1} // [] };
my (@got, @want);
for my $seg (@segs) { push @got, "$seg->[1]-$seg->[2]" }
for my $exon (@exons) { push @want, "$exon->[0]-$exon->[1]" }
my $got_exons = join(' ', @got);
my $want_exons = join(' ', @want);
$want_exons =~ s/-(\d+)$/'-' . $1/e;
check(scalar(($gene->{merge_1} // '') =~ /miniprot_support=SUPPORTED/), 'M6: a split gene is SUPPORTED', $gene->{merge_1});
my $nt = cds_sequence(@segs);
(my $tr = translate($nt)) =~ s/\*$//;
check($tr eq $protein, 'M6: the rebuilt model translates to the reference protein', substr($tr, 0, 60) . '...');
check(@segs == 4 && $segs[0][1] == $exons[0][0] && $segs[1][1] == $exons[1][0] && $segs[1][2] == $exons[1][1]
      && $segs[2][1] == $exons[2][0] && $segs[3][1] == $exons[3][0],
      'M6: the rebuilt model has the gene\'s four exons, with its real splice sites', "$got_exons vs $want_exons");
check(($gene->{merge_1} // '') =~ /start_codon=present/ && ($gene->{merge_1} // '') =~ /stop_codon=present/,
      'M6: start and stop codons present', $gene->{merge_1});

# M7: NOT_SUPPORTED
check(scalar(($gene->{merge_2} // '') =~ /miniprot_support=NOT_SUPPORTED:.*covers_1_of_2_sources/),
      'M7: two different genes merged by mistake: NOT_SUPPORTED (covers 1 of 2)', $gene->{merge_2});

# M8: N-terminus differs -> start completed from fragment A
my ($gene_nd, $cds_nd) = $rebuild->("$dir/merges_nd.txt", 'nd');
check(($gene_nd->{merge_1} // '') =~ /miniprot_support=SUPPORTED/ && ($gene_nd->{merge_1} // '') =~ /start_codon=from_fragment:A/,
      'M8: miniprot starts downstream of the ATG; fragment A\'s start completes the model', $gene_nd->{merge_1});
my @s8 = sort by_start @{ $cds_nd->{merge_1} // [] };
my $nt8 = cds_sequence(@s8);
(my $tr8 = translate($nt8)) =~ s/\*$//;
check($tr8 eq $protein, 'M8: the completed model translates to the whole gene, from its M', substr($tr8, 0, 40) . '...');

# M9: C-terminus differs -> stop completed from fragment B
my ($gene_cd, $cds_cd) = $rebuild->("$dir/merges_cd.txt", 'cd');
check(scalar(($gene_cd->{merge_1} // '') =~ /stop_codon=from_fragment:B/), 'M9: miniprot ends before the stop; fragment B\'s stop completes the model', $gene_cd->{merge_1});
my @s9 = sort by_start @{ $cds_cd->{merge_1} // [] };
my $nt9 = cds_sequence(@s9);
my $tr9 = translate($nt9);
check($tr9 eq "$protein*", 'M9: the completed model translates to the whole gene and ends in its stop codon', substr($tr9, -40));

# M10-M12: several candidate references, tried best-ranked first
my ($gene_mu) = $rebuild->("$dir/merges_multi.txt", 'multi');
my $m10 = $gene_mu->{merge_1} // '';
check($m10 =~ /miniprot_support=SUPPORTED;/ && $m10 =~ /miniprot_ref=REFP;/ && $m10 =~ /miniprot_ref_rank=2;/ && $m10 =~ /miniprot_refs_tried=2;/,
      'M10: the first candidate fails, the second (the gene\'s own protein) SUPPORTS the merge', $m10);
my $m11 = $gene_mu->{merge_2} // '';
check($m11 =~ /miniprot_support=NOT_SUPPORTED:/ && $m11 =~ /miniprot_ref=REFX;/ && $m11 =~ /miniprot_refs_tried=2/,
      'M11: a false merge stays NOT_SUPPORTED with two candidates; the best-ranked one is reported', $m11);
my ($gene_one) = $rebuild->("$dir/merges_multi.txt", 'one', '--max_refs 1');
my $m12 = $gene_one->{merge_1} // '';
check($m12 =~ /miniprot_support=NOT_SUPPORTED:/ && $m12 =~ /miniprot_refs_tried=1/,
      'M12: --max_refs 1 tries only the first candidate', $m12);

print $failed ? "\n$failed FAILED, $passed passed\n" : "all $passed checks passed\n";
exit($failed ? 1 : 0);
