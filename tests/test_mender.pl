#!/usr/bin/perl
use strict;
use warnings;
use File::Temp qw(tempdir);
use FindBin;

# tests/test_mender.pl -- hermetic tests of merge_split_genes.pl and find_split_genes.pl on a
# made-up genome whose gene fragments are designed codon by codon, so the right answer is known.
#   perl tests/test_mender.pl            (exit 0 = all pass)
#
# merge_split_genes.pl
#   M1  + strand, 3 exons per fragment, non-zero phases: stop trimmed, merged CDS translates to
#       exactly fragment A (no stop) + fragment B, and every phase obeys GFF3
#   M2  the same gene on the - strand
#   M3  a fragment that does NOT end in a stop keeps its last codon (no_stop_trimmed)
#   M4  a fragment whose coding length leaves the next one out of frame -> frame_break
#   M5  a stop codon split by an intron is trimmed across the two CDS segments
# find_split_genes.pl (synthetic DIAMOND table)
#   F1  + strand, g1 N-terminal / g2 C-terminal -> a merge candidate
#   F2  + strand, pieces in the wrong order on the reference -> not a candidate
#   F3  - strand, the upstream gene is the C-terminal piece -> a candidate
#   F4  a "fragment" aligning over only 30% of itself (domain sharing) -> not a candidate
#   F5  20 genes inside a long gene are excluded from the rank (a 10-gene look-back missed them)
#   F6  identical output under two Perl hash seeds (ties never depend on hash order)
#   F7  a piece extending 40 aa past its alignment at the junction -> JUNCTION_OVERHANG

my $dir = tempdir(CLEANUP => 1);
my ($passed, $failed) = (0, 0);
# A regex match passed straight to check() is wrapped in scalar(): in list context a FAILED match
# returns an empty list, which shifts the arguments and makes the description the (true) result.
sub check {
  my ($ok, $what, $detail) = @_;
  if ($ok) { $passed++ } else { $failed++; print "FAIL: $what\n" . (defined $detail ? "      got: $detail\n" : '') }
}

my %CODON = (
  TTT=>'F',TTC=>'F',TTA=>'L',TTG=>'L',CTT=>'L',CTC=>'L',CTA=>'L',CTG=>'L',ATT=>'I',ATC=>'I',ATA=>'I',ATG=>'M',
  GTT=>'V',GTC=>'V',GTA=>'V',GTG=>'V',TCT=>'S',TCC=>'S',TCA=>'S',TCG=>'S',CCT=>'P',CCC=>'P',CCA=>'P',CCG=>'P',
  ACT=>'T',ACC=>'T',ACA=>'T',ACG=>'T',GCT=>'A',GCC=>'A',GCA=>'A',GCG=>'A',TAT=>'Y',TAC=>'Y',TAA=>'*',TAG=>'*',
  CAT=>'H',CAC=>'H',CAA=>'Q',CAG=>'Q',AAT=>'N',AAC=>'N',AAA=>'K',AAG=>'K',GAT=>'D',GAC=>'D',GAA=>'E',GAG=>'E',
  TGT=>'C',TGC=>'C',TGA=>'*',TGG=>'W',CGT=>'R',CGC=>'R',CGA=>'R',CGG=>'R',AGT=>'S',AGC=>'S',AGA=>'R',AGG=>'R',
  GGT=>'G',GGC=>'G',GGA=>'G',GGG=>'G');
sub translate { my ($s) = @_; my $p = ''; for (my $i = 0; $i + 3 <= length $s; $i += 3) { $p .= $CODON{substr $s, $i, 3} // 'X' } return $p }
sub revcomp { my ($s) = @_; ($s = reverse $s) =~ tr/ACGT/TGCA/; return $s }
sub write_file { my ($path, $text) = @_; open my $fh, '>', $path or die "$path: $!"; print $fh $text; close $fh }

# ---- a gene built from exons: returns (gene coordinates, GFF lines for gene/mRNA/exon/CDS)
# $cds_pieces: coding sequence split into exon pieces (5'->3' in transcription order).
# Built on + at $pos; for - the whole locus is reverse-complemented into place by the caller.
my $INTRON = 'GT' . ('C' x 36) . 'AG';   # 40 nt, no stop in any frame of the Cs
sub gene_on_plus {
  my ($id, $pos, @pieces) = @_;
  my ($seq, @exons) = ('');
  my $start = $pos;
  for my $i (0 .. $#pieces) {
    push @exons, [$pos + length($seq), $pos + length($seq) + length($pieces[$i]) - 1];
    $seq .= $pieces[$i];
    $seq .= $INTRON if $i < $#pieces;
  }
  return ($seq, \@exons);
}
sub gff_block {
  my ($chr, $id, $strand, $exons, $phases) = @_;
  my ($s, $e) = ($exons->[0][0], $exons->[-1][1]);
  ($s, $e) = ($e, $s) if $s > $e;
  my @lines = (join("\t", $chr, 'test', 'gene', $s, $e, '.', $strand, '.', "ID=$id"),
               join("\t", $chr, 'test', 'mRNA', $s, $e, '.', $strand, '.', "ID=$id-RA;Parent=$id"));
  for my $i (0 .. $#$exons) {
    my ($es, $ee) = sort { $a <=> $b } @{$exons->[$i]};
    push @lines, join("\t", $chr, 'test', 'exon', $es, $ee, '.', $strand, '.', "ID=$id-RA:e$i;Parent=$id-RA");
    push @lines, join("\t", $chr, 'test', 'CDS',  $es, $ee, '.', $strand, $phases->[$i], "ID=$id-RA:c$i;Parent=$id-RA");
  }
  return @lines;
}
# GFF3 phases from piece lengths and the first phase
sub phases_for { my ($first, @lens) = @_; my ($c, @p) = (-$first); for my $l (@lens) { push @p, (3 - $c % 3) % 3; $c += $l } return @p }

# ---- fragment A: M + 20 residues + stop, over 3 exons of 20, 25, 22 nt (66 = 22 codons incl. stop)
my $A_cds = 'ATG' . ('GCT' x 20) . 'TAA';          # MAAAA...A*
my @A_pieces = (substr($A_cds, 0, 20), substr($A_cds, 20, 25), substr($A_cds, 45));
# ---- fragment B: M + 15 residues + stop, over 3 exons of 17, 13, 18 nt
my $B_cds = 'ATG' . ('GGT' x 15) . 'TAG';          # MGGG...G*
my @B_pieces = (substr($B_cds, 0, 17), substr($B_cds, 17, 13), substr($B_cds, 30));
my $expected_protein = translate(substr($A_cds, 0, -3) . $B_cds);   # A without its stop, then B
$expected_protein =~ s/\*$//;

# locus on + : filler, A, spacer, B, filler
my $filler = 'C' x 100;
my ($A_seq, $A_ex) = gene_on_plus('A', 101, @A_pieces);
my $B_start = 101 + length($A_seq) + 60;
my ($B_seq, $B_ex) = gene_on_plus('B', $B_start, @B_pieces);
my $plus = $filler . $A_seq . ('C' x 60) . $B_seq . $filler;

# chrM: the same locus reverse-complemented (so A is still 5' of B in transcription order)
my $L = length $plus;
my $minus = revcomp($plus);
my $flip = sub { my ($ex) = @_; return [ map { [ $L - $_->[0] + 1, $L - $_->[1] + 1 ] } @$ex ] };

# chrN: fragment A without a stop (ends in a sense codon GCT), then B
my $An_cds = 'ATG' . ('GCT' x 21);
my @An_pieces = (substr($An_cds, 0, 20), substr($An_cds, 20, 25), substr($An_cds, 45));
my ($An_seq, $An_ex) = gene_on_plus('An', 101, @An_pieces);
my ($Bn_seq, $Bn_ex) = gene_on_plus('Bn', 101 + length($An_seq) + 60, @B_pieces);
my $chrN = $filler . $An_seq . ('C' x 60) . $Bn_seq . $filler;

# chrF: fragment A whose coding length (before its stop) is not a multiple of 3 -> frame break
my $Af_cds = 'ATG' . ('GCT' x 20) . 'GC' . 'TAA';
my @Af_pieces = (substr($Af_cds, 0, 20), substr($Af_cds, 20, 25), substr($Af_cds, 45));
my ($Af_seq, $Af_ex) = gene_on_plus('Af', 101, @Af_pieces);
my ($Bf_seq, $Bf_ex) = gene_on_plus('Bf', 101 + length($Af_seq) + 60, @B_pieces);
my $chrF = $filler . $Af_seq . ('C' x 60) . $Bf_seq . $filler;

# chrS: fragment A whose stop codon TAA is split by an intron (last exon = "A" only)
my @As_pieces = (substr($A_cds, 0, 30), substr($A_cds, 30, 35), substr($A_cds, 65));   # 30, 35, 1 nt
my ($As_seq, $As_ex) = gene_on_plus('As', 101, @As_pieces);
my ($Bs_seq, $Bs_ex) = gene_on_plus('Bs', 101 + length($As_seq) + 60, @B_pieces);
my $chrS = $filler . $As_seq . ('C' x 60) . $Bs_seq . $filler;

my %chr = (chrP => $plus, chrM => $minus, chrN => $chrN, chrF => $chrF, chrS => $chrS);
my ($fa, $fai, $offset) = ('', '', 0);
for my $name (sort keys %chr) {
  my $header = ">$name\n";
  my $seq = $chr{$name};
  my @lines = $seq =~ /(.{1,60})/g;
  $fai .= join("\t", $name, length $seq, $offset + length $header, 60, 61) . "\n";
  my $record = $header . join("\n", @lines) . "\n";
  $fa .= $record;
  $offset += length $record;
}
write_file("$dir/genome.fa", $fa);
write_file("$dir/genome.fa.fai", $fai);

my @gff = ('##gff-version 3');
push @gff, gff_block('chrP', 'A', '+', $A_ex, [phases_for(0, map { length } @A_pieces)]);
push @gff, gff_block('chrP', 'B', '+', $B_ex, [phases_for(0, map { length } @B_pieces)]);
# on chrM the upstream gene in genome coordinates is B; transcription order still A then B
push @gff, gff_block('chrM', 'MB', '-', $flip->($B_ex), [phases_for(0, map { length } @B_pieces)]);
push @gff, gff_block('chrM', 'MA', '-', $flip->($A_ex), [phases_for(0, map { length } @A_pieces)]);
push @gff, gff_block('chrN', 'An', '+', $An_ex, [phases_for(0, map { length } @An_pieces)]);
push @gff, gff_block('chrN', 'Bn', '+', $Bn_ex, [phases_for(0, map { length } @B_pieces)]);
push @gff, gff_block('chrF', 'Af', '+', $Af_ex, [phases_for(0, map { length } @Af_pieces)]);
push @gff, gff_block('chrF', 'Bf', '+', $Bf_ex, [phases_for(0, map { length } @B_pieces)]);
push @gff, gff_block('chrS', 'As', '+', $As_ex, [phases_for(0, map { length } @As_pieces)]);
push @gff, gff_block('chrS', 'Bs', '+', $Bs_ex, [phases_for(0, map { length } @B_pieces)]);
write_file("$dir/in.gff", join("\n", @gff) . "\n");
write_file("$dir/merges.txt", join("\t", qw(merge_id genes_in_order flag max_tiling_hits best_combined_cov_pct)) . "\n"
  . "merge_1\tA,B\tCLEAN\t3\t95\nmerge_2\tMB,MA\tCLEAN\t3\t95\nmerge_3\tAn,Bn\tCLEAN\t3\t95\n"
  . "merge_4\tAf,Bf\tCLEAN\t3\t95\nmerge_5\tAs,Bs\tCLEAN\t3\t95\n");

my $merge = "$FindBin::Bin/../merge_split_genes.pl";
my $status = system("perl \Q$merge\E --genome_fa \Q$dir/genome.fa\E --removed \Q$dir/removed.gff\E "
                  . "\Q$dir/merges.txt\E \Q$dir/in.gff\E \Q$dir/out.gff\E > \Q$dir/merge.log\E 2>&1");
check($status == 0, 'merge_split_genes.pl runs', `tail -3 \Q$dir/merge.log\E`);

# ---- read the merged genes back: per merge_id, the gene line, and CDS lines
my (%gene_line, %cds);
open my $out, '<', "$dir/out.gff" or die;
my %tid_merge;
while (my $line = <$out>) {
  chomp $line;
  next if $line =~ /^#/ or $line !~ /\t/;
  my @f = split /\t/, $line;
  next unless $f[1] eq 'Mender';
  my ($id) = $f[8] =~ /ID=([^;]+)/;
  my ($parent) = $f[8] =~ /Parent=([^;]+)/;
  if ($f[2] eq 'gene') { my ($m) = $f[8] =~ /merge_id=([^;]+)/; $gene_line{$m} = $line; }
  elsif ($f[2] eq 'mRNA') { my ($m) = $f[8] =~ /merge_id=([^;]+)/; $tid_merge{$id} = $m; }
  elsif ($f[2] eq 'CDS') { push @{$cds{$tid_merge{$parent}}}, \@f; }
}
close $out;

# merged protein and phase check for one merge
sub merged {
  my ($merge_id) = @_;
  my @segs = @{$cds{$merge_id} // []};
  return ('', 'no CDS') unless @segs;
  my $strand = $segs[0][6];
  @segs = $strand eq '-' ? sort { $b->[4] <=> $a->[4] } @segs : sort { $a->[3] <=> $b->[3] } @segs;
  my ($nt, @got_phase, @lens) = ('');
  for my $s (@segs) {
    my $piece = substr($chr{$s->[0]}, $s->[3] - 1, $s->[4] - $s->[3] + 1);
    $nt .= $strand eq '-' ? revcomp($piece) : $piece;
    push @got_phase, $s->[7];
    push @lens, $s->[4] - $s->[3] + 1;
  }
  my @want = phases_for($got_phase[0], @lens);
  (my $protein = translate(substr($nt, $got_phase[0]))) =~ s/\*$//;
  return ($protein, "@got_phase" eq "@want" ? 'ok' : "phases @got_phase, GFF3 wants @want");
}

my ($p1, $ph1) = merged('merge_1');
check($p1 eq $expected_protein, 'M1: + strand merged CDS = A (stop removed) + B, no internal stop', $p1);
check($ph1 eq 'ok', 'M1: + strand phases obey GFF3 across 6 CDS segments', $ph1);
my ($p2, $ph2) = merged('merge_2');
check($p2 eq $expected_protein, 'M2: - strand merged CDS = A (stop removed) + B', $p2);
check($ph2 eq 'ok', 'M2: - strand phases obey GFF3', $ph2);
my ($p3) = merged('merge_3');
check(($gene_line{merge_3} // '') =~ /no_stop_trimmed=An-RA/ && $p3 eq translate($An_cds) . substr($expected_protein, 21),
      'M3: a fragment without a stop keeps its last codon (no_stop_trimmed)', $p3);
check(scalar(($gene_line{merge_4} // '') =~ /frame_break=Af-RA>Bf-RA/), 'M4: coding length leaves B out of frame -> frame_break', $gene_line{merge_4});
check(($gene_line{merge_1} // '') !~ /frame_break|no_stop_trimmed/, 'M1: an in-frame merge carries no frame_break / no_stop_trimmed', $gene_line{merge_1});
my ($p5, $ph5) = merged('merge_5');
check($p5 eq $expected_protein && $ph5 eq 'ok', 'M5: a stop codon split by an intron is trimmed across segments', "$p5 / $ph5");

# ================================================================ find_split_genes.pl
# genes on chr1 (+): g1 g2 adjacent; chr6: g3 g4 adjacent (wrong order on the reference); chr2 (-): h1 h2;
# chr3: d1 d2 domain-sharing (qcov 30); chr4: k1, a long gene L containing 20 genes, k2;
# chr5: o1 o2 with a 40-aa overhang at the junction
my @fgff = ('##gff-version 3');
my $gene = sub { my ($chr, $id, $s, $e, $strand) = @_;
  push @fgff, join("\t", $chr, 't', 'gene', $s, $e, '.', $strand, '.', "ID=$id"),
              join("\t", $chr, 't', 'mRNA', $s, $e, '.', $strand, '.', "ID=$id.t;Parent=$id") };
$gene->('chr1', 'g1', 1000, 2000, '+'); $gene->('chr1', 'g2', 3000, 4000, '+');
$gene->('chr6', 'g3', 10000, 11000, '+'); $gene->('chr6', 'g4', 12000, 13000, '+');
$gene->('chr2', 'h1', 1000, 2000, '-'); $gene->('chr2', 'h2', 3000, 4000, '-');
$gene->('chr3', 'd1', 1000, 2000, '+'); $gene->('chr3', 'd2', 3000, 4000, '+');
$gene->('chr4', 'k1', 500, 900, '+'); $gene->('chr4', 'L', 1000, 900000, '+');
$gene->('chr4', "n$_", 1000 + $_ * 40000, 1000 + $_ * 40000 + 500, '+') for 1 .. 20;
$gene->('chr4', 'k2', 950000, 951000, '+');
$gene->('chr5', 'o1', 1000, 2000, '+'); $gene->('chr5', 'o2', 3000, 4000, '+');
write_file("$dir/f.gff", join("\n", @fgff) . "\n");
# subject: REF (500 aa) and REF2, a second reference with the same tiling (a tie in coverage)
write_file("$dir/subject.fa", ">REF gene:ENSG1 gene_symbol:REF1 description:reference one\nM\n>REF2 gene:ENSG2 gene_symbol:REF2 description:reference two\nM\n");
my @ids = qw(g1 g2 g3 g4 h1 h2 d1 d2 k1 k2 L o1 o2);
write_file("$dir/query.fa", join('', map { ">$_.t desc $_\nM\n" } @ids));
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore scovhsp slen qcovhsp qlen
my $row = sub { my ($q, $s, $ss, $se, $qs, $qe, $qlen) = @_;
  my $scov = int(100 * ($se - $ss + 1) / 500); my $qcov = int(100 * ($qe - $qs + 1) / $qlen);
  return join("\t", "$q.t", $s, 50, $se - $ss + 1, 0, 0, $qs, $qe, $ss, $se, '1e-50', 300, $scov, 500, $qcov, $qlen) . "\n" };
my $blast = join('',
  (map { $row->('g1', $_, 1, 240, 1, 240, 245), $row->('g2', $_, 245, 495, 1, 251, 255) } qw(REF REF2)),   # F1 (+ two tied refs)
  $row->('g3', 'REF', 245, 495, 1, 251, 255), $row->('g4', 'REF', 1, 240, 1, 240, 245),   # F2 wrong order
  $row->('h1', 'REF', 245, 495, 1, 251, 255), $row->('h2', 'REF', 1, 240, 1, 240, 245),   # F3 - strand: h1 (upstream) C-terminal
  $row->('d1', 'REF', 1, 240, 1, 240, 800), $row->('d2', 'REF', 245, 495, 1, 251, 850),   # F4 qcov 30%
  $row->('k1', 'REF', 1, 240, 1, 240, 245), $row->('k2', 'REF', 245, 495, 1, 251, 255),   # F5 around L
  $row->('o1', 'REF', 1, 240, 1, 240, 285), $row->('o2', 'REF', 245, 495, 1, 251, 255),   # F7 o1: 45 aa past its alignment
);
write_file("$dir/blast.tsv", $blast);
my $find = "$FindBin::Bin/../find_split_genes.pl";
my @runs;
for my $seed (1, 2) {
  my $run_dir = "$dir/find$seed";
  mkdir $run_dir;
  local $ENV{PERL_HASH_SEED} = $seed;
  my $st = system("cd \Q$run_dir\E && perl \Q$find\E ../blast.tsv ../f.gff ../subject.fa ../query.fa > find.log 2>&1");
  check($st == 0, "find_split_genes.pl runs (seed $seed)", `tail -3 \Q$run_dir/find.log\E`);
  push @runs, $run_dir;
}
my %cand;
open my $m, '<', "$runs[0]/split_genes_merge.txt" or die;
<$m>;
while (my $line = <$m>) { my @f = split /\t/, $line; $cand{$f[2]} = \@f; }
close $m;
check(exists $cand{'g1,g2'}, 'F1: + strand, N-terminal then C-terminal piece -> a candidate', join(' ', sort keys %cand));
check(!exists $cand{'g3,g4'}, 'F2: + strand, pieces in the wrong order on the reference -> not a candidate');
check(exists $cand{'h1,h2'}, 'F3: - strand, upstream gene is the C-terminal piece -> a candidate');
check(!exists $cand{'d1,d2'}, 'F4: pieces aligning over 30% of themselves (domain sharing) -> not a candidate');
check(exists $cand{'k1,k2'} && ($cand{'k1,k2'}[12] // '') eq '2',
      'F5: 20 genes inside the long gene L are excluded from the rank: k1 .. L .. k2 is dist 2, not 22',
      $cand{'k1,k2'} ? "dist $cand{'k1,k2'}[12]" : 'no candidate');
check(($cand{'g1,g2'}[6] // '') eq 'REF', 'F6: a coverage tie between two references goes to the lower id (REF)', $cand{'g1,g2'}[6]);
check(`cat \Q$runs[0]\E/split_genes_*.txt` eq `cat \Q$runs[1]\E/split_genes_*.txt`, 'F6: identical output under two hash seeds');
check(scalar(($cand{'o1,o2'}[16] // '') =~ /JUNCTION_OVERHANG/), 'F7: a piece 45 aa past its alignment at the junction -> JUNCTION_OVERHANG', $cand{'o1,o2'}[16]);

print $failed ? "\n$failed FAILED, $passed passed\n" : "all $passed checks passed\n";
exit($failed ? 1 : 0);
