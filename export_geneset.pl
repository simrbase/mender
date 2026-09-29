#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Copy qw(copy);

## export_geneset.pl
##
## Step 10 of the Mender pipeline: writes the final annotation as a gene set that can be
## dropped into the genome datastore (genomes/v2/<species>/<assembly>/<geneset>/):
##
##   genes.gff            the annotation (copy of --gff, normally validated.gff)
##   genes.gtf            the same as GTF (gffread -T)
##   protein.aa.fa        proteins            (gffread -y)
##   cds.nt.fa            coding sequences    (gffread -x)
##   transcript.nt.fa     transcripts, UTRs included (gffread -w)
##   protein2gene.txt     sequence id <TAB> gene id, one line per sequence in the FASTA
##   cds2gene.txt
##   transcript2gene.txt
##
## Sequence ids are transcript ids (gffread names every sequence by its transcript), so
## the three maps list the same ids for coding transcripts; each lists exactly the ids
## in its own FASTA. metadata.yaml and DB/ are not written (they are the datastore's).
##
## FASTA descriptions (the text after the id):
##   genes Mender did not touch   the description of that id in --proteome_fa (the input
##                                annotation's own header text), unchanged
##   Mender genes                 merge_id=, merged_from= (source genes), model=miniprot
##                                (rebuilt) or placeholder, rebuild_ref= (the reference
##                                the model was built from), start_codon=/stop_codon=,
##                                utr5=/utr3=
##   transcript.nt.fa also keeps gffread's CDS=<start>-<end> (the coding part of the
##   transcript), before the description.
##
## Usage:
##   perl export_geneset.pl --gff validated.gff --genome_fa genome.fa --out_dir geneset \
##        [--proteome_fa input_proteins.fa] [--gffread_bin gffread]

my ($gff, $genome_fa, $out_dir, $proteome_fa);
my $gffread_bin = "gffread";
GetOptions(
    "gff=s"         => \$gff,
    "genome_fa=s"   => \$genome_fa,
    "out_dir=s"     => \$out_dir,
    "proteome_fa=s" => \$proteome_fa,
    "gffread_bin=s" => \$gffread_bin,
) or die "Error in options\n";
for my $required ([gff => $gff], [genome_fa => $genome_fa], [out_dir => $out_dir]) {
    die "--$required->[0] is required\n" unless defined $required->[1];
}
-d $out_dir or mkdir $out_dir or die "cant create $out_dir: $!\n";

# genes.gff, genes.gtf and the three raw FASTAs
copy($gff, "$out_dir/genes.gff") or die "cant copy $gff to $out_dir/genes.gff: $!\n";
my %raw = (protein => "$out_dir/.protein.raw.fa", cds => "$out_dir/.cds.raw.fa", transcript => "$out_dir/.transcript.raw.fa");
system("$gffread_bin $out_dir/genes.gff -g $genome_fa -y $raw{protein} -x $raw{cds} -w $raw{transcript}") == 0
    or die "gffread FASTA export failed\n";
system("$gffread_bin $out_dir/genes.gff -T -o $out_dir/genes.gtf") == 0 or die "gffread GTF export failed\n";

# transcript -> gene, and the Mender provenance of each Mender gene
my (%is_gene, %tx_gene, %gene_attrs, %mender_desc);
open my $in, "<", "$out_dir/genes.gff" or die "cant open $out_dir/genes.gff: $!\n";
while (my $line = <$in>) {
    next if $line =~ /^#/;
    chomp $line;
    my @f = split /\t/, $line;
    next unless @f >= 9;
    my %attr;
    for my $pair (split /;/, $f[8]) {
        my ($key, $value) = split /=/, $pair, 2;
        $attr{$key} = $value if defined $value;
    }
    next unless defined $attr{ID};
    if ($f[2] eq "gene") {
        $is_gene{ $attr{ID} } = 1;
        $gene_attrs{ $attr{ID} } = \%attr if ($attr{merge_source} // "") eq "Mender";
    } elsif (defined $attr{Parent} && $is_gene{ $attr{Parent} }) {
        $tx_gene{ $attr{ID} } = $attr{Parent};
    }
}
close $in;
for my $tx (keys %tx_gene) {
    my $g = $gene_attrs{ $tx_gene{$tx} } or next;
    my $rebuilt = ($g->{model_source} // "") eq "miniprot";
    my @parts = ("merge_id=" . ($g->{merge_id} // "."), "merged_from=" . ($g->{merged_from} // "."),
                 "model=" . ($rebuilt ? "miniprot" : "placeholder"));
    push @parts, "rebuild_ref=$g->{miniprot_ref}" if $rebuilt && defined $g->{miniprot_ref};
    for my $key (qw(start_codon stop_codon utr5 utr3)) { push @parts, "$key=$g->{$key}" if defined $g->{$key} }
    $mender_desc{$tx} = join(" ", @parts);
}

# the input annotation's own descriptions
my %input_desc;
if (defined $proteome_fa && $proteome_fa ne "") {
    open my $pf, "<", $proteome_fa or die "cant open $proteome_fa: $!\n";
    while (my $line = <$pf>) {
        next unless $line =~ /^>(\S+)[ \t]+(.*?)\s*$/;
        $input_desc{$1} = $2;
    }
    close $pf;
}

# final FASTAs and maps
my %final = (protein => "protein.aa.fa", cds => "cds.nt.fa", transcript => "transcript.nt.fa");
my %map   = (protein => "protein2gene.txt", cds => "cds2gene.txt", transcript => "transcript2gene.txt");
for my $kind (qw(protein cds transcript)) {
    open my $src, "<", $raw{$kind} or die "cant open $raw{$kind}: $!\n";
    open my $fa,  ">", "$out_dir/$final{$kind}" or die "cant write $out_dir/$final{$kind}: $!\n";
    open my $mp,  ">", "$out_dir/$map{$kind}" or die "cant write $out_dir/$map{$kind}: $!\n";
    my ($n, $n_mender, $n_nogene) = (0, 0, 0);
    while (my $line = <$src>) {
        unless ($line =~ /^>/) { print $fa $line; next }
        chomp $line;
        my ($id, $rest) = $line =~ /^>(\S+)\s*(.*)$/;
        my ($cds_range) = $rest =~ /\b(CDS=\d+-\d+)\b/;
        my $desc = $mender_desc{$id} // $input_desc{$id};
        $n_mender++ if defined $mender_desc{$id};
        my @header = (">$id");
        push @header, $cds_range if $kind eq "transcript" && defined $cds_range;
        push @header, $desc if defined $desc && $desc ne "";
        print $fa join(" ", @header), "\n";
        $n++;
        if (defined $tx_gene{$id}) { print $mp "$id\t$tx_gene{$id}\n" } else { $n_nogene++ }
    }
    close $src; close $fa; close $mp;
    unlink $raw{$kind};
    printf "  %-20s %6d sequences (%d Mender)  -> %s\n", $final{$kind}, $n, $n_mender, $map{$kind};
    print "  WARNING: $n_nogene $kind sequences have no parent gene in the GFF (not in $map{$kind})\n" if $n_nogene;
}
print "  genes.gff, genes.gtf written to $out_dir\n";
