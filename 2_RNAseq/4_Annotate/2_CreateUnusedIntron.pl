#!/usr/bin/perl

use strict; use warnings; use mitochy;


my $bedFile = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
my $intronFile = "hg19Gencode_intron.bed";
mitochy::runbash("bedtools_bed_change.pl -b -x 0 -y 5228 -i $bedFile -o terminalRNA.tmp");

open (my $in,  "<", "NT2_rep1_uniq_hg19Gencode_exon.rnacount");
open (my $in2, "<", "NT2_rep2_uniq_hg19Gencode_exon.rnacount");
my %exon;
while (my $line = <$in>) {
	chomp($line);
	my ($gene, $count) = split("\t", $line);
	my @gene = split(",", $gene);
	foreach my $gene (@gene) {
		$exon{$gene} += $count / 2;
	}
}
close $in;
while (my $line = <$in2>) {
	chomp($line);
	my ($gene, $count) = split("\t", $line);
	my @gene = split(",", $gene);
	foreach my $gene (@gene) {
		$exon{$gene} = $count / 2;
	}
}
close $in2;
open (my $in3, "<", "terminalRNA.tmp");
open (my $out, ">", "terminalRNA.bed");
while (my $line = <$in3>) {
	chomp($line);
	my ($chr, $start, $end, $gene, $zero, $strand) = split("\t", $line);
	if (defined($exon{$gene})) {
		print $out "$line\n" if $exon{$gene} > 40;
	}
}
close $in3;
close $out;
system("rm terminalRNA.tmp");

mitochy::runbash("bedtools intersect -u -a hg19Gencode_intron.bed -b terminalRNA.bed > hg19Gencode_BADINTRON.bed");
print "Output = hg19Gencode_BADINTRON.bed\n";
