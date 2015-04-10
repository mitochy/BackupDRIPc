#!/usr/bin/perl
# This script combine rnacount of genes from HTseq counting:
# Before:
# GeneA	155
# GeneA;GeneB	250
# 
# After:
# GeneA 405
# GeneB	250
#
# Then find RPKM (10E9*C)/(N*L)
# Where C = Number of reads mapped to a gene, which is just COUNT
# L = exon length in base-pairs for a gene, found from exon bed
# N = Total mapped reads ONLY TO GENES, which is just total read - total no feature (below)
#   24947538 NT2_rep1_uniq.sam
#   24535583 NT2_rep2_uniq.sam
#   7730539 NT2_rep1_uniq.reads_nofeature.sam
#   7607026 NT2_rep2_uniq.reads_nofeature.sam

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <.rnacount>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");

# Find N
my $N = $input =~ /rep1/ ? 24947538-7730539 : 24535583-7607026;

# Process exon bed file to get length (L)
print "1. Processing exon bed file ../1_*/hg19_gencode19_exon.bed file\n";
my %data;
open (my $in1, "<", "../1_DefineExonIntron/hg19_gencode19_exon.bed") or die "Cannot read from ../1_DefineExonIntron/hg19_gencode19_exon.bed: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	my ($chr, $start, $end, $genes, $zero, $strand) = split("\t", $line);
	my $length = $end - $start;
	my @genes = split(";", $genes);
	foreach my $gene (@genes) {
		$data{$gene}{length} += $length;
	}
}
close $in1;

# Process count file to get count and rpkm
print "2. Processing count file $input\n";
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($genes, $count) = split("\t", $line);
	my @genes = split(";", $genes);
	foreach my $gene (@genes) {
		$data{$gene}{count} += $count;
	}	
}
close $in;

open (my $out, ">", "$name.rpkm") or die "Cannot write to $name.rpkm: $!\n";
print "3. Printing RPKM\n";
foreach my $gene (sort keys %data) {
	my $rpkm = (10E9 * $data{$gene}{count}) / ($data{$gene}{length} * $N);
	print $out "$gene\t$rpkm\n";
}
close $out;

print "Output: $name.rpkm\n";
