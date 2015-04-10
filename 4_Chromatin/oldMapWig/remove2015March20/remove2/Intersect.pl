#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

mkdir "Intersect" if not -d "Intersect";
my @bed = <../Chromatin/peak/*.bed>;
my @feature = qw(promoter terminal genebody);

open (my $out, ">", "Intersect/Result.txt") or die;
print $out "Feature\tChIP\tTotalOrig\tTotalShuf\tTotalGenomic\tTotalChIP\tFeatOrig\tFeatShuf\tFeatChIP\tFeatChIPOrig\tFeatChIPShuf\tPercChIPOrig\tPercChIPShuf\tEnrichment\n";
foreach my $feature (@feature) {
	my $orig = "MappedOriginal/dripc_$feature.bed";
	my $shuf = "MappedShuffled/dripc_$feature.shuffled";
	my $genomic = "../../1_Location/genomic_$feature.bed";
	`bedtools intersect -u -a $genomic -b $orig > Intersect/$feature\_orig.bed` if not -e "Intersect/$feature\_orig.bed";
	`bedtools intersect -u -a $genomic -b $shuf > Intersect/$feature\_shuf.bed` if not -e "Intersect/$feature\_shuf.bed";
	my ($TotalOrig) = `wc -l $orig` =~ /^(\d+) /;
	my ($TotalShuf) = `wc -l $shuf` =~ /^(\d+) /;
	my ($TotalGenomic) = `wc -l $genomic` =~ /^(\d+) /;
	my ($FeatOrig) = `wc -l Intersect/$feature\_orig.bed` =~ /^(\d+) /;
	my ($FeatShuf) = `wc -l Intersect/$feature\_shuf.bed` =~ /^(\d+) /;
	foreach my $bed (@bed) {
		my ($chip) = mitochy::getFilename($bed);
		print "Processing $feature $chip\n";
		my ($TotalChIP) = `wc -l $bed` =~ /^(\d+) /;
		`bedtools intersect -u -a $genomic -b $bed > Intersect/$feature\_$chip.bed` if not -e "Intersect/$feature\_$chip.bed";
		`bedtools intersect -u -a Intersect/$feature\_orig.bed -b Intersect/$feature\_$chip.bed > Intersect/$feature\_$chip\_orig.bed` if not -e "Intersect/$feature\_$chip\_orig.bed";
		`bedtools intersect -u -a Intersect/$feature\_shuf.bed -b Intersect/$feature\_$chip.bed > Intersect/$feature\_$chip\_shuf.bed` if not -e "Intersect/$feature\_$chip\_shuf.bed";
		my ($FeatChIP) = `wc -l Intersect/$feature\_$chip.bed` =~ /^(\d+) /;
		my ($FeatChIPOrig) = `wc -l Intersect/$feature\_$chip\_orig.bed` =~ /^(\d+) /;
		my ($FeatChIPShuf) = `wc -l Intersect/$feature\_$chip\_shuf.bed` =~ /^(\d+) /;
		my $percOrig = int(($FeatChIPOrig+1)/($FeatOrig+1)*1000)/10;
		my $percShuf = int(($FeatChIPShuf+1)/($FeatShuf+1)*1000)/10;
		my $enrichment = int((($FeatChIPOrig+1)/($FeatOrig+1)) / (($FeatChIPShuf+1)/($FeatShuf+1)) * 100) / 100;
		print $out "$feature\t$chip\t$TotalOrig\t$TotalShuf\t$TotalGenomic\t$TotalChIP\t$FeatOrig\t$FeatShuf\t$FeatChIP\t$FeatChIPOrig\t$FeatChIPShuf\t$percOrig\t$percShuf\t$enrichment\n";
	}
}
