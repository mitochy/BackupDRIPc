#!/usr/bin/perl

use strict; use warnings; use mitochy;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

# Parse RNA-seq
my %rna;
open (my $inRNA, "<", "../data/NT2.rpkm") or die "Cannot read from ../data/NT2.rpkm: $!\n";
while (my $line = <$inRNA>) {
        chomp($line);
        next if $line =~ /#/;
        my ($gene, $value) = split("\t", $line);
        $rna{$gene} = $value;
}
close $inRNA;

my $inputtwin = "NT2_RNAtwin.rpkm";
open (my $in1, "<", $inputtwin) or die "Cannot read from $inputtwin: $!\n";
my $last_index = "INIT";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($gene, @twin) = split("\t", $line);
	
	# If number of gene twin is less than 20, don't use it
	print GREEN "\t$gene not used\n" and next if @twin < 20;

	# Check if twin RNA is within 20% of gene's
	my $rna = $rna{$gene};
	foreach my $twin (@twin) {
		my $rna2 = $rna{$twin};
		print GREEN "\t$twin ($rna2) is less than $gene ($rna)\n" if $rna2 < 0.8*$rna and $rna > 1.2*$rna;
	}


	
	# Get index based on gene rpkm
	my $index;
	$index = 10 + int($rna / 10)    if $rna >= 10   and $rna < 100;
	$index = 20 + int($rna / 100)   if $rna >= 100  and $rna < 1000;
	$index = 30 + int($rna / 1000)  if $rna >= 1000 and $rna < 10000;
	$index = 40 + int($rna / 10000) if $rna >= 10000;
	print "\tProcessing index $index\n" if $last_index ne $index;
	$last_index = $index;
}

