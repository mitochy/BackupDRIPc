#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
my $inputrna = "../../data/NT2.rpkm";

open (my $in, "<", $inputrna) or die "Cannot read from $inputrna: $!\n";
my %rna;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($gene, $value) = split("\t", $line);
	next if $value < 20;
	$rna{$gene} = 1;
}

close $in;


open (my $in2, "<", $input) or die;
open (my $out, ">", "$input.out") or die "Cannot write to $input.out: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	my ($chr, $start, $end, $name, $value, $strand, $names) = split("\t", $line);
	my @names = split(";", $names);
	my $check = 0;
	foreach my $namez (@names) {
		$check = 1 if defined($rna{$namez});
	}
	next if $check == 0;
	print $out "$line\n";
}
close $in2;
close $out;
