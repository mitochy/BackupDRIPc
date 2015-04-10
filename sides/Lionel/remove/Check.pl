#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my %good;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $zero, $strand, $chr2, $start2, $end2, $name2, $zero2, $strand2) = split("\t", $line);
	my $first = "$chr\t$start\t$end\t$name\t$zero\t$strand";
	my $second = "$chr2\t$start2\t$end2\t$name2\t$zero2\t$strand2";
	if (not defined($good{$first})) {
		$good{$first} = 1 if $first eq $second;
		$good{$first} = 0 if $first ne $second;
	}
}
close $in1;

my %data;
open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
foreach my $first (sort keys %good) {
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $first);
	$data{$chr}{$start}

}
close $out1;
