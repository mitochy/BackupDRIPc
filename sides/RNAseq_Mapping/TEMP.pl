#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";

while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $val) = split("\t", $line);
	my $newend = $end - 1;
	next if $newend == $start or $newend == $start - 1;
	print $out "$chr\t$start\t$newend\tNAME\t$val\t\+\n";
}

close $in;
close $out;
