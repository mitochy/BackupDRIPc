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
	my ($gene, $val) = split("\t", $line);
	$val = int($val * 100);
	print $out "$gene\t$val\n";
}

close $in;
close $out;
