#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";

while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name, $strand, $chr, $start, $mapq, $map, $start, $zero, $zero2, $seq, $qual) = split("\t", $line);
}

close $in;
close $out;
