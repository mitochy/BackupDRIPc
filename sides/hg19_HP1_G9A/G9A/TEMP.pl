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
	my ($chr, $start, $end, $name, $val, $strand) = split("\t", $line);
	$end += 1;
	print $out "$chr\t$start\t$end\t$name\t$val\t$strand\n";
}

close $in;
close $out;
