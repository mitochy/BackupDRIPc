#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand) = split("\t", $line);
	$end = $end - 1;
	print $out1 "$chr\t$start\t$end\t$name\t$val\t$strand\n";
}
close $in1;
close $out1;
