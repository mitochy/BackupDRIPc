#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;
my $type = $input1 =~ /orig/ ? "orig" : $input1 =~ /shuf/ ? "shuf" : die "Cannot determine type\n";

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
open (my $out1, ">", "$fileName1.TSS") or die "Cannot write to $fileName1.TSS: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
	my ($chr2, $start2, $end2) = $type eq "orig" ? $info =~ /ORIG=ENS\w+\.\d+,(chr\w+),(\d+),(\d+);TWIN=/ : $info =~ /TWIN=ENS\w+\.\d+,(chr\w+),(\d+),(\d+)$/;
	die "Died at $info\n" if not defined ($end2) or $end2 !~ /^\d+$/;
	($start, $end) = ($start2 - 3000, $end2 + 3000);
	print $out1 "$chr2\t$start\t$end\t$name\t$val\t$strand\t$info\n";
}
close $in1;
close $out1;
