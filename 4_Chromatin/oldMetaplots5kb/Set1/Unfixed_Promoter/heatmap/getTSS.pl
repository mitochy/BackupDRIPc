#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
$fileName1 =~ s/dripc/dripcTSS/;
open (my $out1, ">", "$fileName1.input") or die "Cannot write to $fileName1.input: $!\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
	my ($newChr, $newStart, $newEnd) = $info =~ /ORIG=E\w+\.\d+,(chr\w+),(\d+),(\d+);TWIN/ if $input1 =~ /orig/i;
	($newChr, $newStart, $newEnd) = $info =~ /TWIN=E\w+\.\d+,(chr\w+),(\d+),(\d+)$/ if $input1 =~ /shuf/i;
	die "Died at $line\n" if not defined($newStart);
	$newStart -= 3000;
	$newEnd   += 3000;
	print $out1 "$newChr\t$newStart\t$newEnd\t$name\t$val\t$strand\t$info\n";
}
close $in1;

close $out1;
