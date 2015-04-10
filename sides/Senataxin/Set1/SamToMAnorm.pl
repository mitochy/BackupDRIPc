#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name, $flag, $chr, $start, $mapq, $mapped) = split("\t", $line);
	next if $flag != 0 and $flag != 16;
	my @mapped = split("[A-Z]", $mapped);
	my $end = 0;
	foreach my $int (@mapped) {
		next if $mapped =~ /d/i;
		($end) += $mapped =~ /^(\d+)/;
	}
	$end = 50 if $end > 50;
	die "Died at $line\n" if not defined($end);
	$end += $start;
	my $strand = $flag == 0 ? "+" : "-";
	print $out1 "$chr\t$start\t$end\t$strand\n";
}
close $in1;

close $out1;
