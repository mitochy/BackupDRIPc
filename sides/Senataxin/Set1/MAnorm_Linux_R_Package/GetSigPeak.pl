#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
open (my $out1, ">", "$fileName1\_increase.sigpeak") or die "Cannot write to $fileName1\_increase.sigpeak: $!\n";
open (my $out2, ">", "$fileName1\_decrease.sigpeak") or die "Cannot write to $fileName1\_decrease.sigpeak: $!\n";

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $raw1, $raw2, $M, $A, $P) = split("\t", $line);
	next if $line =~ /\tNA/;
	print $out1 "$chr\t$start\t$end\n" if ($M <= 0.5 and $A >= 2 and $P >= 1.3);
	print $out2 "$chr\t$start\t$end\n" if ($M >= 2 and $A >= 2 and $P >= 1.3);
}
close $in1;

close $out1;
