#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my %both;
open (my $in0, "<", "BOTH") or die "Cannot read from BOTH: $!\n";
while (my $line = <$in0>) {
	chomp($line);
	next if $line =~ /#/;
	$both{$line} = 1;
}
close $in0;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
	next if defined($both{$arr[0]});
	print $out1 "$line\n";
}
close $in1;

close $out1;
