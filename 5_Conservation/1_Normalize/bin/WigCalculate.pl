#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";

my %data;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	next if $line =~ /chrom/;
	next if $line =~ /track/;
	my ($pos, $val) = split("\t", $line);
	$data{categorize($val)} ++;
}

close $in;

foreach my $cat (sort keys %data) {
	print "$input $cat $data{$cat}\n";
}
sub categorize {
	my ($input) = @_;
	return("one") if $input == 1;
	return("five") if $input <= 5;
	return("ten") if $input <= 10;
	return("big");
}
