#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <dripc_promoter.shuffled or genebody_shuffled.bed>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

my %data;
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
	my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line) if @arr <= 7;
	($chr, $start, $end, $val, $name, $val, $strand, $info) = split("\t", $line) if @arr == 8;
	my ($orig, $shuf) = $line =~ /ORIG=(\w+\.\d+),chr.+TWIN=(\w+\.\d+),chr/;
	$data{$orig} ++;
}
close $in;

my %shuffled;
foreach my $name (keys %data) {
	$shuffled{$data{$name}} ++;
}

foreach my $count (sort {$b <=> $a} keys %shuffled) {
	print "$count\t$shuffled{$count}\n";
}
