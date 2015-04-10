#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input, $column) = @ARGV;
die "usage: $0 <input> <column [start from 0]>\n" unless @ARGV == 2;

my %col;
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
	my $col = $arr[$column];
	die "Undefined column value at line: $line\n" unless defined($col);
	$col{$col}++;
}
close $in;

foreach my $col (sort {$col{$b} <=> $col{$a}} keys %col) {
	print "$col\t$col{$col}\n";
}
