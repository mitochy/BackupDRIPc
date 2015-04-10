#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";

my $count = 0;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	$count++;
	my ($chr, $start, $end, $name) = split("\t", $line);
	$line =~ s/$name/$name$count/;
	print $out "$line\n";
}

close $in;
close $out;
