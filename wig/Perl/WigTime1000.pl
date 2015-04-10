#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";

my $type = 0;
my $chr = "INIT";
my $lastchr = "INIT";
my $start = 0;
my $count = 0;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /track/;
	if ($line =~ /fixedStep/) {
		$type = 1;
		($chr, $start) = $line =~ /chrom=(.+) start=(\d+) /;
		$count = 0;
		print $out "variableStep chrom=$chr span=1\n" if $chr ne $lastchr;
		$lastchr = $chr;
	}
	elsif ($line =~ /variableStep/) {
		print $out "$line\n";
	}
	else {
		if ($type == 0) {
			my ($pos, $val) = split("\t", $line);
			$val *= 1000;
			print $out "$pos\t$val\n" if $val != 0;
		}
		else {
			my $pos = $start + $count;
			my $val = $line;
			$val *= 1000;
			print $out "$pos\t$val\n" if $val != 0;
			$count ++;
		}
	}
}

close $in;
close $out;
