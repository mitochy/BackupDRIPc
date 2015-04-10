#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";

my %bed;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand) = split("\t", $line);
	$bed{$chr}{$name} = $line;
}

close $in;

my %bad;
my %good;
foreach my $chr (sort keys %bed) {
	print "Processing chr $chr\n";
	foreach my $name (keys %{$bed{$chr}}) {
		foreach my $name2 (keys %{$bed{$chr}}) {
			next if $name eq $name2;
			my ($chr1, $start1, $end1, $name1, $val1, $strand1) = split("\t", $bed{$chr}{$name});
			my ($chr2, $start2, $end2, $name2, $val2, $strand2) = split("\t", $bed{$chr}{$name2});
			next if $end2 < $start1 - 5000;
			next if $start2 > $end1 + 5000;
			next if ($strand1 eq $strand2);
			if ($strand1 eq "+") {
				if ($end2 >= $start1 - 3000 and $end2 <= $start1 + 3000) {
					$bad{$name2} = 1;
				}
			}
			if ($strand1 eq "-") {
				if ($start2 >= $end1 - 3000 and $start2 <= $end1 + 3000) {
					$bad{$name2} = 1;
				}
			}
		}
	}
}

open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";
foreach my $name (keys %bad) {
	print $out "$name\n";
}
close $out;
