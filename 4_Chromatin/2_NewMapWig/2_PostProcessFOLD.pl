#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my %use;
open (my $in0, "<", "USE") or die "Cannot read from USE: $!\n";
while (my $line = <$in0>) {
	chomp($line);
	next if $line =~ /#/;
	$use{$line} = 1;
}
close $in0;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($feat, $chromatin, @arr) = split("\t", $line);
	next if not defined($use{$chromatin});
	next if $use{$chromatin} !~ /^1$/;
	print $out1 "$line\n";
}
close $in1;
close $out1;
