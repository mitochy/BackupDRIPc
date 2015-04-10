#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1, $input2) = @ARGV;
die "usage: $0 <input> <input2>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
my ($folder2, $fileName2) = mitochy::getFilename($input2, "folder");

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $out2, ">", "$fileName2.out") or die "Cannot write to $fileName2.out: $!\n";

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";

my %wig;
my $chr;
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	if ($line =~ /chrom=/) {
		($chr) = $line =~ /chrom=(.+) span/;
		die "Not defined chrom at lien $line\n" if not defined($chr);
	}
	else {
		my ($pos, $val) = split("\t", $line);
		$wig{$chr}{$pos} = $val;
	}
}
close $in1;

open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
my $lastchr = "INIT";

while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	if ($line =~ /chrom=/) {
		($chr) = $line =~ /chrom=(.+) span/;
		die "Not defined chrom at lien $line\n" if not defined($chr);
		if ($chr ne $lastchr) {
			print $out1 "$line\n";
			print $out2 "$line\n";
		}
		$lastchr = $chr;
	}
	else {
		my ($pos, $val) = split("\t", $line);
		next if not defined($wig{$chr}{$pos});
		print $out1 "$pos\t$wig{$chr}{$pos}\n";
		print $out2 "$pos\t$val\n";
	}
}

close $in2;

