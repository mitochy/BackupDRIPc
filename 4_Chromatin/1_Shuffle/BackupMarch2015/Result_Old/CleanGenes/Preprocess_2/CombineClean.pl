#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 outputname\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

my %coor;
open (my $in, "<", "hg19_gencode.bed") or die "Cannot read from hg19_gencode.bed: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	$coor{$name} = $line;
}
close $in;

my %clean;
open (my $in1, "<", "promoter_clean.bed") or die "Cannot read from promoter_clean.bed: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $names, $strand) = split("\t", $line);
	my @names = split(";", $names);
	foreach my $name (@names) {
		$clean{$name} ++;
	}
	
}
close $in1;

open (my $in2, "<", "terminal_clean.bed") or die "Cannot read from terminal_clean.bed: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $names, $strand) = split("\t", $line);
	my @names = split(";", $names);
	foreach my $name (@names) {
		$clean{$name} ++;
	}
	
}
close $in2;

open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";
foreach my $name (keys %clean) {
	print $out "$coor{$name}\n" if $clean{$name} == 2;
}
close $out;

system("sort -k1,1 -k2,2n $fileName.out > $fileName.TEMP && mv $fileName.TEMP $fileName.out");
