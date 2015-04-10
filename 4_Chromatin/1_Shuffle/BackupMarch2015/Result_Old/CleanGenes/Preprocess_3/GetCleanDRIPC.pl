#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <dripcNODRIP_promoter.txt>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my $desiredPeak = 20;
# Get list of clean genes
my %clean;
my $input2 = "gene_clean.out";
open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name) = split("\t", $line);
	$clean{$name} = 1;
}
close $in2;

my %orig;
open (my $in3, "<", "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_promoter.bed") or die "Cannot read from /data/mitochi/Work/Project/DRIPc/1_Location/dripc_promoter.bed: $!\n";
while (my $line = <$in3>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand) = split("\t", $line);
	$orig{$name} = "$chr\t$start\t$end\t$name\t$val\t$strand";
}
close $in3;

# Print shuffles that only have clean original and twin
my %data;
my $totalPeak = 0;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $val1, $name, $val2, $strand, $info) = split("\t", $line);
	my ($twin) = $line =~ /TWIN=(\w+\.\d+),chr/;
	my ($orig) = $line =~ /ORIG=(\w+\.\d+),chr/;
	my ($drip) = $line =~ /DRIP=(.+);ORIG=/;
	$totalPeak ++;
	next if not defined($clean{$twin}) or not defined($clean{$orig});
	push(@{$data{$name}{shuf}}, $line);
	$data{$name}{orig} = "$orig{$name}\t$info";
}
close $in1;


open (my $out1, ">", "$fileName1\_orig.bed") or die "Cannot write to $fileName1\_orig.bed: $!\n";
open (my $out2, ">", "$fileName1\_shuf.bed") or die "Cannot write to $fileName1\_shuf.bed: $!\n";
my $totalPeakClean 	= 0;
my $used		= 0;
foreach my $name (keys %data) {
	$totalPeakClean ++;
	next if (@{$data{$name}{shuf}} < $desiredPeak);
	$used ++;
	print $out1 "$data{$name}{orig}\n";
	foreach my $shuf (@{$data{$name}{shuf}}) {
		print $out2 "$shuf\n";
	}
}
close $out1;
close $out2;

print "Total Peak: $totalPeak, Total Clean Peak: $totalPeakClean, Used: $used\n";
