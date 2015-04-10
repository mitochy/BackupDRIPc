#!/usr/bin/perl
# This script find gene transcript that has similar (+/- 20%) RPKM as the other
# Output: 
# 1. 3T3_RNAtwin.rpkm 
# 2. 3T3_RNAzero.rpkm for all genes with less than 10 RPKM

use strict; use warnings; use mitochy;

my ($input) = "/data/mitochi/Work/Project/DRIPc/data/3T3.rpkm";
print "1. Processing $input\n";
my %data;
my %res;
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($gene, $val)  = split("\t", $line);
	my $index;
	$index = 0 if ($val < 10);
	$index = 10 + int($val / 10)  if $val >= 10 and $val <= 100;
	$index = 20 + int($val / 100) if $val >= 100 and $val < 1000;
	$index = 30 + int($val / 1000) if $val >= 1000 and $val < 10000;
	$index = 40 + int($val / 10000) if $val >= 10000;
	$data{$index}{$gene} = $val;
}
close $in;

print "2. Finding similar exp genes\n";
mkdir "../1_Shuffle/" if not -d "../1_Shuffle/";
open (my $out, ">", "../1_Shuffle/3T3_RNAtwin.rpkm") or die "Cannot write to 3T3_RNAtwin.rpkm: $!\n";
open (my $out2, ">", "../1_Shuffle/3T3_RNAzero.rpkm") or die "Cannot write to 3T3_RNAtzero.rpkm: $!\n";
my $count = 0;
my ($total) = `wc -l $input` =~ /^(\d+)/;
my @zero;
foreach my $index (sort {$a <=> $b} keys %data) {
	print "Processing INDEX $index\n";
	foreach my $gene (sort {$data{$index}{$a} <=> $data{$index}{$b}} keys %{$data{$index}}) {
		#print "Processing $gene\n";
		my $val = $data{$index}{$gene};
		printf "Processed $count / $total (%.2f)\n", $count / $total * 100 if $count % 200 == 0;
		$count++;
		if ($index == 0) {
			print $out2 "$gene\n";
		}
		else {
			my ($indexMin, $indexMax);
			$indexMin = 10 + int($val*0.8 / 10)  if $val >= 10 and $val <= 100;
			$indexMin = 20 + int($val*0.8 / 100) if $val >= 100 and $val < 1000;
			$indexMin = 30 + int($val*0.8 / 1000) if $val >= 1000 and $val < 10000;
			$indexMin = 40 + int($val*0.8 / 10000) if $val >= 10000;
			$indexMax = 10 + int($val*1.2 / 10)  if $val >= 10 and $val <= 100;
			$indexMax = 20 + int($val*1.2 / 100) if $val >= 100 and $val < 1000;
			$indexMax = 30 + int($val*1.2 / 1000) if $val >= 1000 and $val < 10000;
			$indexMax = 40 + int($val*1.2 / 10000) if $val >= 10000;
			print $out "$gene";
			for (my $i = $indexMin; $i <= $indexMax; $i++) {
				next if not defined($data{$i}) or (keys $data{$i}) == 0;
				foreach my $gene2 (sort {$data{$i}{$a} <=> $data{$i}{$b}} keys %{$data{$i}}) {
					next if $gene eq $gene2;
					next if $data{$i}{$gene2} <= $val * 0.8;
					last if $data{$i}{$gene2} >= $val * 1.2;
					print $out "\t$gene2";
				}
			}
			print $out "\n";
		}
	}
}
