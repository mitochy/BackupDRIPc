#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1, $input2) = @ARGV;
die "usage: $0 <orig> <shuf>\n" unless @ARGV == 2;

my $linecount = 0;

open (my $in0, "<", "../../../../data/NT2.rpkm") or die "Cannot read from NT2.rpkm: $!\n";
my %rna;
while (my $line = <$in0>) {
	chomp($line);
	my ($name, $val) = split("\t", $line);
	$rna{$name} = $val;
}
close $in0;


my ($folder2, $fileName2) = mitochy::getFilename($input2, "folder");
my %shuf;
open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
$linecount = 0;
while (my $line = <$in2>) {
	chomp($line);
	$linecount ++;
	my ($shufgene, @arr) = split("\t", $line);
	$shufgene =~ s/>//;
	my $rna = defined($rna{$shufgene}) ? $rna{$shufgene} : 0;
	my $index = getIndex($rna);
	$shuf{$index}{$linecount}{name} = $shufgene;
	$shuf{$index}{$linecount}{rna} = $rna;
	$shuf{$index}{$linecount}{dens}[1] = mean(@arr[18..27]);
	$shuf{$index}{$linecount}{dens}[2] = mean(@arr[28..37]);
	$shuf{$index}{$linecount}{dens}[3] = mean(@arr[38..47]);
	$shuf{$index}{$linecount}{dens}[4] = mean(@arr[48..57]);
}
close $in2;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
$linecount = 0;
my %orig;
while (my $line = <$in1>) {
	chomp($line);
	$linecount ++;
	my ($origgene, @arr) = split("\t", $line);
	my $rna = defined($rna{$origgene}) ? $rna{$origgene} : 0;
	my $index = getIndex($rna);
	my ($indexMin, $indexMax) = (0,0);
        $indexMin = 10 + int($rna*0.8 / 10)  if $rna >= 10 and $rna <= 100;
        $indexMin = 20 + int($rna*0.8 / 100) if $rna >= 100 and $rna < 1000;
        $indexMin = 30 + int($rna*0.8 / 1000) if $rna >= 1000 and $rna < 10000;
        $indexMin = 40 + int($rna*0.8 / 10000) if $rna >= 10000;
        $indexMax = 10 + int($rna*1.2 / 10)  if $rna >= 10 and $rna <= 100;
        $indexMax = 20 + int($rna*1.2 / 100) if $rna >= 100 and $rna < 1000;
        $indexMax = 30 + int($rna*1.2 / 1000) if $rna >= 1000 and $rna < 10000;
        $indexMax = 40 + int($rna*1.2 / 10000) if $rna >= 10000;
	my @dens;
	#$orig{$linecount}{name} = $origgene;
	#$orig{$linecount}{rna} = $rna;
	$dens[1] = mean(@arr[18..27]);
	$dens[2] = mean(@arr[28..37]);
	$dens[3] = mean(@arr[38..47]);
	$dens[4] = mean(@arr[48..57]);

	for (my $i = $indexMin; $i <= $indexMax; $i++) {
		next if not defined($shuf{$i}) or (keys $shuf{$i}) == 0;
		foreach my $shuflinecount (sort {$shuf{$i}{$a}{rna} <=> $shuf{$i}{$b}{rna}} keys %{$shuf{$i}}) {
        		next if $shuf{$i}{$shuflinecount} <= $rna * 0.8;
        		last if $shuf{$i}{$shuflinecount} >= $rna * 1.2;
			my $check = 0;
			for (my $j = 1; $j < 5; $j++) {
				$check = 1 if $shuf{$linecount}{dens}[$i] < $dens[$i] - 0.1 or $shuf{$linecount}{dens}[$i] > $dens[$i] + 0.1;
			}
			print $out1 "$linecount\_$shuflinecount\n" if $check != 1;
		}
        }
}
close $in1;
close $out1;

sub mean {
	my @arr = @_;
	my $total;
	my $count;
	foreach my $arr (@arr) {
		next if $arr =~ /NA/i;
		$total += $arr;
		$count ++;
	}
	return($total / $count);
}

sub getIndex {
	my ($val) = @_;
	my $index = 0;
        $index = 10 + int($val / 10)  if $val >= 10 and $val <= 100;
        $index = 20 + int($val / 100) if $val >= 100 and $val < 1000;
        $index = 30 + int($val / 1000) if $val >= 1000 and $val < 10000;
        $index = 40 + int($val / 10000) if $val >= 10000;
	return($index);
}


