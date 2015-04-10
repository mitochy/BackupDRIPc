#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;

die "

This script parse RNAseq data (NT2.rpkm) and dripc_promoter/genebody/terminal.input file and put RNA into value

usage: $0 <input1>

" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

# Parse RNA
my %rna;
open (my $in1, "<", "../../data/NT2.rpkm") or die "Cannot read from ../../data/NT2.rpkm: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($gene, $rna) = split("\t", $line);
	$rna{$gene} = $rna;
}
close $in1;

my %rna2;
open (my $in2, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $value, $strand, $info) = split("\t", $line);
	my $rna;
	my $namez = "$chr\_$start\_$end\_$strand";
	my $gene;
	if ($name =~ /^ENS/) {
		$gene = $name;
		$rna2{$namez} = defined($rna{$gene}) ? $rna{$gene} : 0;
	}
	else {
		($gene) = $info =~ /TWIN=(ENST\w+\.\d+),chr/;
		$rna2{$namez} = defined($rna{$gene}) ? $rna{$gene} : 0;
	}
}
close $in2;

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $in3, "<", "dripc_main.input") or die "Cannot read from dripc_main.input: $!\n";
while (my $line = <$in3>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $value, $strand) = split("\t", $line);
	my $namez = "$chr\_$start\_$end\_$strand";
	my $rna = defined($rna2{$namez}) ? $rna2{$namez} : 0;
	print $out1 "$chr\t$start\t$end\t$rna\t$name\t$value\t$strand";
}
close $in3;
close $out1;
