#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $in, "<", "../../bed/hg19_gencode.bed") or die "Cannot read from ../../bed/hg19_gencode.bed: $!\n";
my %gene;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	$gene{$name}{start} = $start;
	$gene{$name}{end} = $end;
	$gene{$name}{strand} = $strand;
}

close $in;
open (my $in2, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";

while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $zero, $strand, $info) = split("\t", $line);
	my ($dripOrig, $geneOrig, $twinOrig) = split(";", $info);
	my ($chrdrip, $startdrip, $enddrip) = split(",", $dripOrig);
	my ($nameOrig, $chrOrig, $startOrig, $endOrig) = split(",", $geneOrig);
	$nameOrig =~ s/ORIG=//;
	my ($nameTwin, $chrTwin) = split(",", $twinOrig);
	$nameTwin =~ s/TWIN=//;
	my ($startpos, $endpos) = ($startdrip - $startOrig, $endOrig - $enddrip);
	die "Undefined $nameOrig\n" unless defined($gene{$nameOrig}{strand});
	die "Undefined $nameTwin\n" unless defined($gene{$nameTwin}{strand});
	my $strandOrig = $gene{$nameOrig}{strand};
	my $strandTwin = $gene{$nameTwin}{strand};
	my ($startTwin, $endTwin) = $strandTwin eq "+" ? ($gene{$name}{start}-2000,$gene{$name}{start}+2000) : ($gene{$name}{end}-2000,$gene{$name}{end}+2000);
	($start, $end) = $strandTwin eq $strandOrig ? ($startTwin + $startpos, $endTwin - $endpos) : ($startTwin + $endpos, $endTwin - $startpos);
	$twinOrig = "$name,$chr,$startTwin,$endTwin";
	my $newStrand = $strandTwin eq "+" ? "-" : "+";
	my ($newstart, $newend);
	if ($input =~ /CENTER/) {
		$newstart = $start + int(($start + $end) / 2) - 5000;
		$newend   = $start + int(($start + $end) / 2) + 5000;
	}
	else {
		$newstart = $startTwin;
		$newend   = $endTwin;
	}
	print $out "$chr\t$newstart\t$newend\t$name\t$zero\t$newStrand\t$dripOrig;$geneOrig;$twinOrig\n";
}

close $in2;

close $out;
