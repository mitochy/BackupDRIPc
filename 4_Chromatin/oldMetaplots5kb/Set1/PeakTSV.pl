#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
die if $input1 =~ /TSS/;

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
open (my $outTSS, ">", "TSSPeak\_$fileName1.tsv") or die "Cannot write to: $!\n";
open (my $outDRIPc, ">", "DRIPcPeak\_$fileName1.tsv") or die "Cannot write to: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
	#my ($TSSstart, $TSSend) = $info =~ /ORIG=E\w+\.\d+,chr\w+,(\d+),(\d+);/ if $input1 !~ /TSS/;
	my ($TSSstart, $TSSend) = $info =~ /ORIG=E\w+\.\d+,chr\w+,(\d+),(\d+);/ if $input1 =~ /orig/;
	($TSSstart, $TSSend) = $info =~ /TWIN=E\w+\.\d+,chr\w+,(\d+),(\d+)/ if $input1 =~ /shuf/;
	$TSSstart -= 3000;
	$TSSend += 3000;
	my ($DRIPstart, $DRIPend) = $info =~ /DRIP=chr\w+,(\d+),(\d+);/;
	my $length = $DRIPend - $DRIPstart;
	my $newstart = int(($start + $end) / 2) - int($length / 2);
	my $newend = int(($start + $end) / 2) + int($length / 2);
	my $newTSSstart = int(($TSSstart + $TSSend) / 2) - 500;
	my $newTSSend   = int(($TSSstart + $TSSend) / 2) + 500;
	# Now print TSV over TSS
	print $outTSS "$name\t";
	my @temp;
	for (my $i = $TSSstart + 100; $i <= $TSSend - 100; $i += 100) {
		push(@temp,  "0") if $i < $newstart;
		push(@temp,  "1") if $i >= $newstart and $i <= $newend;
		push(@temp,  "0") if $i > $newend;
	}
	@temp = reverse(@temp) if $strand eq "-";
	my $result = join("\t", @temp); @temp = ();
	print $outTSS "$result\n";
	# Now print TSV over dripc centered region
	print $outDRIPc "$name\t";
	for (my $i = $start + 100; $i <= $end - 100; $i += 100) {
		push(@temp, "0") if $i < $newTSSstart;
		push(@temp, "1") if $i >= $newTSSstart and $i <= $newTSSend;
		push(@temp, "0") if $i > $newTSSend;
	}
	@temp = reverse(@temp) if $strand eq "-";
	$result = join("\t", @temp); @temp = ();
	print $outDRIPc "$result\n";
}
close $in1;
