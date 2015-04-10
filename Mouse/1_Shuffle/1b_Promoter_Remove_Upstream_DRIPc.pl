#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
my ($pwd) = `pwd`;
die "\nusage: $0 <anything_promoter.bed>\n\nTHIS SCRIPT CAN BE USED FOR TXT (H3K4me3_dripc_promoter.txt) or .bed (dripc_promoter.input1)\n\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
open (my $out1, ">", "$fileName1.promfixed") or die "Cannot write to $fileName1.promfixed: $!\n";
open (my $out2, ">", "$fileName1.removed") or die "Cannot write to $fileName1.removed: $!\n";
print "This is original\n" if $fileName1 =~ /orig/i;
print "This is shuffled\n" if $fileName1 =~ /shuf/i;

while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;

	# Parse
	my @arr = split("\t", $line); my $column = @arr;
	die "Wrong format: Only has $column columns (tab separated) - has to be 7 or 8 columns in format: chr start end [value] name value strand DRIP=chr1,etc\n" if @arr != 7 and @arr != 8;
	my ($chr, $start, $end, $value1, $name, $value, $strand, $info) = split("\t", $line) if @arr == 8;
	   ($chr, $start, $end, $name, $value, $strand, $info) = split("\t", $line) if @arr == 7;

	# Get center of TSS
	my ($orig, $TSSchr, $TSSstart, $TSSend) = $info =~ /ORIG=(ENS\w+\.\d+),(chr\w+),(\d+),(\d+)/ if $fileName1 =~ /orig/i;
	my $shuf;
	   ($shuf, $TSSchr, $TSSstart, $TSSend) = $info =~ /TWIN=(ENS\w+\.\d+),(chr\w+),(\d+),(\d+)/ if $fileName1 =~ /shuf/i;
	die "Undefined TSSchr/start/end (input1 = $input1 pwd = $pwd)\n" if not defined($TSSend);
	#die "Undefined initial Origin or twin at input1 $input1 info: $info\nline:\n\n$line\n\n" if not defined($orig) or not defined($shuf);
	my $TSSCenter = int(($TSSstart + $TSSend) / 2);
	#($start, $end) = $info =~ /DRIP=chr\w+,(\d+),(\d+)/;
	# Get start of peak
	#my $PeakStart = $strand eq "+" ? $start : $end;
	my $PeakStart = ($start + $end) /2;
	
	# Print out peak which start (or end if - strand) is DOWNSTREAM of TSS
	# + strand: PeakStart > TSSCenter ------TSS->>>-DRIPc->>>
	# - strand: PeakStart < TSSCenter <<<-DRIPc-<<<-TSS------
	
	my $lineout = "$chr\t$start\t$end\t$orig\t$value\t$strand\t$info" if @arr == 7 and ($fileName1=~ /orig/i);
	$lineout = "$chr\t$start\t$end\t$value1\t$orig\t$value\t$strand\t$info" if @arr == 8 and ($fileName1=~ /orig/i);
	$lineout = "$chr\t$start\t$end\t$shuf\t$value\t$strand\t$info" if @arr == 7 and ($fileName1=~ /shuf/i);
	$lineout = "$chr\t$start\t$end\t$value1\t$shuf\t$value\t$strand\t$info" if @arr == 8 and ($fileName1=~ /shuf/i);

	print $out1 "$lineout\n" and next if $PeakStart > $TSSCenter and $strand eq "+";
	print $out1 "$lineout\n" and next if $PeakStart < $TSSCenter and $strand eq "-";
	# Get name of ORIG and SHUF
	my ($origName, $shufName) = $info =~ /ORIG=(ENS\w+\.\d+),chr.+TWIN=(ENS\w+\.\d+),chr/; 
	die "Undefined Origin or twin at input1 $input1 line:\n\n$line\n\n" if not defined($origName) or not defined($shufName);
	print $out2 "$lineout\n";
}
close $in1;
close $out1;
close $out2;

print "Output; $fileName1.promfixed\nRemoved peaks: $fileName1.removed\n";
