#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$fileName.dist") or die "Cannot write to $fileName.dist: $!\n";

while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
	my ($origStart, $origEnd) = $info =~ /ORIG=\w+\.\d+,chr\w+,(\d+),(\d+);/;
	my $centerOrig = ($origEnd - $origStart) / 2 + $origStart;
	my $centerDRIP = ($end - $start) / 2 + $start;
	my $dist =  $centerOrig - $centerDRIP;
	$dist = -1 * $dist if ($strand eq "-");
	#print $out "$dist\tORIG: ($origEnd - $origStart) / 2 + $origStart = $centerOrig, DRIP ($end - $start) / 2 + $start = $centerDRIP\n";
	print $out "$dist\n";
}

close $in;
close $out;
print "Output: $fileName.dist\n";

__END__


10         5010      10010
A----------C----------B
   D-----T-----E
       4020


dist = 4020 - 5010
dist = (B - A)/2 + A - (E - D) / 2 + D
