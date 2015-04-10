#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($bedFile, $outputFile) = @ARGV;
die "usage: $0 <Bed> <output File>\n" unless @ARGV == 2;

my ($folder, $fileName) = mitochy::getFilename($bedFile, "folder");
die "$outputFile exist: will not overwrite\n" if -e $outputFile;
my %data;
open (my $in1, "<", $bedFile) or die "Cannot read from $bedFile: $!\n";
open (my $out, ">", "$outputFile") or die "Cannot write to $outputFile: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand) = split("\t", $line);
	$strand = 0 if $strand eq "+";
	$strand = 16 if $strand eq "-";
	my $pos = $end - $start;
	$pos .= "M";
	print $out "$name\t$strand\t$chr\t$start\t255\t$pos\t\*\t0\t0\tNA\tNA\n";
}
close $in1;
close $out;

__END__
