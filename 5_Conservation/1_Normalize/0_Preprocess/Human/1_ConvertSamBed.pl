#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($bedFile) = @ARGV;
die "usage: $0 <Bed>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($bedFile, "folder");

my %data;
open (my $in1, "<", $bedFile) or die "Cannot read from $bedFile: $!\n";
open (my $out, ">", "$fileName.sam.out") or die "Cannot write to $fileName.sam.out: $!\n";
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
open (my $in2, "<", $samFile) or die "Cannot read from $samFile: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}
close $in2;


