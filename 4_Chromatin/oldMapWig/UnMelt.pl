#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
my %data;
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($feature, $type, $junk1, $junk2, $fold) = split("\t", $line);
	$data{$feature}{$type} = $fold;
}
close $in1;

open (my $out1, ">", "$fileName1.fixed") or die "Cannot write to $fileName1.fixed: $!\n";
print $out1 "id\tpromoter\tgenebody\tterminal\n";
foreach my $feature (sort keys %data) {
	foreach my $type (sort keys %{$data{$feature}}) {
		next if not defined($data{promoter}{$type}) or not defined($data{terminal}{$type}) or not defined($data{genebody}{$type});
		print $out1 "$type\t$data{promoter}{$type}\t$data{genebody}{$type}\t$data{terminal}{$type}\n";
	}
	last;
}
close $out1;
