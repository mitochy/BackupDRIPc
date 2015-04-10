#!/usr/bin/perl
# This script averages RNAseq rpkm data


use strict; use warnings; use mitochy;

my ($input1, $input2) = @ARGV;
die "usage: $0 <.rpkm rep1> <.rpkm rep2>\n" unless @ARGV == 2;

my ($folder1, $name1) = mitochy::getFilename($input1, "folder");
my ($folder2, $name2) = mitochy::getFilename($input2, "folder");

# Process count file to get count and rpkm
my %data;
print "1. Processing count file $input1\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($genes, $count) = split("\t", $line);
	my @genes = split(";", $genes);
	foreach my $gene (@genes) {
		$data{$gene} += $count/2;
	}	
}
close $in1;
print "2. Processing count file $input2\n";
open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($genes, $count) = split("\t", $line);
	my @genes = split(";", $genes);
	foreach my $gene (@genes) {
		$data{$gene} += $count/2;
	}	
}
close $in2;

open (my $out, ">", "NT2.rpkm") or die "Cannot write to NT2.rpkm: $!\n";
print "3. Printing RPKM\n";
foreach my $gene (sort keys %data) {
	my $rpkm = $data{$gene};
	print $out "$gene\t$rpkm\n";
}
close $out;

print "Output: NT2.rpkm\n";
