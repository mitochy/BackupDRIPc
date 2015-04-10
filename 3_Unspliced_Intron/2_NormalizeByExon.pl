#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input, $input2) = @ARGV;
die "usage: $0 <gtf> <rnaseq count>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input2, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";
print "Output = $name.out\n";
my %gene;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $source, $type, $start, $end, $zero, $strand, $dot, $name) = split("\t", $line);
	($name) = $name =~ /gene_id \"(.+)\"; transcrip/;
	my $length = $end - $start;
	$gene{$name} = $length;
}

close $in;

open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";

while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name, $count) = split("\t", $line);
	my $length = $gene{$name};
	die "Died at $name\n" unless defined($length);
	$count = $count * 100 / $length; # times 100 coz each read is 100
	print $out "$name\t$count\n";
}

close $in2;
close $out;
