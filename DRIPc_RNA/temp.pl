#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input, $gene) = @ARGV;
die "usage: $0 <input> <gene>\n" unless @ARGV == 2;

my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";
my %data;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
	$data{$arr[4]} = $arr[3];
}

close $in;

open (my $in2, "<", $gene) or die;
while (my $line = <$in2>) {
	chomp($line);
	my ($chr, $start, $end, $name, $junk, $strand) = split("\t", $line);
	if (defined($data{$name}) and $data{$name} =~ /\d+/) {
		print $out "$chr\t$start\t$end\t$name\t$data{$name}\t$strand\n";
	}
}
close $in2;
close $out;
