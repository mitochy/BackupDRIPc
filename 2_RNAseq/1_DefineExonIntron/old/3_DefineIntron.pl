#!/usr/bin/perl
# This file convert bed to gtf

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <output of 1_MergeGene.pl>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name\_intron.bed") or die "Cannot write to $name\_intron.bed: $!\n";

my %group;
my $curr_name = "INIT";
my $curr_start;
my $curr_end;
my $curr_intron_number;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	if ($curr_name ne $name) {
		$curr_name = $name;
		$curr_start = $start;
		$curr_end   = $end;
		$curr_intron_number = 1;
	}
	else {
		my $start2 = $curr_end + 1;
		my $end2   = $start - 1;
		my $name2  = "$name\_$curr_intron_number";
		print $out "$chr\t$start2\t$end2\t$name2\t$zero\t$strand\n" if $start2 < $end2;
		$curr_name = $name;
		$curr_start = $start;
		$curr_end   = $end;
		$curr_intron_number ++;
	}
}

close $in;
close $out;
