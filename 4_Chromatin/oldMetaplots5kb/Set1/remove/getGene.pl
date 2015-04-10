#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std;
use vars qw($opt_o);
getopts("o:");

my ($geneList, $input, $column) = @ARGV;
die "usage: $0 <bed genelist> <bed all gene> <column (def: 0 or 3)> -o output\n" unless @ARGV == 3;

my ($folder, $name) = mitochy::getFilename($geneList, "folder");

open (my $in, "<", $geneList) or die "Cannot read from $geneList: $!\n";
my %data;
my $linecount = 0;
while (my $line = <$in>) {
	chomp($line);
	$linecount++;
	next if $line =~ /^#/;
	my @arr = split("\t", $line);
	if (defined($column) and $column =~ /^\d+/) {
		$data{$arr[$column]}{line} = $linecount;
		$data{$arr[$column]}{count} ++;
	}
	else { 
		if (@arr > 1) {
			$data{$arr[3]}{line} = $linecount;
			$data{$arr[3]}{count} ++;
		}
		else {
			$data{$arr[0]}{line} = $linecount;
			$data{$arr[0]}{count} ++;
		}
	}
}

close $in;
my $output = defined($opt_o) ? $opt_o : "$name.out";
open (my $in2, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$output") or die "Cannot write to $output: $!\n";
print "Output = $output\n";
my %bed;
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /^#/;
	my @arr = split("\t", $line);
	if (defined($data{$arr[3]}{line}) and $data{$arr[3]}{line} =~ /\d+/) {
		my $count = $data{$arr[3]}{count};
		$bed{$arr[3]}{line}{$count} = $line;
		$bed{$arr[3]}{cnt}{$count} = $data{$arr[3]};
		$data{$arr[3]}{count} --;
	}
}

close $in2;

foreach my $gene (sort keys %bed) {
	foreach my $count (keys %{$bed{$gene}{line}}) {
		print $out "$bed{$gene}{line}{$count}\n";
	}
}

close $out;
