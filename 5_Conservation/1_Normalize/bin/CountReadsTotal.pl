#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");
my @chr = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY);
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
#open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";
my $linecount = 0;
my $total = 0;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /^#/;
	my ($read, $number, $chr, $pos, $test, $mapped) = split("\t", $line);
	next if not grep(/^$chr$/, @chr);
	#print $out "$line\n" if $chr !~ /\*/;
	#while ($mapped =~ /[A-Z]/g) {
	#	my $last = $`;
	#	my $curr = $&;
	#	my $next = $';
	#	my ($num) = $last =~ /[A-Z]*(\d+)$/;
	#	$num = 0 if not defined($num);
	#	$total += $num;
	#}
	my ($num) = $mapped =~ /^(\d+)M$/;
	$total += $num;
	$linecount++;
}

close $in;
#close $out;

print "$input: $total\n";
