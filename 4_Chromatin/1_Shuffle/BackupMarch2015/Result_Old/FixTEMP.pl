#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1, $input2) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input2, "folder");

open (my $in, "<", $input1) or die "Cannot read from $input1: $!\n";
my %data;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $dripname, $val2, $strand, $info) = split("\t", $line);
	my ($orig) = $info =~ /ORIG=(\w+.\d+),chr/;
	die "Died at $input1 line $line\n" unless defined($orig);
	$data{$dripname} = $orig;
	
}

close $in;

open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";

while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $val1, $dripname, $val2, $strand, $info) = split("\t", $line);
	my $name2 = $data{$dripname};
	die "$input2: Died undefined name2 (dripname $dripname) at line $line\n" unless defined($name2);
	$info =~ s/ORIG=\w+,chr/ORIG=$name2,chr/;#die "NEW INFO $name2 $info\n";
	print $out "$chr\t$start\t$end\t$val1\t$dripname\t$val2\t$strand\t$info\n";
	
}

close $in2;
close $out;
