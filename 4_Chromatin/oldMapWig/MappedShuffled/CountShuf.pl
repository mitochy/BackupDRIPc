#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my %data;
open (my $out1, ">", "$fileName1.fixed") or die "Cannot write to $fileName1.fixed: $!\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
        chomp($line);
        next if $line =~ /#/;
        my @arr = split("\t", $line);
        die "Wrong format col has to be 7 or 8\n" if @arr != 7 and @arr != 8;
        my ($chr, $start, $end, $value1, $gene, $value, $strand, $info) = split("\t", $line) if @arr == 8;
           ($chr, $start, $end, $gene, $value, $strand, $info) = split("\t", $line) if @arr == 7;
        my ($orig, $shuf) = $info =~ /ORIG=(ENS\w+\.\d+),chr\w+,\d+,\d+;TWIN=(ENS\w+\.\d+),chr\w+,\d+,\d+/;
        die "Undefined orig shuf at lien $line\n" if not defined($orig) or not defined($shuf);
	$data{$orig} ++;
        #my $name = "$orig\_$shuf";
        #next if (defined($set{$name}));
        #print $out "$line\n";

}
close $in1;

#close $out1;

my %count;
foreach my $orig (sort {$data{$b} <=> $data{$a}} keys %data) {
	$count{$data{$orig}} ++;
}

foreach my $count (sort {$a <=> $b} keys %count) {
	last if $count == 125;
	print "$count: $count{$count}\n";
}
