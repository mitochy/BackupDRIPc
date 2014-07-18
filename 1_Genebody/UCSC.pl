#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");

http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1%3A101427998-101471719&hgsid=382626139_wUtK1DqTeJcoTym83mu8Tz2lOuHM

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";

while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}

close $in;
close $out;
