#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input, $outputFolder) = @ARGV;
die "usage: $0 <input.BED> <output Folder>\n" unless @ARGV == 2;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

my @check_if_prev_result_exist = <$outputFolder\/$fileName\_*.out>;
if (@check_if_prev_result_exist > 0) {
	print "WARNING: Remove all $outputFolder\/$fileName\_*.out? (Y to continue) ";
	my $respond = <STDIN>;
	print "\n";
	if ($respond =~ /y/i) {
		system("rm $outputFolder\/$fileName\_*.out");
	}
	else {
		die;
	}
}

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	open (my $out, ">>", "$outputFolder\/$fileName\_$chr.out") or die;
	print $out "$line\n";
	close $out;

}
close $in;
