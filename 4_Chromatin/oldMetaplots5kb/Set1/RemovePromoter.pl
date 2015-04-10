#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my $removeFile = $input1 =~ /orig/i ? "../../dripc_promoter_orig.removed" : $input1 =~ /shuf/i ? "../../dripc_promoter_shuf.removed" : die "Undefined file type\n";
print "Remvfile = $removeFile\n";
my %set;
open (my $in, "<", $removeFile) or die;
while (my $line = <$in>) {
	chomp($line);
	my ($id, $orig, $shuf) = split("\t", $line);
	my $name = "$orig\_$shuf";
	$set{$name} = 1;
}
close $in;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
open (my $out, ">", "$fileName1.promfixed") or die;

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
	die "Wrong format col has to be 7 or 8\n" if @arr != 7 and @arr != 8;
	my ($chr, $start, $end, $value1, $gene, $value, $strand, $info) = split("\t", $line);
	my ($orig, $shuf) = $info =~ /ORIG=(ENS\w+\.\d+),chr\w+,\d+,\d+;TWIN=(ENS\w+\.\d+),chr\w+,\d+,\d+/;
	die "Undefined orig shuf at lien $line\n" if not defined($orig) or not defined($shuf);
	my $name = "$orig\_$shuf";
	next if (defined($set{$name}));
	print $out "$line\n";

}
close $in1;

