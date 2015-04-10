#!/usr/bin/perl
#
# Total Reads:
# NT2_hg192.sam.out       24697838
# Fibro_hg192.sam.out     10168581
# E14_hg19.sam.out         3056096
# 3T3_hg19.sam.out        23118079
# NT2_mm10.sam.out        24777637
# Fibro_mm10.sam.out      10207869
# E14_mm102.sam.out        3046511
# 3T3_mm102.sam.out       23047346
#
# Total Basepair Reads:
# NT2_hg192.sam.out:
# Fibro_hg192.sam.out:     848794144
# E14_hg19.sam.out:        345987161
# 3T3_hg19.sam.out:       2627116591
# NT2_mm10.sam.out:       1269862239
# Fibro_mm10.sam.out:      949043234
# E14_mm102.sam.out:       279130387
# 3T3_mm102.sam.out:      2118141611
# hg19 normalize to E14 (345987161)
# mm10 normalize to E14 (279130387)

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input.sam>\n" unless @ARGV;

# EDIT: For some reason Fibro and 3T3 are not normalized even though their reads are. WIG file shows 3T3 are 0.5x NT2 and E14,
# while Fibro 1.3. Therefore these are going to be further normalized by 2x (3T3) or 0.75x (Fibro)

my $totalBasepair = 0;
$totalBasepair = (345987161/1175665022)*24697838 if $input =~ /NT2_hg192/;
$totalBasepair = (345987161/ 848794144)*10168581*0.75 if $input =~ /Fibro_hg192/;
$totalBasepair = (345987161/ 345987161)* 3056096 if $input =~ /E14_hg19/;
$totalBasepair = (345987161/2627116591)*23118079*2.05 if $input =~ /3T3_hg19/;
$totalBasepair = (279130387/1269862239)*24777637 if $input =~ /NT2_mm10/;
$totalBasepair = (279130387/ 949043234)*10207869*0.75 if $input =~ /Fibro_mm10/;
$totalBasepair = (279130387/ 279130387)* 3046511 if $input =~ /E14_mm102/;
$totalBasepair = (279130387/2118141611)*23047346*2.05 if $input =~ /3T3_mm102/;

my $totalRead = 0;
$totalRead = 24697838 if $input =~ /NT2_hg192/;
$totalRead = 10168581 if $input =~ /Fibro_hg192/;
$totalRead =  3056096 if $input =~ /E14_hg19/;
$totalRead = 23118079 if $input =~ /3T3_hg19/;
$totalRead = 24777637 if $input =~ /NT2_mm10/;
$totalRead = 10207869 if $input =~ /Fibro_mm10/;
$totalRead =  3046511 if $input =~ /E14_mm102/;
$totalRead = 23047346 if $input =~ /3T3_mm102/;

my %read;
while ((keys %read) < $totalBasepair) {
	my $rand = int(rand($totalRead))+1;
	$rand = $totalRead if $rand > $totalRead;
	$read{$rand} = 1;
}
print "$input: Randomizing done\n";
my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.norm") or die "Cannot write to $name.norm: $!\n";

my $linecount = 0;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /^#/;
	next if $line =~ /chrM/;
	$linecount ++;
	print $out "$line\n" if defined($read{$linecount} and $read{$linecount} == 1);
}

close $in;
close $out;
