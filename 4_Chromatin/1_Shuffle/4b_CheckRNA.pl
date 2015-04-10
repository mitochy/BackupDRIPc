#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input.shuffled>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $in, "<", "/data/mitochi/Work/Project/DRIPc/data/NT2.rpkm") or die "Cannot read from ../../data/NT2.rpkm: $!\n";
my %rna;
while (my $line = <$in>) {
        chomp($line);
        next if $line =~ /#/;
        my ($gene, $rna) = split("\t", $line);
        $rna{$gene} = $rna;
}
close $in;

open (my $in2, "<", $input) or die "Cannot read from $input: $!\n";
my $diff;
my $total;
my @diff;
my $countRNAHigher = 0;
my $countRNALower = 0;
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
	my ($gene, $randGene) = $info =~ /ORIG=(\w+\.\d+),chr.+TWIN=(\w+\.\d+),chr/;
	my $rna1 = $rna{$gene};
	my $rna2 = $rna{$randGene};
	next if ($rna1 < 10);
	push(@diff, ($rna2)/($rna1));
	$diff += ($rna2)/($rna1);
	if ($rna2 / $rna1 < 0.9 or $rna2 / $rna1 > 1.1) {
		my $dif = int($rna1 / $rna2*100) / 100;
		$countRNAHigher ++ if $rna2 / $rna1 > 1.1;
		$countRNALower ++ if $rna2 / $rna1 < 0.9;
		print "$dif Gene $gene randgene $randGene has rna1 $rna1 rna2 $rna2\n" if $countRNAHigher + $countRNALower > 10;
		print "More than 10 RNA shuffles higher than 1.1x or lower than 0.9x of original\n" if $countRNAHigher + $countRNALower > 10;
	}
	$total++;
}
close $in2;
@diff = sort {$a <=> $b} @diff;
my $median = $diff[int(@diff/2)];
$diff /= $total;
open (my $out, ">", "TEMP.dat") or die;
print $out "NAME";
for (my $i = 0; $i < @diff; $i += 100) {
	print $out "\t$diff[$i]";
}
close $out;
print "Diff = $diff, median = $median\n";

print "RNA of shuffles higher than 1.1x of original: $countRNAHigher\n";
print "RNA of shuffles lower than 0.9x of original: $countRNALower\n";
