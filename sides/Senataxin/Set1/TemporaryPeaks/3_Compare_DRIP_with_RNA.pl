#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <SIGDOWN.BED or SIGUP.BED>\n" unless @ARGV;

# Process RNA
my $RNA = "/data/mitochi/Work/Project/DRIPc/sides/Senataxin/Set1/RNA.rpkm";
my %rna;
open (my $in1, "<", $RNA) or die "Cannot read from $RNA: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name, $gene, $val, $qval) = split("\t", $line);
	$rna{$gene} = $val;
}
close $in1;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
open (my $in2, "<", $input1) or die "Cannot read from $input1: $!\n";
open (my $out1, ">", "$fileName1.DRIP_with_RNA") or die "Cannot write to $fileName1.DRIP_with_RNA: $!\n";
print $out1 "Gene Chr\tGene Start\tGene End\tGene\tmRNA Fold Change\tDRIP Log2 Fold Change\tRank\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $names, $fold, $strand, $C, $D, $E, $F, $meanCD, $meanEF, $pval, $chr2, $start2, $end2) = split("\t", $line);
	my @names = split(";", $names);
	my $bestGene;
	my $bestRNA = 0;
	foreach my $name (@names) {
		if (defined($rna{$name})) {
			if (not defined($bestGene)) {
				$bestGene = $name;
				$bestRNA = $rna{$name};
			}
			elsif ($bestRNA > $rna{$name}) {
				die "WEIRD GENE $name BEST GENE $bestGene $bestRNA CURRENT GENE $name $rna{$name}\n";
				$bestGene = $name;
				$bestRNA = $rna{$name};
			}
			else {
				next;
			}
		}
	}
	my $rank;
	if ($fold =~ /Inf/i or $fold =~ /^\-?99$/) {
		$rank = $meanCD;
	}
	else {
		$rank = $meanCD * $fold;
	}
	$bestRNA = $bestRNA == 0 ? 0 : int($bestRNA * 100)/100;
	if (not defined($bestGene)) {
		print $out1 "$chr2\t$start2\t$end2\t$names\t0\t$fold\t$rank\n";
	}
	else {
		print $out1 "$chr2\t$start2\t$end2\t$bestGene\t$bestRNA\t$fold\t$rank\n";
	}
}
close $in2;

close $out1;
