#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input2, $gencode, $rnaseq) = @ARGV;
die "usage: $0 <dripc_promoter.name> <mm9_gencode.bed> <rnaseq.rpkm>\n" unless @ARGV == 3;

my %data;
open (my $in1, "<", $gencode) or die "Cannot read from $gencode: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	$data{$name} = $strand;
}
close $in1;

my %rna;
open (my $in0, "<", $rnaseq) or die "Cannot read from $rnaseq: $!\n";
while (my $line = <$in0>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name, $val) = split("\t", $line);
	$rna{$name} = $val;
}
close $in0;


my ($folder2, $fileName2) = mitochy::getFilename($input2, "folder");
open (my $out2, ">", "$fileName2.out") or die "Cannot write to $fileName2.out: $!\n";

open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $junk, $zero, $strand, $names) = split("\t", $line);
	my @names = split(";", $names);
	if (@names == 1) {
		$strand = $data{$names[0]};
		print $out2 "$chr\t$start\t$end\t$junk\t$zero\t$strand\t$names\n";
	}
	else {
		my ($best_name, $best_rna);
		foreach my $name (@names) {
			if (not defined($best_name)) {
				$best_name = $name;
				$best_rna  = defined($rna{$name}) ? $rna{$name} : 0;
			}
			else {
				if (defined($rna{$name})) {
					$best_name = $name if $best_rna < $rna{$name};
					$best_rna  = $rna{$name} if $best_rna < $rna{$name};
				}
			}
		}
		$strand = $data{$best_name};die "Died at $best_name\n" if not defined($strand);
		print $out2 "$chr\t$start\t$end\t$junk\t$zero\t$strand\t$best_name\n";
	}
}
close $in2;
close $out2;

system("mv $fileName2.out $fileName2.name");
