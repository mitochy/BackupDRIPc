#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my %rna = %{parse_rna("../data/NT2.rpkm")};

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
open (my $out1, ">", "$fileName1\_zero.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $out2, ">", "$fileName1\_low.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $out3, ">", "$fileName1\_med.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $out4, ">", "$fileName1\_high.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $out5, ">", "$fileName1\_super.out") or die "Cannot write to $fileName1.out: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $names, $val, $strand) = split("\t", $line);
	my @names = split(";", $names);
	my ($bestname, $bestrna);
	foreach my $name (@names) {
		my $rna = defined($rna{$name}) ? $rna{$name} : 0;
		if (defined($bestname)) {
			$bestname = $name if $bestrna < $rna;
			$bestrna = $rna if $bestrna < $rna;
		}
		else {
			$bestname = $name;
			$bestrna = $rna;
		}
	}
	if ($bestrna >= 200) {
		print $out5 "$chr\t$start\t$end\t$bestname\t$bestrna\t$strand\n";
	}
	elsif ($bestrna >= 100) {
		print $out4 "$chr\t$start\t$end\t$bestname\t$bestrna\t$strand\n";
	}
	elsif ($bestrna >= 50) {
		print $out3 "$chr\t$start\t$end\t$bestname\t$bestrna\t$strand\n";
	}
	elsif ($bestrna >= 10) {
		print $out2 "$chr\t$start\t$end\t$bestname\t$bestrna\t$strand\n";
	}
	else {
		print $out1 "$chr\t$start\t$end\t$bestname\t$bestrna\t$strand\n";
	}
}
close $in1;
close $out1;


sub parse_rna {
        my ($input) = @_;
        my %data;
        open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
        while (my $line = <$in>) {
                chomp($line);
                next if $line =~ /#/;
                my ($gene, $val) = split("\t", $line);
                $data{$gene} = $val;
        }
        close $in;
        return(\%data);
}
