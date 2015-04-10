#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my $NT2_rnaFile = "/data/mitochi/Work/Project/DRIPc/data/NT2.rpkm";
my $E14_rnaFile = "/data/mitochi/Work/Project/DRIPc/data/E14.rpkm";
my $X3T3_rnaFile = "/data/mitochi/Work/Project/DRIPc/data/3T3.rpkm";
my $Fib_rnaFile = "/data/mitochi/Work/Project/DRIPc/data/Fibro.rpkm";
my %NT2_rna = %{parse_rna($NT2_rnaFile)};
my %E14_rna = %{parse_rna($E14_rnaFile)};
my %X3T3_rna = %{parse_rna($X3T3_rnaFile)};
my %Fib_rna = %{parse_rna($Fib_rnaFile)};
my %orth;
open (my $in2, "<", "human_mouse.orth") or die;
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /Ensembl/;
	my ($hGID, $hTID, $mGID, $id) = split("\t", $line);
	next if not defined($mGID) or $mGID !~ /ENS/;
	if (not defined($orth{human}{$hTID}) or $orth{human}{$hTID}{id} < $id) {
		$orth{human}{$hTID}{name} = $mGID;
		$orth{human}{$hTID}{id}   = $id;
		$orth{mouse}{$mGID}{name} = $hTID;
	}
}
close $in2;

open (my $in3, "<", "mouse.geneinfo") or die;
while (my $line = <$in3>) {
	chomp($line);
	next if $line =~ /Ensembl/;
	my ($mGID, $mTID) = split("\t", $line);
	if (defined($orth{mouse}{$mGID})) {
		my $hTID = $orth{mouse}{$mGID}{name};
		push(@{$orth{human}{$hTID}{orth}}, $mTID);
	}
}
close $in3;
my %data;

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $val, $val2, $strand, $chr2, $start2, $end2, $gene) = split("\t", $line);
	($gene) = $gene =~ /^(.+)\.\d+/ if $gene =~ /.+\.\d+/;
	my $name = "$chr\t$start\t$end\t$val\t$val2\t$strand";
	push(@{$data{$chr}{$start}{$end}{name}}, $gene);
	$data{$chr}{$start}{$end}{line} = $name;
}
close $in1;

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
foreach my $chr (sort keys %data) {
	foreach my $start (sort {$a <=> $b} keys %{$data{$chr}}) {
		foreach my $end (sort {$a <=> $b} keys %{$data{$chr}{$start}}) {
			my $line = $data{$chr}{$start}{$end}{line};
			my $bestNT2rna = 0;
			my $bestNT2gene = 0;
			my $bestFibrna = 0;
			my $bestFibgene = 0;
			my $bestE14rna = 0;
			my $bestE14gene = 0;
			my $best3T3rna = 0;
			my $best3T3gene = 0;
			foreach my $gene (@{$data{$chr}{$start}{$end}{name}}) {
				my $NT2rna = defined($NT2_rna{$gene}) ? $NT2_rna{$gene} : 0;
				if ($NT2rna > $bestNT2rna) {
					$bestNT2rna = $NT2rna;	
					$bestNT2gene = $gene;
				}
				my $Fibrna = defined($Fib_rna{$gene}) ? $Fib_rna{$gene} : 0;
				if ($Fibrna > $bestFibrna) {
					$bestFibrna = $Fibrna;	
					$bestFibgene = $gene;
				}
				foreach my $mouse (@{$orth{human}{$gene}{orth}}) {
					#die "MOUSE $mouse RNA $E14_rna{$gene}\n";
					my $E14rna = defined($E14_rna{$mouse}) ? $E14_rna{$mouse} : 0;
					if ($E14rna > $bestE14rna) {
						$bestE14rna = $E14rna;	
						$bestE14gene = $gene;
					}
					my $X3T3rna = defined($X3T3_rna{$mouse}) ? $X3T3_rna{$mouse} : 0;
					if ($X3T3rna > $best3T3rna) {
						$best3T3rna = $X3T3rna;	
						$best3T3gene = $gene;
					}
				}
			}
			my ($chr, $start, $end, $val, $val2, $strand) = split("\t", $line);
			next if $bestNT2rna == 0;
			next if $bestFibrna == 0;
			#next if $bestE14rna == 0;
			#next if $best3T3rna == 0;
			print $out1 "$chr\t$start\t$end\t$val\t$val2\t$strand\t$bestNT2rna\t$bestFibrna\t$bestE14rna\t$best3T3rna\t$bestNT2gene\t$bestFibgene\t$bestE14gene\t$best3T3gene\n";
		}
	
	}
}
close $out1;

sub parse_rna {
        my ($input) = @_;
        my %data;
        open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
        while (my $line = <$in>) {
                chomp($line);
                next if $line =~ /#/;
                my ($gene, $val) = split("\t", $line);
		($gene) = $gene =~ /^(.+)\.\d+/ if $gene =~ /.+\.\d+/;
                $data{$gene} = $val;
		#print "GENE $gene VAL $val\n" if $gene =~ /ENSMUST00000030949/;
        }
        close $in;
        return(\%data);
}
