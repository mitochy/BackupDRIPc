#!/usr/bin/perl

use strict; use warnings; use mitochy;

my $orthTable = "UCSC_hg19mm10_orthtable.txt";
my $humanFile = "human_ordered.bed";
my $mouseFile = "mouse_ordered.bed";

open (my $out, ">", "RESULT.txt") or die "Cannot write to RESULT.txt: $!\n";

my %orth;
my $linecount = 0;
open (my $in1, "<", $orthTable) or die "Cannot read from $orthTable: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	$linecount++;
	next if $line =~ /#/;
	my @arr = split("\t", $line);
	$orth{count}{human}{$arr[0]} = $linecount;
	$orth{count}{mouse}{$arr[1]} = $linecount;
	$orth{revcount}{$linecount}{human} = $arr[0];
	$orth{revcount}{$linecount}{mouse} = $arr[1];
}
close $in1;

my %data;
open (my $in2, "<", $humanFile) or die "Cannot read from $humanFile: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $DRIP, $strand) = split("\t", $line);
	next if not defined($orth{count}{human}{$name});
	my $count = $orth{count}{human}{$name};
	$data{DRIP}{$count}{human}{val} = $DRIP;
	$data{DRIP}{$count}{human}{line} = $line;
	#print "HUMAN NAME $name = $count\n";
}
close $in2;

open (my $in3, "<", $mouseFile) or die "Cannot read from $mouseFile: $!\n";
while (my $line = <$in3>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $DRIP, $strand) = split("\t", $line);
	next if not defined($orth{count}{mouse}{$name});
	my $count = $orth{count}{mouse}{$name};
	$data{DRIP}{$count}{mouse}{val} = $DRIP;
	$data{DRIP}{$count}{mouse}{line} = $line;
}
close $in3;

foreach my $count (sort {$a <=> $b} keys %{$data{DRIP}}) {
	my $humanDRIP = $data{DRIP}{$count}{human}{val};
	my $mouseDRIP = $data{DRIP}{$count}{mouse}{val};
	next if not defined($humanDRIP) or not defined($mouseDRIP);
	my $humanline = $data{DRIP}{$count}{human}{line};
	my $mouseline = $data{DRIP}{$count}{mouse}{line};
	my $DRIPmult = $humanDRIP * $mouseDRIP;
	$data{FINAL}{$DRIPmult} = "$humanline\t$mouseline";
}

print $out "HUMAN_chr\tstart\tend\tgene\tDRIP\tstrand\tMOUSE_chr\tstart\tend\tgene\tDRIP\tstrand\tDRIP_BOTH\n";
foreach my $DRIPmult (sort {$b <=> $a} keys %{$data{FINAL}}) {
	print $out "$data{FINAL}{$DRIPmult}\t$DRIPmult\n";
}
close $out;
