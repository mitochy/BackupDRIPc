#!/usr/bin/perl
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use strict; use warnings; use mitochy;

my ($input1, $input2) = @ARGV;
die "usage: $0 <input1> <input2>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
my ($folder2, $fileName2) = mitochy::getFilename($input2, "folder");
my %data;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($cell, $tf, $junk, $junk1, $val, $junk3, $val2, $junk4, $val3) = split("\t", $line);
	$data{$cell}{$tf}{DRIP1}{1} = $val if $junk eq "DRIP1";
	$data{$cell}{$tf}{DRIP2}{1} = $val if $junk eq "DRIP2";
	$data{$cell}{$tf}{DRIP1}{shuf} = $val2 if $junk eq "DRIP1";
	$data{$cell}{$tf}{DRIP2}{shuf} = $val2 if $junk eq "DRIP2";

}
close $in1;

open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($cell, $tf, $junk, $junk1, $val, $junk3, $val2, $junk4, $val3) = split("\t", $line);
	$data{$cell}{$tf}{DRIP1}{2} = $val if $junk eq "DRIP1";
	$data{$cell}{$tf}{DRIP2}{2} = $val if $junk eq "DRIP2";
	$data{$cell}{$tf}{DRIP1}{peak} = $val3 if $junk eq "DRIP1";
	$data{$cell}{$tf}{DRIP2}{peak} = $val3 if $junk eq "DRIP2";
}
close $in2;

my %stat;
foreach my $cell (keys %data) {
	foreach my $tf (keys %{$data{$cell}}) {
		my $fold  = int(100*$data{$cell}{$tf}{DRIP1}{2}**3 / $data{$cell}{$tf}{DRIP1}{1}**2);
		my $fold2 = int(100*$data{$cell}{$tf}{DRIP2}{2}**3 / $data{$cell}{$tf}{DRIP2}{1}**2);
		#my $fold2 = int(10000*$data{$cell}{$tf}{DRIP1}{2} / $data{$cell}{$tf}{DRIP1}{1})/100;
		$stat{"$cell\t$tf"}{fold} = $fold;
		$stat{"$cell\t$tf"}{fold2} = $fold2;
		$stat{"$cell\t$tf"}{val11} = $data{$cell}{$tf}{DRIP1}{1};
		$stat{"$cell\t$tf"}{val12} = $data{$cell}{$tf}{DRIP1}{2};
		$stat{"$cell\t$tf"}{val21} = $data{$cell}{$tf}{DRIP2}{1};
		$stat{"$cell\t$tf"}{val22} = $data{$cell}{$tf}{DRIP2}{2};
		$stat{"$cell\t$tf"}{shuf1} = $data{$cell}{$tf}{DRIP1}{shuf};
		$stat{"$cell\t$tf"}{shuf2} = $data{$cell}{$tf}{DRIP2}{shuf};
		$stat{"$cell\t$tf"}{peak1} = $data{$cell}{$tf}{DRIP1}{peak};
		$stat{"$cell\t$tf"}{peak2} = $data{$cell}{$tf}{DRIP2}{peak};
		
	}
}

foreach my $name (sort {$stat{$b}{fold} <=> $stat{$a}{fold}} keys %stat) {
	print RED "$name\t$stat{$name}{fold}\t$stat{$name}{val11}\t$stat{$name}{val12}\t$stat{$name}{shuf1}\t$stat{$name}{peak1}\t$stat{$name}{fold2}\t$stat{$name}{val21}\t$stat{$name}{val22}\t$stat{$name}{shuf2}\t$stat{$name}{peak2}\n" and next if $name =~ /ctcf/i;
	print BLUE "$name\t$stat{$name}{fold}\t$stat{$name}{val11}\t$stat{$name}{val12}\t$stat{$name}{shuf1}\t$stat{$name}{peak1}\t$stat{$name}{fold2}\t$stat{$name}{val21}\t$stat{$name}{val22}\t$stat{$name}{shuf2}\t$stat{$name}{peak2}\n" and next if $name =~ /\tpol/i;
	print GREEN "$name\t$stat{$name}{fold}\t$stat{$name}{val11}\t$stat{$name}{val12}\t$stat{$name}{shuf1}\t$stat{$name}{peak1}\t$stat{$name}{fold2}\t$stat{$name}{val21}\t$stat{$name}{val22}\t$stat{$name}{shuf2}\t$stat{$name}{peak2}\n" and next if $name =~ /rad21/i;
	print YELLOW "$name\t$stat{$name}{fold}\t$stat{$name}{val11}\t$stat{$name}{val12}\t$stat{$name}{shuf1}\t$stat{$name}{peak1}\t$stat{$name}{fold2}\t$stat{$name}{val21}\t$stat{$name}{val22}\t$stat{$name}{shuf2}\t$stat{$name}{peak2}\n" and next if $name =~ /runx/i;
	print MAGENTA "$name\t$stat{$name}{fold}\t$stat{$name}{val11}\t$stat{$name}{val12}\t$stat{$name}{shuf1}\t$stat{$name}{peak1}\t$stat{$name}{fold2}\t$stat{$name}{val21}\t$stat{$name}{val22}\t$stat{$name}{shuf2}\t$stat{$name}{peak2}\n" and next if $name =~ /p300/i;
	print "$name\t$stat{$name}{fold}\t$stat{$name}{val11}\t$stat{$name}{val12}\t$stat{$name}{shuf1}\t$stat{$name}{peak1}\t$stat{$name}{fold2}\t$stat{$name}{val21}\t$stat{$name}{val22}\t$stat{$name}{shuf2}\t$stat{$name}{peak2}\n" 

}

