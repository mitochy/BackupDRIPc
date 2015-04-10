#!/usr/bin/perl

use strict; use warnings; use mitochy; use Cache::FileCache;

my ($input1, $input2) = @ARGV;
die "usage: $0 <DRIPc1> <DRIPc2>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my $cache = new Cache::FileCache;
$cache->set_cache_root("/data/mitochi/Work/Cache/");
open (my $out, ">", "Combined.wig") or die;
my ($index) = "/data/Genome/hg19.fa.fai";
open (my $in, "<", $index) or die;
while (my $line = <$in>) {
	chomp($line);
	my ($chr, $end) = split("\t", $line);
	print $out "variableStep chrom=$chr span=10\n";
	for (my $i = 0; $i < int($end / 1000000); $i++) {
		printf "Procesing $chr: %.2f %%\n", int($i * 1000000 / $end * 10000)/100;
		my %result;
		my (%wig1, %wig2);
		my ($fileName1, $fileName2) = (mitochy::getFilename($input1, "fullname"), mitochy::getFilename($input2, "fullname"));
		my $DRIPc1 = $cache->get("$fileName1\.$chr\.$i.cache");
		my $DRIPc2 = $cache->get("$fileName2\.$chr\.$i.cache");
		#print "DOESNT WORK $fileName1.$chr.$i.cache\n" if not defined($DRIPc1);
		#print "WORK $fileName1.$chr.$i.cache\n" if defined($DRIPc1);
		%wig1 = %{$DRIPc1} if defined($DRIPc1);
		%wig2 = %{$DRIPc2} if defined($DRIPc2);
		for (my $j = 0; $j < int($end / 1000); $j++) {
			foreach my $pos (keys %{$wig1{$chr}{$i}{$j}}) {
				my $value1 = $wig1{$chr}{$i}{$j}{$pos};
				my $value2 = defined($wig2{$chr}{$i}{$j}{$pos}) ? $wig2{$chr}{$i}{$j}{$pos} : 0;
				$result{$pos} = $value1 + $value2;
			}
			foreach my $pos (keys %{$wig2{$chr}{$i}{$j}}) {
				my $value2 = $wig2{$chr}{$i}{$j}{$pos};
				my $value1 = defined($wig1{$chr}{$i}{$j}{$pos}) ? $wig1{$chr}{$i}{$j}{$pos} : 0;
				$result{$pos} = $value1 + $value2;
			}
		}
		foreach my $pos (sort {$a <=> $b} keys %result) {
			print $out "$pos\t$result{$pos}\n";
		}

	}
}
__END__

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
my ($line1, $line2);
my ($chr1, $chr2);
my ($pos1, $pos2);
my ($val1, $val2);
while (1) {
	if (not defined($line1)) {
		$line1 = <$in1>;
		$line2 = <$in2>;
	}
	$line1 = <$in1> if $line1 =~ /track/;
	$line2 = <$in2> if $line2 =~ /track/;
	($chr1) = $line1 =~ /chrom=(\w+) span/ and $line1 = <$in1> if ($line1 =~ /variable/);
	($chr2) = $line2 =~ /chrom=(\w+) span/ and $line2 = <$in2> if ($line2 =~ /variable/);
	
	($pos1, $val1) = split("\t", $line1) if ($line1 =~ /^\d+\t\d+/);
	($pos2, $val2) = split("\t", $line2) if ($line2 =~ /^\d+\t\d+/);
	
	if ($pos1 == $pos2 and $chr1 eq $chr2) {
		
	}

}
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}
close $in1;

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;
