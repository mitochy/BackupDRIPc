#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input, $rnaFile, $genomicFile) = @ARGV;
die "usage: $0 <dripc_promoter.shuffled> <NT2.rpkm> <genomic_genebody or hg19_gencode_promoter.bed>\n" unless @ARGV == 3;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";

my $desiredPeak = 125;
my %count;
my %stat;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $value, $strand, $info) = split("\t", $line);
	my ($DRIPCHR, $DRIPstart, $DRIPend) = $info =~ /DRIP=(\w+),(\d+),(\d+)/;
	$name = "$DRIPCHR\_$DRIPstart\_$DRIPend";
	$count{$name}++;

	$stat{$name}{total} ++;

	# CHR unbiased
	$stat{$name}{chr}{$chr}++;

	# DRIP length correct
	my ($gene) = $info =~ /TWIN=(\w+\.\d+)/;
	my ($twinstart, $twinend) = $info =~ /TWIN=ENS\w+\.\d+,\w+,(\d+),(\d+)$/;
	my $dist = $end - $start;
	my $DRIPdist = $DRIPend - $DRIPstart;
	print "DRIPLENGTH ERROR: ($gene $chr $start $end, $name $DRIPCHR $DRIPstart $DRIPend, $dist vs $DRIPdist)\n" if ($DRIPend - $DRIPstart != $end - $start);
		
	# Gene unbiased
	$stat{$name}{gene}{$gene} ++;

	# Strand unbiased
	$stat{$name}{strand}{$strand}++;

	# Startpos is unbiased
	my $startpos = $strand eq "+" ? ($start - $twinstart) / ($twinend - $twinstart) : ($twinend - $end) / ($twinend - $twinstart);
	push(@{$stat{$name}{startpos}}, $startpos);
	
}
close $in;

my %count_dist;
foreach my $name (keys %count) {
	my $count = $count{$name};
	$count_dist{$count} ++;
}

my $max_count = 0;
foreach my $count (sort {$b <=> $a} keys %count_dist) {
	$max_count = $count;
	last;
}
$max_count = $desiredPeak;
my $printed = 0;
my $total_count = 0;
foreach my $name (keys %count) {
	my $count = $count{$name};
	print "$name has count $count while max is $max_count\n" if $count < $max_count and $printed < 10;
	$printed ++ if $count < $max_count;
	$total_count ++ if $count < $max_count;
}
print "Has less than $max_count: $total_count\n";
my $drip_peak_number = (keys %count);
print "$input: Total peak = $drip_peak_number, shuffled $max_count times\n";

foreach my $name (keys %stat) {
	foreach my $chr (sort keys %{$stat{$name}{chr}}) {
		print "BIASEDCHR: $name\t$chr\t$stat{$name}{chr}{$chr} / $stat{$name}{total}\n" if $stat{$name}{chr}{$chr} / $stat{$name}{total} > 0.3;
	}
}

foreach my $name (keys %stat) {
	foreach my $gene (sort keys %{$stat{$name}{gene}}) {
		print "BIASEDgene: $name\t$gene\t$stat{$name}{gene}{$gene} / $stat{$name}{total}\n" if $stat{$name}{gene}{$gene} / $stat{$name}{total} > 0.3;
	}
}


foreach my $name (keys %stat) {
	foreach my $strand (sort keys %{$stat{$name}{strand}}) {
		print "BIASEDstrand: $name\t$strand\t$stat{$name}{strand}{$strand}\n" if $stat{$name}{strand}{$strand} / $stat{$name}{total} > 0.7;
	}
}

