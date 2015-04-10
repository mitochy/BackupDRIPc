#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

mkdir "Intersect" if not -d "Intersect";
my @bed = <../Chromatin/peak/*.bed>;
my @feature = qw(promoter terminal genebody);

my %done;
open (my $in, "<", "FoldIntersectPeak.txt") or die;
while (my $line = <$in>) {
	my ($feature, $chip) = split("\t", $line);
	$done{$feature}{$chip} = 1;
}
close $in;

open (my $out, ">>", "Intersect/Result.txt") or die;
open (my $out2, ">>", "Intersect/ResultPeak.txt") or die;
open (my $out3, ">>", "FoldIntersectPeak.txt") or die;
print $out "Feature\tChIP\tTotalOrig\tTotalShuf\tTotalGenomic\tTotalChIP\tFeatOrig\tFeatShuf\tFeatChIP\tFeatChIPOrig\tFeatChIPShuf\tPercChIPOrig\tPercChIPShuf\tEnrichment\n";
print $out2 "Feature\tChIP\tTotalOrig\tTotalShuf\tTotalChIP\tOrigChIP\tShufChIP\tOrigChIP_Perc\tShufChIP_Perc\tPeakEnrichment\tOrigChIPRev\tShufChIPRev\tOrigChIPRev_Perc\tShufChIPRev_Perc\tPeakEnrichmentRev\n";
foreach my $feature (@feature) {
	my $orig = "dripc_$feature\_orig.peak";
	my $shuf = $orig; $shuf =~ s/orig/shuf/;
	my $genomic = "../../1_Location/genomic_$feature.bed";
	`bedtools intersect -u -a $genomic -b $orig > Intersect/$feature\_orig.bed` if not -e "Intersect/$feature\_orig.bed";
	`bedtools intersect -u -a $genomic -b $shuf > Intersect/$feature\_shuf.bed` if not -e "Intersect/$feature\_shuf.bed";
	my ($TotalOrig) = `wc -l $orig` =~ /^(\d+) /;
	my ($TotalShuf) = `wc -l $shuf` =~ /^(\d+) /;
	my ($TotalGenomic) = `wc -l $genomic` =~ /^(\d+) /;
	my ($FeatOrig) = `wc -l Intersect/$feature\_orig.bed` =~ /^(\d+) /;
	my ($FeatShuf) = `wc -l Intersect/$feature\_shuf.bed` =~ /^(\d+) /;

	foreach my $bed (@bed) {
		my ($chip) = mitochy::getFilename($bed);
		next if defined($done{$feature}{$chip}) and $done{$feature}{$chip} eq 1;
		`bedtools intersect -u -b $genomic -a $bed > Intersect/$chip\_at_$feature.bed` if not -e "Intersect/$chip\_at_$feature.bed";
		my ($ChIPAtFeat) = `wc -l Intersect/$chip\_at_$feature.bed` =~ /^(\d+) /;
		print "Processing $feature $chip\n";
		my ($TotalChIP) = `wc -l $bed` =~ /^(\d+) /;
		`bedtools intersect -u -a $orig -b Intersect/$chip\_at_$feature.bed > Intersect/$orig\_$chip.bed` if not -e "Intersect/$orig\_$chip.bed";
		`bedtools intersect -u -a $shuf -b Intersect/$chip\_at_$feature.bed > Intersect/$shuf\_$chip.bed` if not -e "Intersect/$shuf\_$chip.bed";
		`bedtools intersect -u -b $orig -a Intersect/$chip\_at_$feature.bed > Intersect/$orig\_rev\_$chip.bed` if not -e "Intersect/$orig\_rev\_$chip.bed";
		`bedtools intersect -u -b $shuf -a Intersect/$chip\_at_$feature.bed > Intersect/$shuf\_rev\_$chip.bed` if not -e "Intersect/$shuf\_rev\_$chip.bed";
		`bedtools intersect -u -a $genomic -b $bed > Intersect/$feature\_$chip.bed` if not -e "Intersect/$feature\_$chip.bed";
		`bedtools intersect -u -a Intersect/$feature\_orig.bed -b Intersect/$feature\_$chip.bed > Intersect/$feature\_$chip\_orig.bed` if not -e "Intersect/$feature\_$chip\_orig.bed";
		`bedtools intersect -u -a Intersect/$feature\_shuf.bed -b Intersect/$feature\_$chip.bed > Intersect/$feature\_$chip\_shuf.bed` if not -e "Intersect/$feature\_$chip\_shuf.bed";
		my ($FeatChIP) = `wc -l Intersect/$feature\_$chip.bed` =~ /^(\d+) /;
		my ($FeatChIPOrig) = `wc -l Intersect/$feature\_$chip\_orig.bed` =~ /^(\d+) /;
		my ($FeatChIPShuf) = `wc -l Intersect/$feature\_$chip\_shuf.bed` =~ /^(\d+) /;
		my ($OrigChIP) = `wc -l Intersect/$orig\_$chip.bed` =~ /^(\d+) /;
		my ($ShufChIP) = `wc -l Intersect/$shuf\_$chip.bed` =~ /^(\d+) /;
		my ($OrigChIPRev) = `wc -l Intersect/$orig\_rev\_$chip.bed` =~ /^(\d+) /;
		my ($ShufChIPRev) = `wc -l Intersect/$shuf\_rev\_$chip.bed` =~ /^(\d+) /;
		# Feature with Peak intersect
		my $percOrig = int(($FeatChIPOrig+1)/($FeatOrig+1)*1000)/10;
		my $percShuf = int(($FeatChIPShuf+1)/($FeatShuf+1)*1000)/10;
		my $enrichment = int((($FeatChIPOrig+1)/($FeatOrig+1)) / (($FeatChIPShuf+1)/($FeatShuf+1)) * 100) / 100;
		# Peak with peak intersect
		my $percOrigChIP = int(($OrigChIP+1)/($TotalOrig+1)*1000)/10;
		my $percShufChIP = int(($ShufChIP+1)/($TotalShuf+1)*1000)/10;
		my $percOrigChIPRev = int(($OrigChIPRev+1)/($ChIPAtFeat+1)*1000)/10;
		my $percShufChIPRev = int(($ShufChIPRev+1)/($ChIPAtFeat+1)*1000)/10;
		my $enrichmentPeak = $percShufChIP == 0 ? 0 : int(($percOrigChIP / $percShufChIP)*100)/100;
		my $enrichmentPeakRev = $percShufChIPRev == 0 ? 0 : int(($percOrigChIPRev / $percShufChIPRev)*100)/100;
		print $out "$feature\t$chip\t$TotalOrig\t$TotalShuf\t$TotalGenomic\t$TotalChIP\t$FeatOrig\t$FeatShuf\t$FeatChIP\t$FeatChIPOrig\t$FeatChIPShuf\t$percOrig\t$percShuf\t$enrichment\n";
		print $out2 "$feature\t$chip\t$TotalOrig\t$TotalShuf\t$TotalChIP\t$OrigChIP\t$ShufChIP\t$percOrigChIP\t$percShufChIP\t$enrichmentPeak\t$OrigChIPRev\t$ShufChIPRev\t$percOrigChIPRev\t$percShufChIPRev\t$enrichmentPeakRev\n";
		
		# Calculate intersect percentage and p value
		my %data;
		open (my $in1, "<", "$orig") or die;
		while (my $line = <$in1>) {
			chomp($line);
		        my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
			my ($drip) = $info =~ /DRIP=(.+);ORIG=/; die if not defined($drip);
			$data{$drip}{orig} = 1;
		}
		close $in1;

		open (my $in2, "<", "Intersect/$shuf\_$chip.bed") or die;
		while (my $line = <$in2>) {
			chomp($line);
		        my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
			my ($drip) = $info =~ /DRIP=(.+);ORIG=/; die if not defined($drip);
			push(@{$data{$drip}{shuf}}, 1);
		}
		my $count = 0;
		my $origTotal = keys %data;
		foreach my $drip (keys %data) {
			$count ++;
			#print "$count / $origTotal\n" if $count % 100 == 0;
			for (my $i = 0; $i < 125; $i++) {
				$data{$drip}{shuf}[$i] = 0 if not defined($data{$drip}{shuf}[$i]);
			}
			@{$data{$drip}{shuf}} = shuffle(@{$data{$drip}{shuf}});
		}
		close $in2;

		my ($nge, $nle, $nrun) = (0,0,0);
		my @shuf; my @fold;
		my $origCount = $percOrigChIP;
		for (my $i = 0; $i < 125; $i++) {
			my $shufCount = 0;
			foreach my $drip (keys %data) {
				$shufCount ++ if $data{$drip}{shuf}[$i] == 1;
			}
			$shufCount = int($shufCount / $origTotal*1000)/10;
			push(@shuf, $shufCount);
			my $fold = $shufCount == 0 ? $origCount + 1 : $origCount / $shufCount;
			push(@fold, $fold);
			$nge ++ if $origCount <= $shufCount;
			$nle ++ if $origCount >= $shufCount;
			$nrun ++;
		}
		my $medianShuf = median(@shuf);
		my $medianFold = median(@fold);
		my $foldAll = join("\t", @fold);
		my $pvalHi = int(($nge + 1) / ($nrun + 1) * 1000)/1000;
		my $pvalLo = int(($nle + 1) / ($nrun + 1) * 1000)/1000;
		my $pval = $pvalHi <= 0.05 ? $pvalHi : $pvalLo <= 0.05 ? $pvalLo : 1;
		print $out3 "$feature\t$chip\t$medianFold\t$pval\t$medianFold\t$foldAll\n";
	}
}


sub process_bed {
	my ($input) = @_;
	my %data;

}

sub shuffle {
        my (@value) = @_;
        #print "Before: @value\n";
        for (my $i = 0; $i < 1000; $i++) {
                my $rand1 = int(rand(@value));
                my $rand2 = int(rand(@value));
                my $val1 = $value[$rand1];
                my $val2 = $value[$rand2];
                $value[$rand1] = $val2;
                $value[$rand2] = $val1;
        }
        #print "After: @value\n";
        return(@value);
}

sub median {
        my (@value) = @_;
        @value = sort {$a <=> $b} (@value);
        my $median = $value[int(@value/2)];
        return($median);
}
