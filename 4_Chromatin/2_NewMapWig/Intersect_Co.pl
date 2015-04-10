#!/usr/bin/perl

use strict; use warnings; use mitochy; use Thread; use Thread::Queue;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

mkdir "Intersect" if not -d "Intersect";
my @feature = qw(promoter terminal genebody);

my %done;
#my @chips = qw(CTCF Rad21 Znf143 H3K4me1 P300);
#my @chips = qw(H3K4me3 H3K4me1 H4K20me1);
open (my $in, "<", "Intersect/ResultPeak.txt") or die;
while (my $line = <$in>) {
	next if $line =~ /Feature\tChIP/;
	my @arr = split("\t", $line);
	my ($feature, $chip, $totalorig, $totalshuf, $origchip, $shufchip) = ($arr[0], $arr[1], $arr[2], $arr[3], $arr[7], $arr[8]);
	$done{$feature}{$chip}{orig} = $origchip;
	$done{$feature}{$chip}{shuf} = $shufchip;
	$done{$feature}{$chip}{totalorig} = $totalorig;
	$done{$feature}{$chip}{totalshuf} = $totalshuf;
}
close $in;


my $Q = new Thread::Queue;
foreach my $feature (@feature) {
	print "Processing $feature\n";
	my @used;
	my @bed = <Intersect/dripc_$feature\_orig.peak_*.bed>;
	for (my $i = 0;$i < @bed; $i++) {
		my $bed1 = $bed[$i];
		my ($chip1) = $bed1 =~ /dripc_$feature\_orig.peak_(\w+).bed/; die "Died at $bed1 undef chip\n" unless defined($chip1);
		my $orig1 = $bed1;
		my $shuf1 = $bed1; $shuf1 =~ s/orig/shuf/;
		for (my $j = 0;$j < @bed; $j++) {
			my $bed2 = $bed[$j];
			my ($chip2) = $bed2 =~ /dripc_$feature\_orig.peak_(\w+).bed/; die "Died at $bed2 undef chip\n" unless defined($chip2);
			next if grep(/^$chip1\_$chip2$/, @used);
			next if grep(/^$chip2\_$chip1$/, @used);
			next if $chip1 eq $chip2;
			#next if not grep(/^$chip1$/, @chips);
			#next if not grep(/^$chip2$/, @chips);
			#print "\t$i.$j. $chip1 with $chip2\n";
			my $orig2 = $bed2;
			my $shuf2 = $bed2; $shuf2 =~ s/orig/shuf/;
			my $cmd = "bedtools intersect -u -a $orig1 -b $orig2 > Intersect/$chip1\_$chip2\_$feature\_orig.COINT";# if not -e "Intersect/$chip1\_$chip2\_$feature\_orig.COINT";
			$Q -> enqueue($cmd) if not -e "Intersect/$chip1\_$chip2\_$feature\_orig.COINT";
			$cmd = "bedtools intersect -u -a $shuf1 -b $shuf2 > Intersect/$chip1\_$chip2\_$feature\_shuf.COINT";# if not -e "Intersect/$chip1\_$chip2\_$feature\_shuf.COINT";
			$Q -> enqueue($cmd) if not -e "Intersect/$chip1\_$chip2\_$feature\_shuf.COINT";
			push(@used, "$chip1\_$chip2");
		}
	}
}
$Q->end;
my @threads;
for (my $i = 0; $i < 20; $i++) {
	$threads[$i] = threads->create(\&worker, $i, $Q);
}
for (my $i = 0; $i < 20; $i++) {
	$threads[$i] ->join();
}
open (my $out, ">", "Intersect/ResultPeakCoIntersect.txt") or die;
print $out "Feature\tChIP1\tChIP2\tEnrichment\tOrigEnrichment\tShufEnrichment\tOrigChIP1\tOrigChIP2\tOrigExp\tOrigInt\tShufChIP1\tShufChIP2\tShufExp\tShufInt\n";
foreach my $feature (@feature) {
	print "Processing $feature\n";
	my @used;
	my @bed = <Intersect/dripc_$feature\_orig.peak_*.bed>;
	for (my $i = 0;$i < @bed; $i++) {
		my $bed1 = $bed[$i];
		my ($chip1) = $bed1 =~ /dripc_$feature\_orig.peak_(\w+).bed/; die "Died at $bed1 undef chip\n" unless defined($chip1);
		my $orig1 = $bed1;
		my $shuf1 = $bed1; $shuf1 =~ s/orig/shuf/;
		for (my $j = 0;$j < @bed; $j++) {
			my $bed2 = $bed[$j];
			my ($chip2) = $bed2 =~ /dripc_$feature\_orig.peak_(\w+).bed/; die "Died at $bed2 undef chip\n" unless defined($chip2);
			next if grep(/^$chip1\_$chip2$/, @used);
			next if grep(/^$chip2\_$chip1$/, @used);
			next if $chip1 eq $chip2;
			#next if not grep(/^$chip1$/, @chips);
			#next if not grep(/^$chip2$/, @chips);
			#print "\t$i.$j. $chip1 with $chip2\n";
			my $orig2 = $bed2;
			my $shuf2 = $bed2; $shuf2 =~ s/orig/shuf/;
			`bedtools intersect -u -a $orig1 -b $orig2 > Intersect/$chip1\_$chip2\_$feature\_orig.COINT` if not -e "Intersect/$chip1\_$chip2\_$feature\_orig.COINT";
			`bedtools intersect -u -a $shuf1 -b $shuf2 > Intersect/$chip1\_$chip2\_$feature\_shuf.COINT` if not -e "Intersect/$chip1\_$chip2\_$feature\_shuf.COINT";
			my $totalorig1 = $done{$feature}{$chip1}{totalorig}; # Total DRIPc chip 1 (same as chip2) 
			my $origchip1  = $done{$feature}{$chip1}{orig}; # DRIPc with chip1 (percentage)
			my $totalorig2 = $done{$feature}{$chip2}{totalorig};
			my $origchip2  = $done{$feature}{$chip2}{orig};
			my $origexp    = $origchip1 * $origchip2 / 10000 * 100; # Expected intersect between chip1 and chip2 (percent)
			my ($origcoint)  = `wc -l Intersect/$chip1\_$chip2\_$feature\_orig.COINT` =~ /^(\d+) /; $origcoint = $origcoint * 100 / $totalorig1; # Observed intersect
			#my $origenrich = $origchip1 == 0 ? $origcoint - $origexp : ($origcoint - $origexp) / $origchip1 * 100;
			my $origenrich = $origexp == 0 ? ($origcoint + 1) : (1+$origcoint) / (1+$origexp);

			my $totalshuf1 = $done{$feature}{$chip1}{totalshuf}; # Total DRIPc chip 1 (same as chip2)
			my $shufchip1  = $done{$feature}{$chip1}{shuf}; # DRIPc with chip1 (percentage)
			my $totalshuf2 = $done{$feature}{$chip2}{totalshuf};
			my $shufchip2  = $done{$feature}{$chip2}{shuf};
			my $shufexp    = $shufchip1 * $shufchip2 / 10000 * 100; # Expected intersect between chip1 and chip2 (percent)
			my ($shufcoint)  = `wc -l Intersect/$chip1\_$chip2\_$feature\_shuf.COINT` =~ /^(\d+) /; $shufcoint = $shufcoint * 100 / $totalshuf1; # Observed intersect
			my $shufenrich = $shufexp == 0 ? ($shufcoint + 1) : (1+$shufcoint) / (1+$shufexp);
			#my $shufenrich = $shufchip1 == 0 ? $shufcoint - $shufexp : ($shufcoint - $shufexp) / $shufchip1 * 100;
			#die "totalshuf1 $totalshuf1 shufchip1 $shufchip1 totalshuf2 $totalshuf2 shufchip2 $shufchip2 shufexpected $shufexp shufcoint $shufcoint shufenrich $shufenrich\n" if $shufenrich == 0;
			my $enrichment = int(($origenrich + 1) / (1+$shufenrich)*100)/100;
			printf $out "$feature\t$chip1\t$chip2\t$enrichment\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"
			,$origenrich,$shufenrich,$origchip1,$origchip2,$origexp,$origcoint,$shufchip1,$shufchip2,$shufexp,$shufcoint;
			push(@used, "$chip1\_$chip2");
		}
	}
}

sub worker {
	my ($thread, $queue) = @_;
	while ($queue->pending) {
		my $command = $queue->dequeue;
		`$command`;
	}
	return;

}


__END__






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
