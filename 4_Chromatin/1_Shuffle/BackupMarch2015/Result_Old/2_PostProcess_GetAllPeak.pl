#!/usr/bin/perl
# This script filter out dripc peaks that has less than 125 shuffles
# The difference with 2_PostProcess.pl is this prints out all peaks, not just 125

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <dripc_promoter.txt>\n" unless @ARGV;

my $desiredPeak = 125;
my ($folder, $dripcName) = mitochy::getFilename($input, "folder");
my $randomTooBig    = 0;
my $LessThanDesiredPeakNumberd = 0;
my $totalPeak_withLowDRIPcValue       = 0;
my %data;
my %totalPeak;

# Parse .txt file
open (my $in5, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in5>) {
        chomp($line);
        next if $line =~ /#/;
        my ($chr, $start, $end, $value, $name, $valueOrig, $strand, $info) = split("\t", $line);
	$totalPeak{$name}++;

        # Next if DRIPc value too big (more than 2)
        if ($value > 2) {
                $randomTooBig ++;
                next;
        }
        else {
                $data{$name}{count} ++;
                my $valuez = "$chr\t$start\t$end\t$name\t$valueOrig\t$strand\t$info";
                push(@{$data{$name}{line}}, $valuez);
        }
}
close $in5;

# Get dripc peaks that has $desiredPeak number of shuffles
my $mean = 0;
my %good;
open (my $out5, ">", "$dripcName.tempshuf") or die "Cannot write to $dripcName.tempshuf: $!\n";
foreach my $name (keys %data) {
        $LessThanDesiredPeakNumberd ++ if $data{$name}{count} < $desiredPeak;
        $mean += $data{$name}{count};
        $totalPeak_withLowDRIPcValue ++;
        if ($data{$name}{count} >= $desiredPeak) {
                for (my $i = 0; $i < $data{$name}{count}; $i++) {
                        push(@{$good{$name}}, "$data{$name}{line}[$i]");
                }
        }
}

my %randomI;
foreach my $name (keys %good) {
	for (my $i = 0; $i < @{$good{$name}}; $i++) {
		print $out5 "$good{$name}[$i]\n";
	}
}

close $out5;

system("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcName.tempshuf > $dripcName.shuffled && rm $dripcName.tempshuf");

my $totalPeakTotal = (keys %totalPeak);
printf "
Total DRIPc peak in this file = $totalPeakTotal
Total DRIPc peak after filtering for high DRIPc value = $totalPeak_withLowDRIPcValue (%.2f %%)
Desired Number of Peak = $desiredPeak
Not used = $LessThanDesiredPeakNumberd
", int($totalPeak_withLowDRIPcValue / $totalPeakTotal * 10000)/100;
