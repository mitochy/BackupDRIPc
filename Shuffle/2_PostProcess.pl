#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <dripc_promoter.txt>\n" unless @ARGV;

my $desiredPeak = 125;
my ($folder, $dripcName) = mitochy::getFilename($input, "folder");
my $randomTooBig    = 0;
my $lessThanThousand = 0;
my $totalPeak       = 0;
my %data;
my %totalPeak;
open (my $in5, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in5>) {
        chomp($line);
        next if $line =~ /#/;
        my ($chr, $start, $end, $value, $name, $valueOrig, $strand, $info) = split("\t", $line);
	if ($value !~ /^\-?\d+\.?\d*$/) {
	        ($chr, $start, $end, $name, $valueOrig, $strand, $info) = split("\t", $line);
		$value = 0; 
	}
	
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

my $mean = 0;
my %good;
open (my $out5, ">", "$dripcName.tempshuf") or die "Cannot write to $dripcName.tempshuf: $!\n";
foreach my $name (keys %data) {
        $lessThanThousand ++ if $data{$name}{count} < $desiredPeak;
        $mean += $data{$name}{count};
        $totalPeak ++;
        if ($data{$name}{count} >= $desiredPeak) {
                for (my $i = 0; $i < $data{$name}{count}; $i++) {
                        push(@{$good{$name}}, "$data{$name}{line}[$i]");
                }
        }
}

my %randomI;
foreach my $name (keys %good) {
	while ((keys %{$randomI{$name}}) < $desiredPeak) {
		my $random = int(rand(@{$good{$name}}));
		$randomI{$name}{$random} = 1;
	}
}

foreach my $name (keys %randomI) {
	foreach my $random (keys %{$randomI{$name}}) {
		print $out5 "$good{$name}[$random]\n";
	}
}

close $out5;

system("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcName.tempshuf > $dripcName.shuffled && rm $dripcName.tempshuf");

my $totalPeakTotal = (keys %totalPeak);
print "Desired Peak = $desiredPeak\nNot used = $lessThanThousand\nTotal Peak = $totalPeak (total = $totalPeakTotal)\n";
