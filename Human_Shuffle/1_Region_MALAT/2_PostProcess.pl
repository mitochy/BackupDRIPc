#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <dripName.txt>\n" unless @ARGV;

my ($folder, $dripcName) = mitochy::getFilename($input, "folder");
($dripcName) = $dripcName =~ /^(\w+)\_Shuffled$/;
# c. Print out values 
my $randomTooBig    = 0;
my $lessThan50 = 0;
my $totalPeak       = 0;
my %data;
open (my $in5, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in5>) {
        chomp($line);
        next if $line =~ /#/;
        my ($chr, $start, $end, $name, $valueOrig, $strand, $info) = split("\t", $line);
	my ($shuffled) = $line =~ /SHUFFLED=(\d+)/;

        # Next if DRIPc value too big (more than 2)
       # if ($value > 2) {
       #         $randomTooBig ++;
       #         next;
       # }
       # else {
                $data{$name}{count} ++;
                my $valuez = "$chr\t$start\t$end\t$name\t$valueOrig\t$strand\t$info";
                push(@{$data{$name}{line}}, $valuez);
       # }
}
close $in5;

my $mean = 0;
open (my $out5, ">", "$dripcName.shuffled") or die "Cannot write to $dripcName.shuffled: $!\n";
foreach my $name (keys %data) {
        $lessThan50 ++ if $data{$name}{count} < 100;
        $totalPeak ++;
        if ($data{$name}{count} >= 100) {
                for (my $i = 0; $i < 100; $i++) {
                        print $out5 "$data{$name}{line}[$i]\t$i\n";
                }
        }
        $mean += $data{$name}{count};
}
close $out5;
my $remainingPeak = $totalPeak - $lessThan50;
printf "MeanPeak = $mean / $remainingPeak %.2f\n", $mean / $totalPeak;
print "Peak Used = $remainingPeak\nRandomTooBig = $randomTooBig\nlessThan50 = $lessThan50\nTotal Peak = $totalPeak\n";
# e. sort
#mitochy::runbash("cat $dripcName.shuffled | sort -k1,1 -k2,2n > $dripcName.tmp && mv $dripcName.tmp $dripcName.shuffled");
mitochy::runbash("cat $dripcName.shuffled | sort -k4,4 -k8,8n > $dripcName.tmp && mv $dripcName.tmp $dripcName.shuffled");

