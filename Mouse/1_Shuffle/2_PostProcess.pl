#!/usr/bin/perl

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
        my ($chr, $start, $end, $name, $valueOrig, $strand, $info) = split("\t", $line);
	my ($DRIP) = $info =~ /DRIP=(.+);ORIG/;
	my ($orig) = $info =~ /ORIG=(E\w+\.\d+),chr/; die if not defined($orig);
	$totalPeak{$DRIP}++;
        $data{$DRIP}{count} ++;
	$data{$DRIP}{value} = $valueOrig;
        my $valuez = "$chr\t$start\t$end\t$name\t$valueOrig\t$strand\t$info";
        push(@{$data{$DRIP}{line}}, $valuez);
	push(@{$data{$DRIP}{info}}, $info);
}
close $in5;
my $keys = keys %totalPeak;
print "Done, printing out. Total peak = $keys\n";
# Get dripc peaks that has $desiredPeak number of shuffles
my $mean = 0;
my %good;
my $dripcTSSname = $dripcName;
$dripcTSSname =~ s/dripc_/dripcTSS_/ if $dripcTSSname =~ /dripc/;
$dripcTSSname =~ s/drip/dripTSS/ if $dripcTSSname !~ /dripc/;
die if not defined($dripcTSSname);
open (my $out1, ">", "$dripcName\_orig.peak") or die "Cannot write to $dripcName\_orig.peak: $!\n";
open (my $out2, ">", "$dripcName\_shuf.peak") or die "Cannot write to $dripcName\_shuf.peak: $!\n";
open (my $out3, ">", "$dripcName\_orig.metaplot") or die "Cannot write to $dripcName\_orig.metaplot: $!\n";
open (my $out4, ">", "$dripcName\_shuf.metaplot") or die "Cannot write to $dripcName\_shuf.metaplot: $!\n";
open (my $out5, ">", "$dripcTSSname\_orig.metaplot") or die "Cannot write to $dripcTSSname\_orig.metaplot: $!\n";
open (my $out6, ">", "$dripcTSSname\_shuf.metaplot") or die "Cannot write to $dripcTSSname\_shuf.metaplot: $!\n";
foreach my $name (keys %data) {
        $LessThanDesiredPeakNumberd ++ if $data{$name}{count} < $desiredPeak;
        $mean += $data{$name}{count};
        $totalPeak_withLowDRIPcValue ++;
        if ($data{$name}{count} >= $desiredPeak) {
		my $info = $data{$name}{info}[0];
		my ($dripchr, $dripstart, $dripend, $dripstrand, $orig, $origchr, $origstart, $origend, $shuf, $shufchr, $shufstart, $shufend) = $info =~ /DRIP=(chr\w+),(\d+),(\d+),([\-\+]);ORIG=(ENS\w+\.\d+),(chr\w+),(\d+),(\d+);TWIN=(ENS\w+\.\d+),(chr\w+),(\d+),(\d+)/;
		die "Died at $name\n" if not defined($origend) or not defined($shufend);

		# Print ORIG PEAK
		print $out1 "$dripchr\t$dripstart\t$dripend\t$orig\t$data{$name}{value}\t$dripstrand\t$info\n";
		# Print ORIG METAPLOT (get center then +/- 5k)
		my $newdripstart = int(($dripstart + $dripend)/2) - 5000;
		my $newdripend   = int(($dripstart + $dripend)/2) + 5000;
		print $out3 "$dripchr\t$newdripstart\t$newdripend\t$orig\t$data{$name}{value}\t$dripstrand\t$info\n";
		# PRINT ORIG METAPLOT AT TSS (get origstart and origend then +/- 3000
		my $neworigstart = $origstart - 3000;
		my $neworigend   = $origend   + 3000;
		print $out5 "$origchr\t$neworigstart\t$neworigend\t$orig\t$data{$name}{value}\t$dripstrand\t$info\n";

		my @random = shuffle(@{$data{$name}{line}});
                for (my $i = 0; $i < $desiredPeak; $i++) {
			# PRINT SHUF PEAK
			print $out2 "$random[$i]\n";
			# Print SHUF METAPLOT (get center then +/- 5k)
			my ($chr, $start, $end, $name, $valueOrig, $strand, $info2) = split("\t", $random[$i]);
			($dripchr, $dripstart, $dripend, $dripstrand, $orig, $origchr, $origstart, $origend, $shuf, $shufchr, $shufstart, $shufend) = $info2 =~ /DRIP=(chr\w+),(\d+),(\d+),([\-\+]);ORIG=(ENS\w+\.\d+),(chr\w+),(\d+),(\d+);TWIN=(ENS\w+\.\d+),(chr\w+),(\d+),(\d+)/;
			my $newstart = int(($start + $end)/2) - 5000;
			my $newend   = int(($start + $end)/2) + 5000;
			print $out4 "$chr\t$newstart\t$newend\t$shuf\t$valueOrig\t$strand\t$info2\n";
			# PRINT ORIG METAPLOT AT TSS (get origstart and origend then +/- 3000
			my $newshufstart = $shufstart - 3000;
			my $newshufend   = $shufend   + 3000;
			print $out6 "$shufchr\t$newshufstart\t$newshufend\t$shuf\t$valueOrig\t$strand\t$info2\n";
                }
        }
}

print "Done, sorting\n";
system("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcName\_orig.peak > $dripcName\_orig.peak.temp && mv $dripcName\_orig.peak.temp $dripcName\_orig.peak");
system("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcName\_shuf.peak > $dripcName\_shuf.peak.temp && mv $dripcName\_shuf.peak.temp $dripcName\_shuf.peak");
system("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcName\_orig.metaplot > $dripcName\_orig.metaplot.temp && mv $dripcName\_orig.metaplot.temp $dripcName\_orig.metaplot");
system("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcName\_shuf.metaplot > $dripcName\_shuf.metaplot.temp && mv $dripcName\_shuf.metaplot.temp $dripcName\_shuf.metaplot");
system("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcTSSname\_orig.metaplot > $dripcTSSname\_orig.metaplot.temp && mv $dripcTSSname\_orig.metaplot.temp $dripcTSSname\_orig.metaplot");
system("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcTSSname\_shuf.metaplot > $dripcTSSname\_shuf.metaplot.temp && mv $dripcTSSname\_shuf.metaplot.temp $dripcTSSname\_shuf.metaplot");
my $totalPeakTotal = (keys %totalPeak);
print "
Total DRIPc peak in this file = $totalPeakTotal
(DEPRECATED) Total DRIPc peak after filter for high DRIPc value = $totalPeak_withLowDRIPcValue
Desired Number of Peak = $desiredPeak
Not used = $LessThanDesiredPeakNumberd
";


sub shuffle {
        my (@value) = @_;
        #print "Before: @value\n";
        for (my $i = 0; $i < 10000; $i++) {
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
