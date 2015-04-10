#!/usr/bin/perl
# Version: 15 Oct 2014
# - Parse Original
# - Parse Shuffled
# - Get the first 150 shuffled (otherwise discard)
use strict; use warnings; use mitochy; use R_toolbox; use Statistics::Test::WilcoxonRankSum; use Math::CDF qw (:all); use Statistics::Multtest qw(:all); use Statistics::Normality qw(:all);
my ($origInput) = @ARGV;
die "usage: $0 <Histone_origInput.txt>\n" unless @ARGV == 1;

my ($folder, $fileName) = mitochy::getFilename($origInput, "folder");
my $shuffleInput = "MappedShuffled/$fileName.txt"; die "Shuffled $shuffleInput does not exist\n" if not -e $shuffleInput;
my ($histone, $feature) = $fileName =~ /^(\w+)_drip\w?_(\w+)$/; 
($histone, $feature) = $fileName =~ /^(\w+)_enh\w*_(\w+)$/ if not defined($histone);
die "Histone and/or Feature undefined\n" unless defined($histone) and defined($feature);
my ($rnaFile) = "/data/mitochi/Work/Project/DRIPc/data/NT2.rpkm";
my %rna = %{parse_rna($rnaFile)};
my $usd = 0;
##########################
# 1. Parse Shuffled
print "1. Processing $shuffleInput\n";
my %data;
my %names;
open (my $in1, "<", $shuffleInput) or die "Cannot read from $shuffleInput: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $value, $name, $val, $strand, $info) = split("\t", $line);
	my ($gene) = $info =~ /TWIN=(\w+\.\d+),chr/;
	($gene) = $info =~ /TWIN=(\w+\.\d+)/ if not defined($gene);
	my ($orig) = $info =~ /ORIG=(\w+\.\d+),chr/;
	($orig) = $info =~ /ORIG=(\w+\.\d+)/ if not defined($orig);
	$names{$name}{orig} = $orig;
	$names{$name}{gene} = $gene;
	die "Undefined RNA seq for gene $gene shuffled\n" unless defined($rna{$gene});
	my $rna = $rna{$orig};
	next if $rna < 10;
	#next if $rna <95 or $rna > 105;
	my @arr = split("\t", $line);
	#push(@{$data{$name}{random}}, $value);
	$value = int($value * 100)/100;
	push(@{$data{$orig}{random}}, $value);
}
close $in1;
##########################

##########################
# 2. Parse Original
print "2. Processing $origInput\n";
open (my $in2, "<", $origInput) or die "Cannot read from $origInput: $!\n";
my ($totalDRIPc, $undefined) = (0,0);
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $value, $name) = split("\t", $line);
	# If the name isn't the same as orig
	$name = $names{$name}{orig} if $name !~ /ENS/;
	print "Undefined orig for name $name\n" and next if (not defined($name));
	$value = int($value * 100)/100;
	my @arr = split("\t", $line);	
	$totalDRIPc ++;
	#next if $rna{$name} >10;# and $rna{$name} > 10;
	# Next if shuffled of this data does not exist
	# Because DRIPc too low or number of twin is less than 50 or less than 100 shuffles
	if (not defined($data{$name})) {
		$undefined ++;
		next;
	}
	else {
		$data{$name}{orig} = $value;
	}
}
close $in2;
printf "DRIPc peak not used due to no shuffles: $undefined / $totalDRIPc (%.2f %%)\n", 100 * $undefined / $totalDRIPc;
##########################

##########################
# 3. 
print "3. Do Statistics!\n";
my (@origR, @randomR, @foldAll, @random);
my (@origR2, @randomR2);
my @fold;
my ($sigUp, $sigDown) = (0,0);
foreach my $name (keys %data) {
	my $valueOrig = $data{$name}{orig};
	next if not defined($valueOrig);
	next if @{$data{$name}{random}} < 25;
	$usd++;

	@{$data{$name}{random}} = shuffle(@{$data{$name}{random}});
	for (my $i = 0; $i < 25; $i++) {
		my $valueRandom = $data{$name}{random}[$i];
		my $fold = $valueRandom == 0 ? $valueOrig : $valueOrig / $valueRandom;
		push(@random, $valueRandom);
		push(@{$origR[$i]}, $valueOrig);
		push(@{$randomR[$i]}, $valueRandom);
		#push(@{$origR2[$i]}, $valueRandom);
		#push(@{$randomR2[$i]}, $valueOrig);
	}
}
die "There is no peak with more than 25 shuffles\n" if $usd == 0;

# Calculate monte carlo p value
my %stat;
my @test = (1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.25,2.5,2.75,3);
my ($nge1, $nge2, $nrun) = (0,0,0);
for (my $j = 0; $j < @test; $j++) {
	for (my $i = 0; $i < 25; $i++) {
		my $medianOrig = median(@{$origR[$i]});
		my $medianShuf = median(@{$randomR[$i]});
		my $fold = $medianShuf == 0 ? 1+$medianOrig : $medianOrig / $medianShuf;
		#print "$histone $feature shuffle is 0: $medianShuf orig = $medianOrig test $i\n";
		push(@foldAll, $fold) if $j == 0;
		$nge1 ++ if $medianOrig <= $medianShuf*$test[$j];
		$nge2 ++ if $medianOrig <= $medianShuf*1/$test[$j];
		$nrun ++;
	}
	my $origbigger = int(($nge1 + 1) / ($nrun + 1)*1000)/1000;
	my $origsmaller = 1-(int(($nge2 + 1) / ($nrun + 1)*1000)/1000);
	my $test2 = int(1 / $test[$j]*100)/100;
	$stat{$test2} = $origsmaller;
	$stat{$test[$j]} = $origbigger;
}
my $medianAll = int(100*median(@foldAll))/100; 
print "MEDIAN = $medianAll\n";
my $foldAll = join("\t", @foldAll);
open (my $outFold, ">>", "FOLD.fold") or die;
print $outFold "$feature\t$histone\t";
my $lastP = "1";
my $lastTest = 0;
foreach my $test (sort {$a <=> $b} keys %stat) {
	if ($medianAll < 1 and $stat{$test} <= 0.05 and $lastTest < 1) {
		print $outFold "$lastTest\t$lastP\t";
		last;
	}
	elsif ($medianAll < 1 and $stat{$test} > 0.05 and $lastTest >= 1) {
		print $outFold "$lastTest\t$lastP\t";
		last;
	}
	elsif ($medianAll >= 1 and $stat{$test} > 0.05 and $lastTest >= 1) {
		print $outFold "$lastTest\t$lastP\t";
		last;
	}
	elsif ($medianAll >= 1 and $stat{$test} > 0.05 and $lastTest >= 1) {
		print $outFold "$lastTest\t$lastP\t";
		last;
	}
	else {
		$lastP = $stat{$test};
		$lastTest = $test;
	}
}
print $outFold "$medianAll\t$foldAll\n";
close $outFold;

sub mean {
	my (@value) = @_;
	my $mean;
	foreach my $value (@value) {
		print "MEAN: UNDEFINED VALUE\n" and next if not defined($value);
		if ($value eq "NA" or $value =~ /inf/i) {
			$value = 0;
			print "NON DIGIT VALUE FOUND AT $value\n";
		}
		$mean += $value / @value;
	}
	return($mean);
}

sub median {
	my (@value) = @_;
	@value = sort {$a <=> $b} (@value);
	my $median = $value[int(@value/2)];
	return($median);
}

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

sub parse_rna {
        my ($input) = @_;
        my %data;
        open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
        while (my $line = <$in>) {
                chomp($line);
                next if $line =~ /#/;
                my ($gene, $val) = split("\t", $line);
                $data{$gene} = $val;
        }
        close $in;
        return(\%data);
}
__END__
