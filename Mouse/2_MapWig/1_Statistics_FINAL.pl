#!/usr/bin/perl
# Version: 15 Oct 2014
# - Parse Original
# - Parse Shuffled
# - Get the first 150 shuffled (otherwise discard)
use strict; use warnings; use mitochy; use R_toolbox; use Statistics::Test::WilcoxonRankSum; use Math::CDF qw (:all); use Statistics::Multtest qw(:all); use Statistics::Normality qw(:all);
my ($origInput) = @ARGV;
die "usage: $0 <Histone_origInput.txt>\n" unless @ARGV == 1;

my ($folder, $fileName) = mitochy::getFilename($origInput, "folder");
my $shufInput = "$fileName.txt"; $shufInput =~ s/orig/shuf/; die "Shuffled $shufInput does not exist\n" if not -e $shufInput;
#die "INPUT: $origInput, $shufInput\n";
my ($histone, $feature, $tag) = $fileName =~ /^(\w+)_drip\w*_(\w+)_(\w+)$/; 
($histone, $feature, $tag) = $fileName =~ /^(\w+)_enh\w*_(\w+)_(\w+)$/ if not defined($histone);
die "Histone and/or Feature and/or Tag undefined\n" unless defined($histone) and defined($feature) and defined($tag);
my ($rnaFile) = $origInput =~ /E14/ ? "/data/mitochi/Work/Project/DRIPc/data/E14.rpkm" : "/data/mitochi/Work/Project/DRIPc/data/3T3.rpkm";

my %rna = %{parse_rna($rnaFile)};
my $usd = 0;
##########################
# 1. Parse Shuffled
print "1. Processing $shufInput\n";
my %data;
my %names;
open (my $in1, "<", $shufInput) or die "Cannot read from $shufInput: $!\n";
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
	next if @{$data{$name}{random}} < 125;
	$usd++;

	@{$data{$name}{random}} = shuffle(@{$data{$name}{random}});
	for (my $i = 0; $i < 125; $i++) {
		my $valueRandom = $data{$name}{random}[$i];
		my $fold = $valueRandom == 0 ? $valueOrig : $valueOrig / $valueRandom;
		push(@random, $valueRandom);
		push(@{$origR[$i]}, $valueOrig);
		push(@{$randomR[$i]}, $valueRandom);
	}
}
die "There is no peak with more than 125 shuffles\n" if $usd == 0;

# Calculate median of Fold
for (my $i = 0; $i < 125; $i++) {
	my $medianOrig = median(@{$origR[$i]});
	my $medianShuf = median(@{$randomR[$i]});
	my $fold = $medianShuf == 0 ? 1+$medianOrig : $medianOrig / $medianShuf;
	push(@foldAll, $fold);
}
my $medianAll = int(100*median(@foldAll))/100; 

# Calculate monte carlo p value
my %stat;
my @test = getTest(median(@foldAll));
my ($nge, $nle, $nrun) = (0,0,0);
for (my $j = 0; $j < @test; $j++) {
	for (my $i = 0; $i < 125; $i++) {
		my $medianOrig = median(@{$origR[$i]});
		my $medianShuf = median(@{$randomR[$i]});
		my $fold = $medianShuf == 0 ? 1+$medianOrig : $medianOrig / $medianShuf;
		$nge ++ if $medianOrig >= $medianShuf*$test[$j];
		$nle ++ if $medianOrig <= $medianShuf*$test[$j];
		$nrun ++;
	}
	my $origbigger  = int(($nge + 1) / ($nrun + 1)*1000)/1000;
	my $origsmaller = int(($nle + 1) / ($nrun + 1)*1000)/1000;
	$stat{$test[$j]} = $origsmaller if $medianAll < 1;
	$stat{$test[$j]} = $origbigger if $medianAll >= 1;
}
print "MEDIAN = $medianAll\n";
my $foldAll = join("\t", @foldAll);
open (my $outFold, ">>", "FOLD.fold") or die;
print $outFold "$feature\t$histone\t";
my $lastP = "1";
my $lastTest = 0;

my $check = 0;
foreach my $test (sort {$a <=> $b} keys %stat) {
	next if $test < 1 and $medianAll > 1;
	next if $test > 1 and $medianAll < 1;
	print "Test on $test: $stat{$test}\n";
	if ($medianAll < 1 and $stat{$test} <= 0.05) {
		print "\t$test: $test $stat{$test}\n";
		print $outFold "$test\t$stat{$test}\t";
		$check = 1;
		last;
	}
	elsif ($medianAll >= 1 and $stat{$test} <= 0.05) {
		$lastP = $stat{$test};
		$lastTest = $test;
	}
}
if ($check == 0) {
	print $outFold "1\t1\t" if $medianAll < 1;
	print $outFold "1\t1\t" if $medianAll >= 1 and $lastP >= 0.05;
	print $outFold "$lastTest\t$lastP\t" if $medianAll >= 1 and $lastP < 0.05;
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

sub getTest {
	my ($median) = @_;
	my @test;
	if ($median < 1) {
		# Test with 0.001, 0.01, 0.1 + median
		push(@test, $median + 0.001, $median + 0.01, $median + 0.1, 1);
	}
	else {
		# Test with 0.001, 0.01, 0.1 + median
		push(@test, 1, $median - 0.001, $median - 0.01, $median - 0.1);
	}
	return(@test);
}


__END__
