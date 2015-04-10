#!/usr/bin/perl
# - Parse Original
# - Parse Shuffled
# - Get the first 50 shuffled (otherwise discard)
use strict; use warnings; use mitochy; use R_toolbox;

my ($origInput) = @ARGV;
die "usage: $0 <Histone_origInput.txt>\n" unless @ARGV == 1;

my ($folder, $fileName) = mitochy::getFilename($origInput, "folder");
my $shuffleInput = "MappedShuffled/$fileName.txt"; die "Shuffled $shuffleInput does not exist\n" if not -e $shuffleInput;
my ($histone, $feature) = $fileName =~ /^(\w+)_drip_(\w+)$/; die "Histone and/or Feature undefined\n" unless defined($histone) and defined($feature);

##########################
# 1. Parse Shuffled
my %data;
open (my $in1, "<", $shuffleInput) or die "Cannot read from $shuffleInput: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $value, $name) = split("\t", $line);
	my @arr = split("\t", $line);
	push(@{$data{$name}{random}}, $value);
}
close $in1;
##########################

##########################
# 2. Parse Original
open (my $in2, "<", $origInput) or die "Cannot read from $origInput: $!\n";
my ($totalDRIPc, $undefined) = (0,0);
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $value, $name) = split("\t", $line);
	my @arr = split("\t", $line);	
	$totalDRIPc ++;
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
#printf "DRIPc peak not used: $undefined / $totalDRIPc (%.2f %%)\n", 100 * $undefined / $totalDRIPc;
##########################

##########################
# 3. 
my (@origR, @randomR, @fold);
my @random;
my ($nge, $nle, $nrun) = (0,0,0);
my ($sigUp, $sigDown) = (0,0);
foreach my $name (keys %data) {
	my $valueOrig = $data{$name}{orig};
	push(@origR, $valueOrig);

	# Pool random values for boxplot purpose.
	# If there are too many data points, then sample.
	# This is for boxplot purpose so it's ok
	# Les than 500k: Use everything
	# 500k-5 Mil: Get every 10
	# 5Mil above: Get every 100
	my $totalRandom = @{$data{$name}{random}};
	my $iter = $totalRandom <= 500000 ? 1 : $totalRandom < 5000000 ? 10 : 100;
	for (my $i = 0; $i < @{$data{$name}{random}}; $i++) {
		my $valueRandom = $data{$name}{random}[$i];
		push(@randomR, $valueRandom) if $i % $iter == 0;
		push(@{$random[$i]}, $valueRandom);
	}
}

my $medianOrig = median(@origR);
for (my $i = 0; $i < @random; $i++) {
	my $medianRandom = median(@{$random[$i]});
	my $meanRandom   = mean(@{$random[$i]});
	#$medianRandom = $meanRandom if $medianRandom == 0;
	next if $medianRandom == 0;
	push(@fold, $medianOrig / $medianRandom);
	$nge++ if $medianOrig <= $medianRandom;
	$nle++ if $medianOrig >= $medianRandom;
	#print "nge = $nge, orgi = $medianOrig, rand = $medianRandom\n" if $medianOrig <= $medianRandom;
	#print "nle = $nle, orgi = $medianOrig, rand = $medianRandom\n" if $medianOrig >= $medianRandom;
	$nrun++;
}
my $fold = join("\t", @fold);
my $pvalUp = ($nge + 1) / ($nrun + 1);
my $pvalDown = ($nle + 1) / ($nrun + 1);

open (my $outFold, ">>", "FOLD.fold") or die;
print $outFold "$feature\t$histone\t$pvalUp\tHigher\t$fold\n" if $pvalUp <= 0.05;
print $outFold "$feature\t$histone\t$pvalDown\tLower\t$fold\n" if $pvalDown <= 0.05;
print $outFold "$feature\t$histone\t1\tNone\t$fold\n" if $pvalUp > 0.05 and $pvalDown > 0.05;
close $outFold;

my $orig = R_toolbox::newRArray(\@origR, "orig");
my $random = R_toolbox::newRArray(\@randomR, "random");

my $Rscript = "
$orig
$random
df = rbind(df,data.frame(histone=\"$histone\",type=\"Original\",val=orig))
df = rbind(df,data.frame(histone=\"$histone\",type=\"Random\",val=random))
";

open (my $ou, ">>", "$feature.R") or die;
print $ou "$Rscript\n";
close $ou;

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
__END__
die;

