#!/usr/bin/perl
# Version: 15 Oct 2014
# - Parse Original
# - Parse Shuffled
# - Get the first 150 shuffled (otherwise discard)
use strict; use warnings; use mitochy; use R_toolbox; use Statistics::Test::WilcoxonRankSum; use Math::CDF qw (:all); use Statistics::Multtest qw(:all);
my ($origInput) = @ARGV;
die "usage: $0 <Histone_origInput.txt>\n" unless @ARGV == 1;

my ($folder, $fileName) = mitochy::getFilename($origInput, "folder");
my $shuffleInput = "MappedShuffled/$fileName.txt"; die "Shuffled $shuffleInput does not exist\n" if not -e $shuffleInput;
my ($histone, $feature) = $fileName =~ /^(\w+)_drip\w?_(\w+)$/; die "Histone and/or Feature undefined\n" unless defined($histone) and defined($feature);
my ($rnaFile) = "/data/mitochi/Work/Project/DRIPc/data/NT2.rpkm";
my %rna = %{parse_rna($rnaFile)};
my $usd = 0;
##########################
# 1. Parse Shuffled
print "1. Processing $shuffleInput\n";
my %data;
open (my $in1, "<", $shuffleInput) or die "Cannot read from $shuffleInput: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $value, $name, $val, $strand, $info) = split("\t", $line);
	my ($gene) = $info =~ /TWIN=(\w+\.\d+),chr/;
	my ($orig) = $info =~ /ORIG=(\w+\.\d+),chr/;
	die "Undefined RNA seq for gene $gene shuffled\n" unless defined($rna{$gene});
	my $rna = $rna{$orig};
	#next if $rna <50 or $rna > 150;
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
my (@origR, @randomR, @fold);
my @random;
my ($nge, $nle, $nrun) = (0,0,0);
my ($sigUp, $sigDown) = (0,0);
foreach my $name (keys %data) {
	my $valueOrig = $data{$name}{orig};
	next if not defined($valueOrig);
	next if @{$data{$name}{random}} < 125;
	$usd++;

	# Pool random values for boxplot purpose.
	# If there are too many data points, then sample.
	# This is for boxplot purpose so it's ok
	# Les than 500k: Use everything
	# 500k-5 Mil: Get every 10
	# 5Mil above: Get every 100
	my $iter = 20;
	@{$data{$name}{random}} = shuffle(@{$data{$name}{random}});
	for (my $i = 0; $i < 125; $i++) {
		my $valueRandom = $data{$name}{random}[$i];
		my $fold = $valueRandom == 0 ? $valueOrig : $valueOrig / $valueRandom;
		push(@random, $valueRandom);
		push(@origR, $valueOrig);
		push(@randomR, $fold) if $i % $iter == 0;
		push(@fold, $fold);
	}
}
my $total = (keys %data);
print "USED (shuffled is 125x) = $usd / $total\n";
my $medianOrig = median(@origR);
my $meanOrig   = mean(@origR);
my $check = 0;
my $higher = 0;
my $lower  = 0;
my $notsig = 0;
my @pval;
# NOTES
# ORIGIN array is @origR
# RANDOM array is @random (ordered)
my @dataset_1 = @origR;
my @dataset_2 = @random;
my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
$wilcox_test->load_data(\@dataset_1, \@dataset_2);
my $prob = $wilcox_test->probability();
my %wilcox = %{$wilcox_test->as_hash()};
my $zval = $wilcox{probability_normal_approx}{z};
my $pval2 = $zval > 0 ? 1 - pnorm($zval) : pnorm($zval);
my $sample_size = @random;
my $pval = $pval2 * sqrt($sample_size/100);
print "pval = $pval2, sample size = $sample_size, new pval = $pval\n";
my $summary = $wilcox_test->summary();
my @summary = split("\n", $summary);
my $status = not defined($summary[10]) ? "NONE" : $summary[10] =~ /higher/i ? "HIGHER" : "LOWER";
my $fold = median(@fold);
open (my $outFold, ">>", "FOLD_ONEWILCOXON.fold") or die;
print $outFold "$feature\t$histone\t$pval\t$status\t$fold\n";
close $outFold;
mkdir "Rdata" if not -d "Rdata";
#open (my $outValueR, ">", "Rdata/$histone\_$feature\_value.tsv") or die "Cannot write to Rdata/$histone\_$feature\_value.tsv: $!\n";
#for (my $i = 0; $i < @random; $i++) {
#	print $outValueR "$origR[$i]\t$random[$i]\n";
#}
#open (my $outFoldR, ">", "Rdata/$histone\_$feature.tsv") or die "Cannot write to Rdata/$histone\_$feature.tsv: $!\n";
#for (my $i = 0; $i < @fold; $i++) {
#	print $outFoldR "$fold[$i]\n";
#}
#close $outFoldR;

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
