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
my $shuffleNum = 125;
my $iter = 10;
my @foldAll;
for (my $h = 0; $h <= 1000; $h+=$iter) {
my $usd = 0;
$iter = 50 if $h >= 150;
$iter = 200 if $h >= 400;
print "Processing $h\n";
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
	if ($h < 1000) {
		next if $rna >= $h+$iter or $rna <$h;
	}
	else {
		next if $rna < 1000;
	}
	$shuffleNum = 50 if $h < 10;
	$shuffleNum = 125 if $h >= 10;
	my $rna2 = $rna{$gene};
	next if ($rna2 < $rna *0.9 or $rna2 > $rna *1.1);
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
	next if @{$data{$name}{random}} < $shuffleNum;
	$usd++;

	# Pool random values for boxplot purpose.
	# If there are too many data points, then sample.
	# This is for boxplot purpose so it's ok
	# Les than 500k: Use everything
	# 500k-5 Mil: Get every 10
	# 5Mil above: Get every 100
	my $totalRandom = @{$data{$name}{random}};
	my $iter = 5;#$totalRandom <= 50 ? 1 : $totalRandom < 500 ? 10 : 100;
	@{$data{$name}{random}} = shuffle(@{$data{$name}{random}});
	#($nge, $nle, $nrun) = (0,0,0);
	#for (my $i = 0; $i < @{$data{$name}{random}}; $i++) {
	for (my $i = 0; $i < $shuffleNum; $i++) {
		my $valueRandom = $data{$name}{random}[$i];
		push(@{$random[$i]}, $valueRandom);
		my $fold = $valueRandom == 0 ? $valueOrig : $valueOrig / $valueRandom;
		#push(@random, $valueRandom);
		push(@origR, $valueOrig);
		push(@randomR, $fold) if $i % $iter == 0;
		push(@fold, $fold);
		#$nge ++ if $valueRandom >= $valueOrig;
		#$nle ++ if $valueRandom <= $valueOrig;
		#$nrun++;
	}
	#my $pvalUp = int(100*($nge + 1) / ($nrun + 1))/100;
	#my $pvalDown = int(100*($nle + 1) / ($nrun + 1))/100;
	#$sigUp ++ if $pvalUp <= 0.05;
	#$sigDown ++ if $pvalDown <= 0.05;
	#print "$name\tpUp $pvalUp\tpDown $pvalDown\tNGE $nge\tNLE $nle\t$rna{$name}\n";
}
#print "Sigup $sigUp Down $sigDown\n";
my $total = (keys %data);
print "USED (shuffled is $shuffleNum) = $usd / $total\n";
#die;
my $medianOrig = median(@origR);
my $meanOrig   = mean(@origR);
my $check = 0;
#open (my $outS, ">", "TEST.tsv") or die;
my $higher = 0;
my $lower  = 0;
my $notsig = 0;
my @pval;
push(@foldAll, @fold);
for (my $i = 0; $i < $shuffleNum; $i++) {#@random; $i++) {

	# NOTES
	# ORIGIN array is @origR
	# RANDOM array is @{$random[$i]}
	
	my @dataset_1 = @origR;
	my @dataset_2 = @{$random[$i]};
	#my @dataset_2 = @random;
	my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
	$wilcox_test->load_data(\@dataset_1, \@dataset_2);
	my $prob = $wilcox_test->probability();
	#my $status = $wilcox_test->probability_status();
	#my $smaller = $wilcox_test->get_smaller_rank_sum();
	#my $summary = $wilcox_test->summary();
	#my @summary = split("\n", $summary);
	#print "$summary[8]\t$summary[9]\t$summary[10]\n";
	my %wilcox = %{$wilcox_test->as_hash()};
	my $zval = $wilcox{probability_normal_approx}{z};
	my $pval = $zval > 0 ? 1 - pnorm($zval) : pnorm($zval);
	
	#print "PVAL $pval = pnorm($zval)\n";
	my $summary = $wilcox_test->summary();
	my @summary = split("\n", $summary);
	#print "$summary[8]\t$summary[9]\t$summary[10]\n";
	if ($pval <= 0.05 and defined($summary[10])) {
		$higher ++ if $summary[10] =~ /higher/i;
		$lower ++ if $summary[10] =~ /lower/i;
	}
	else {
		$notsig ++;
	}
	push(@pval, $pval);
	#print "ZVAL $zval PVAL $pval $summary[10]\n" if $pval <= 0.05 and defined($summary[10]);
	#print "ZVAL $zval PVAL $pval NS\n" if $pval <= 0.05 and not defined($summary[10]);
	#print "ZVAL $zval PVAL $pval NS\n" if $pval > 0.05;
}
my $pval = 0;
#$pval = 0;
for (my $i = 0; $i < @pval; $i++) {
	$pval ++ if $pval[$i] <= 0.05;
}

my @adjpval = @{BH(\@pval)};
my $adjpval = 0;
for (my $i = 0; $i < @adjpval; $i++) {
	$adjpval ++ if $adjpval[$i] <= 0.05;
}
my $meanadjpval = mean(@adjpval);
my $medianadjpval = median(@adjpval);
my $fold = median(@fold);
#print "\nhigher $higher lower $lower notsig $notsig FOLD $fold\n";
print "pval <= 0.05 = $pval, adjpval <= 0.05 = $adjpval mean_adjpval = $meanadjpval median_adjpval = $medianadjpval\n";

open (my $outFold, ">>", "FOLD.fold") or die;
print $outFold "$feature\t$histone\t$meanadjpval\t$h\t$higher,$lower,$notsig\t$undefined\t$totalDRIPc\t$usd\t$total\t$fold\n" if $meanadjpval <= 0.05;
print $outFold "$feature\t$histone\t$meanadjpval\t$h\tNS\t$undefined\t$totalDRIPc\t$usd\t$total\t$fold\n" if $meanadjpval > 0.05;
}

open (my $outFoldAll, ">", "Rdata/$feature\_$histone.fold") or die;
foreach my $foldAll (@foldAll) {
	print $outFoldAll "$foldAll\n";
}
=comment
my $orig = R_toolbox::newRArray(\@origR, "orig");
my $random = R_toolbox::newRArray(\@randomR, "random");
#my $foldR = R_toolbox::newRArray(\@fold, "fold");
my $Rscript = "
#$orig
#$random
#df = rbind(df,data.frame(histone=\"$histone\",type=\"Original\",val=orig))
#df = rbind(df,data.frame(histone=\"$histone\",type=\"Random\",val=random))
df = rbind(df,data.frame(histone=\"$histone\",val=fold))
";

open (my $ou, ">>", "$feature.R") or die;
print $ou "$Rscript\n";
close $ou;
=cut
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
