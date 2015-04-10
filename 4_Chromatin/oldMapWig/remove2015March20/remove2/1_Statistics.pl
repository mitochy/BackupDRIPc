#!/usr/bin/perl
# Version: 15 Oct 2014
# - Parse Original
# - Parse Shuffled
# - Get the first 150 shuffled (otherwise discard)
use strict; use warnings; use mitochy; use R_toolbox;

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
	my $rna = $rna{$gene};
	next if $rna < 50;
	my @arr = split("\t", $line);
	#push(@{$data{$name}{random}}, $value);
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
printf "DRIPc peak not used: $undefined / $totalDRIPc (%.2f %%)\n", 100 * $undefined / $totalDRIPc;
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
	push(@origR, $valueOrig);

	# Pool random values for boxplot purpose.
	# If there are too many data points, then sample.
	# This is for boxplot purpose so it's ok
	# Les than 500k: Use everything
	# 500k-5 Mil: Get every 10
	# 5Mil above: Get every 100
	my $totalRandom = @{$data{$name}{random}};
	my $iter = 5;#$totalRandom <= 50 ? 1 : $totalRandom < 500 ? 10 : 100;
	@{$data{$name}{random}} = shuffle(@{$data{$name}{random}});
	for (my $i = 0; $i < @{$data{$name}{random}}; $i++) {
		my $valueRandom = $data{$name}{random}[$i];
		push(@randomR, $valueRandom) if $i % $iter == 0;
		push(@{$random[$i]}, $valueRandom);
	}
}
my $total = (keys %data);
print "USED = $usd / $total\n";
my $medianOrig = median(@origR);
my $meanOrig   = mean(@origR);
my $check = 0;
open (my $outS, ">", "TEST.tsv") or die;
for (my $i = 0; $i < 124; $i++) {#@random; $i++) {

	my $medianRandom = median(@{$random[$i]});
	my $meanRandom   = mean(@{$random[$i]});
	if ($medianRandom == 0) {
		$medianRandom = $meanRandom;
		$medianOrig   = $meanOrig;
		$check = 1;
	}
	#next if $medianRandom == 0;
	#my $folds = $medianOrig < 0 and $medianRandom < 0 ? -
	print $outS "$origR[$i]";
	for (my $j = 0; $j < 5; $j++) {
		print $outS "\t$random[$i][$j]";
	}
	print $outS "\n";
	my $foldVal = $medianOrig / $medianRandom;
	#$foldVal .= "ORIG:$medianOrig,RAND:$medianRandom";
	push(@fold, $foldVal);
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
if ($check == 0) {
	print $outFold "$feature\t$histone\t$pvalUp\tHigher\t$fold\n" if $pvalUp <= 0.05;
	print $outFold "$feature\t$histone\t$pvalDown\tLower\t$fold\n" if $pvalDown <= 0.05;
	print $outFold "$feature\t$histone\t1\tNone\t$fold\n" if $pvalUp > 0.05 and $pvalDown > 0.05;
}
else {
	print $outFold "$feature\t$histone\t$pvalUp\tmean_Higher\t$fold\n" if $pvalUp <= 0.05;
	print $outFold "$feature\t$histone\t$pvalDown\tmean_Lower\t$fold\n" if $pvalDown <= 0.05;
	print $outFold "$feature\t$histone\t1\tmean_None\t$fold\n" if $pvalUp > 0.05 and $pvalDown > 0.05;

}
close $outFold;

my $orig = R_toolbox::newRArray(\@origR, "orig");
my $random = R_toolbox::newRArray(\@randomR, "random");
my $foldR = R_toolbox::newRArray(\@fold, "fold");
my $Rscript = "
$foldR
#$orig
#$random
#df = rbind(df,data.frame(histone=\"$histone\",type=\"Original\",val=orig))
#df = rbind(df,data.frame(histone=\"$histone\",type=\"Random\",val=random))
df = rbind(df,data.frame(histone=\"$histone\",val=fold))
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

sub shuffle {
	my (@value) = @_;
	for (my $i = 0; $i < @value; $i++) {
		my $rand1 = int(rand(@value));
		my $rand2 = int(rand(@value));
		my $val1 = $value[$rand1];
		my $val2 = $value[$rand2];
		$value[$rand1] = $val2;
		$value[$rand2] = $val1;
	}
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
