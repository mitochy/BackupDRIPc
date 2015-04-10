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
	#$value = int($value * 100)/100;
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
my (@origR, @randomR, @foldAll);
my @fold;
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
	push(@origR, $valueOrig);
	for (my $i = 0; $i < 125; $i++) {
		my $valueRandom = $data{$name}{random}[$i];
		my $fold = $valueRandom == 0 ? $valueOrig : $valueOrig / $valueRandom;
		push(@randomR, $valueRandom);
		#push(@origR, $valueOrig);
		#push(@randomR, $fold) if $i % $iter == 0;
		push(@{$fold[$i]}, $valueRandom);
	}
}
for (my $i = 0; $i < 125; $i++) {
	my $medianOrig = median(@origR);
	my $medianShuf = median(@{$fold[$i]});
	push(@foldAll, ($medianOrig+1) / ($medianShuf+1));
}
my $total = (keys %data);
print "USED (shuffled is 125x) = $usd / $total\n";
my $medianAll = median(@foldAll);
# ORIGIN array is @origR
# RANDOM array is @random (ordered)
mkdir "Rdata" if not -d "Rdata";
#if (not -e "Rdata/$histone\_$feature\_orig.tsv") {
	open (my $outOrigR, ">", "Rdata/$histone\_$feature\_orig.tsv") or die "Cannot write to Rdata/$histone\_$feature\_orig.tsv: $!\n";
	for (my $i = 0; $i < @origR; $i++) {
		print $outOrigR "$origR[$i]\n";
	}
	close $outOrigR;
#}
#if (not -e "Rdata/$histone\_$feature\_shuf.tsv") {
	open (my $outShufR, ">", "Rdata/$histone\_$feature\_shuf.tsv") or die "Cannot write to Rdata/$histone\_$feature\_shuf.tsv: $!\n";
	for (my $i = 0; $i < @randomR; $i++) {
		print $outShufR "$randomR[$i]\n";
	}
	close $outShufR;
#}
#if (not -e "Rdata/$histone\_$feature\_fold.tsv") {
	open (my $outFoldR, ">", "Rdata/$histone\_$feature\_fold.tsv") or die "Cannot write to Rdata/$histone\_$feature\_fold.tsv: $!\n";
	for (my $i = 0; $i < @foldAll; $i++) {
		print $outFoldR "$foldAll[$i]\n";
	}
	close $outFoldR;
#}

my $Rscript = "
library(ggplot2)

orig = read.table(\"Rdata/$histone\_$feature\_orig.tsv\")\$V1;
shuf = read.table(\"Rdata/$histone\_$feature\_shuf.tsv\")\$V1;
mu = median(shuf)
origbigger0 = wilcox.test(shuf,orig,paired=FALSE,mu=-mu,alternative=\"less\")\$p.value
origsmaller0 = wilcox.test(shuf,orig,paired=FALSE,mu=mu,alternative=\"greater\")\$p.value
mu = median(shuf) * 0.15
origbigger1 = wilcox.test(shuf,orig,paired=FALSE,mu=-mu,alternative=\"less\")\$p.value
origsmaller1 = wilcox.test(shuf,orig,paired=FALSE,mu=mu,alternative=\"greater\")\$p.value
mu = median(shuf) * 0.3
origbigger2 = wilcox.test(shuf,orig,paired=FALSE,mu=-mu,alternative=\"less\")\$p.value
origsmaller2 = wilcox.test(shuf,orig,paired=FALSE,mu=mu,alternative=\"greater\")\$p.value
mu = median(shuf) * 0.6
origbigger3 = wilcox.test(shuf,orig,paired=FALSE,mu=-mu,alternative=\"less\")\$p.value
origsmaller3 = wilcox.test(shuf,orig,paired=FALSE,mu=mu,alternative=\"greater\")\$p.value
pval = paste(origbigger0, origsmalle0, origbigger1,origsmaller1,origbigger2,origsmaller2,origbigger3,origsmaller3,sep=\"\t\")
write.table(pval,file=\"Rdata/$histone\_$feature\_pval.txt\");
#fold = read.table(\"Rdata/$histone\_$feature\_fold.tsv\");
";

R_toolbox::execute_Rscript($Rscript);
my (@pval) = `cat Rdata/$histone\_$feature\_pval.txt`;
my ($pvalUp, $pvalDown) = $pval[1] =~ / "(.+)\t(.+)"$/;
open (my $outz, ">>", "FOLD_ONEWILCOXON_R.txt") or die;
print $outz "$histone\t$feature\t$pvalUp\t$pvalDown\t$medianAll\n";
close $outz;

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
