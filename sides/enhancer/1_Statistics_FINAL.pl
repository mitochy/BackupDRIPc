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
my ($histone, $feature) = $fileName =~ /^(\w+)_enh\w?_(\w+)$/; die "Histone and/or Feature undefined\n" unless defined($histone) and defined($feature);
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
	my ($gene) = $info =~ /TWIN=(\w+\.\d+)/;
	my ($orig) = $info =~ /ORIG=(\w+\.\d+)/;
	$names{$name}{orig} = $orig;
	$names{$name}{gene} = $gene;
	die "Undefined RNA seq for gene $gene shuffled\n" unless defined($rna{$gene});
	my $rna = $rna{$orig};
	#next if $rna < 10;
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
	# If the name isn't the same as orig
	my $orig = $names{$name}{orig};
	print "Undefined orig for name $name\n" and next if (not defined($orig));
	$name = $orig;
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
		#push(@origR, $valueOrig);
		#push(@randomR, $fold) if $i % $iter == 0;
		push(@{$origR[$i]}, $valueOrig);
		push(@{$origR2[$i]}, $valueRandom);
		push(@{$randomR[$i]}, $valueRandom);
		push(@{$randomR2[$i]}, $valueOrig);
		#push(@{$fold[$i]}, $fold);
	}
}
die "There is no peak with more than 125 shuffles\n" if $usd == 0;

# Calculate monte carlo p value
my %stat;
my @test = (1,1.15,1.3,1.6);
my $pval2;
for (my $j = 0; $j < @test; $j++) {
	my ($nge1, $nle1, $nge2, $nle2, $nrun) = (0,0,0,0,0);
	for (my $i = 0; $i < 125; $i++) {
		my $medianOrig = median(@{$origR[$i]});
		my $medianShuf = median(@{$randomR[$i]});
		push(@foldAll, ($medianOrig) / ($medianShuf)) if $j == 0;
		$nge1 ++ if $medianOrig <= $medianShuf*$test[$j];
		#$nle1 ++ if $medianOrig >= $medianShuf*$test[$j];
		$nge2 ++ if $medianOrig <= $medianShuf*1/$test[$j];
		#$nle2 ++ if $medianOrig >= $medianShuf*1/$test[$j];

		$nrun ++;
	}
	#my $pval1Up = int(($nge1 + 1) / ($nrun + 1)*1000)/1000;
	my $origbigger = int(($nge1 + 1) / ($nrun + 1)*1000)/1000;
	#my $pval2Up = int(($nge2 + 1) / ($nrun + 1)*1000)/1000;
	my $origsmaller = int(($nge2 + 1) / ($nrun + 1)*1000)/1000;
	#$pval2 .= "\t$pval2Up\t$pval2Down";
	my $test2 = int(1 / $test[$j]*100)/100;
	#$pval2 .= "FOLD ORIG / SHUF = $test[$j]=$origbigger\t$test2=$origsmaller\n";#$pval1Up\t$pval1Down\t$pval2Up\t$pval2Down\n";
	$pval2 .= "\t>$test[$j]=$origbigger\t>$test2=$origsmaller";#$pval1Up\t$pval1Down\t$pval2Up\t$pval2Down\n";
}
my $medianAll = median(@foldAll);
#print "\t\tOrig > Shuf * 1.5\t\tOrig < Shuf * 0.67\n";
#die "MEDIAN ORIG / SHUF = $medianAll\n$pval2\n";
open (my $outFold, ">>", "FOLD2.fold") or die;
#print $outFold "$feature\t$histone\t$pvalUp\tHigher\t$fold\n" if $pvalUp <= 0.05;
#print $outFold "$feature\t$histone\t$pvalDown\tLower\t$fold\n" if $pvalDown <= 0.05;
#print $outFold "$feature\t$histone\t1\tNone\t$fold\n" if $pvalUp > 0.05 and $pvalDown > 0.05;
print $outFold "$feature\t$histone$pval2\t$medianAll\n";# if $pval2Up <= 0.05;
#print $outFold "$feature\t$histone$pval2\tLower\t$medianAll\n" if $pval2Down <= 0.05;
#print $outFold "$feature\t$histone$pval2\tNone\t$medianAll\n" if $pval2Up > 0.05 and $pval2Down > 0.05;
close $outFold;
die;
my $total = (keys %data);
print "USED (shuffled is 125x) = $usd / $total\n";
# ORIGIN array is @origR
# RANDOM array is @random (ordered)
mkdir "Rdata" if not -d "Rdata";
#if (not -e "Rdata/$histone\_$feature\_orig.tsv") {
	open (my $outData, ">", "Rdata/$histone\_$feature\_data.tsv") or die "Cannot write to Rdata/$histone\_$feature\_data.tsv: $!\n";
	for (my $i = 0; $i < @origR; $i++) {
		for (my $j = 0; $j < @{$origR[$i]}; $j++) {
			print $outData "$origR[$i][$j]\t$randomR[$i][$j]\n";
		}
	}
	close $outData;
#if (not -e "Rdata/$histone\_$feature\_fold.tsv") {
	open (my $outFoldR, ">", "Rdata/$histone\_$feature\_fold.tsv") or die "Cannot write to Rdata/$histone\_$feature\_fold.tsv: $!\n";
	for (my $i = 0; $i < @foldAll; $i++) {
		print $outFoldR "$foldAll[$i]\n";
	}
	close $outFoldR;
#}

my $Rscript = "
library(ggplot2)

df = read.table(\"Rdata/$histone\_$feature\_data.tsv\")
orig = df\$V1
shuf = df\$V2
#shuf = read.table(\"Rdata/$histone\_$feature\_shuf.tsv\")\$V1;
# any diff
mu = 0
origbigger0 = wilcox.test(orig,shuf,paired=FALSE,mu=median(shuf)*mu,alternative=\"greater\")\$p.value
origsmaller0 = wilcox.test(orig,shuf,paired=FALSE,mu=median(shuf)*1,alternative=\"greater\")\$p.value
# 1.15x fold enrichment or depletion (1/1.15)
mu = 1.15
origbigger1 = wilcox.test(orig,shuf,paired=FALSE,mu=median(shuf)*(mu-1),alternative=\"greater\")\$p.value
origsmaller1 = wilcox.test(orig,shuf,paired=FALSE,mu=median(shuf)*(1/mu-1),alternative=\"greater\")\$p.value
mu = 1.3
origbigger2 = wilcox.test(orig,shuf,paired=FALSE,mu=median(shuf)*(mu-1),alternative=\"greater\")\$p.value
origsmaller2 = wilcox.test(orig,shuf,paired=FALSE,mu=median(shuf)*(1/mu-1),alternative=\"greater\")\$p.value
mu = 1.6
origbigger3 = wilcox.test(orig,shuf,paired=FALSE,mu=median(shuf)*(mu-1),alternative=\"greater\")\$p.value
origsmaller3 = wilcox.test(orig,shuf,paired=FALSE,mu=median(shuf)*(1/mu-1),alternative=\"greater\")\$p.value
pval = paste(origbigger0, origsmaller0, origbigger1,origsmaller1,origbigger2,origsmaller2,origbigger3,origsmaller3,sep=\"\t\")
write.table(pval,file=\"Rdata/$histone\_$feature\_pval.txt\");
#fold = read.table(\"Rdata/$histone\_$feature\_fold.tsv\");
";

R_toolbox::execute_Rscript($Rscript);
my (@pval) = `cat Rdata/$histone\_$feature\_pval.txt`;
my ($pvals) = $pval[1] =~ / "(.+)"$/;
open (my $outz, ">>", "FOLD_ONEWILCOXON_R.txt") or die;
printf $outz "$histone\t$feature";
@pval = split("\t", $pvals);
for (my $i = 0; $i < @pval; $i++) {
	if ($pval[$i] > 0.001) {
		$pval[$i] = int($pval[$i] * 1000) / 1000;
		print $outz "\t$pval[$i]";
	}
	else {
		printf $outz "\t%.2e", $pval[$i];
	}
}
print $outz "\t$medianAll\n";
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
