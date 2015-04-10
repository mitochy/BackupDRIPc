#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox; use Statistics::Basic qw(:all);

my ($input) = @ARGV;
die "$0 <input>\n" unless @ARGV == 1;
my ($folder, $fileName) = mitochy::getFilename($input, "folder");

my $NT2_DRIP = `bedtools intersect -u -a $input -b NT2_DRIP.bed`;
my ($NT2_DRIPOrig)  = `bedtools intersect -u -a $fileName.name -b NT2_DRIP.bed | wc -l` =~ /^(\d+)/;
my ($NT2_DRIPTotal) = `wc -l $fileName.name` =~ /^(\d+)/;
my $NT2_DRIPPerc = $NT2_DRIPOrig / $NT2_DRIPTotal;

my @NT2_DRIP = split("\n", $NT2_DRIP);

my %NT2_DRIP;
my %temp;
my @resG;
my $totalgene = 100;
my %names;
foreach my $line (@NT2_DRIP) {
	my ($chr, $start, $end, $name, $junk1, $junk2, $junk3, $count) = split("\t", $line);
	$count ++;
	$names{$name} = 1;
	$NT2_DRIP{$count}{val} = 0 if not defined($NT2_DRIP{$count}{val});
	$NT2_DRIP{$count}{val} ++;
	$NT2_DRIP{$count}{tot} ++;
	#print "$name\t$count\n";
}
my $total = (keys %names);
my $value_ave = 0;
my ($nle, $nge, $nrun) = (0,0,0);
foreach my $count (sort {$a <=> $b} keys %NT2_DRIP) {
	my $value = $NT2_DRIP{$count}{val} / $total;
	#print "Value = $NT2_DRIP{$count}{val} / $total = $value\n";
	$nge ++ if $value >= $NT2_DRIPPerc;
	$nle ++ if $value <= $NT2_DRIPPerc;
	$nrun++;
	push(@resG, $value);
	#print "Value = $value NT2_DRIP = $NT2_DRIPOrig / $NT2_DRIPTotal = $NT2_DRIPPerc\n";
}
my $mean = mean(@resG);
my $pvalUp = ($nge + 1) / ($nrun + 1);
my $pvalDown = ($nle + 1) / ($nrun + 1);
print "$input Pval (Original > Shuffled) = $pvalUp, Shuffled_Intersect = $mean, Original_Intersect = $NT2_DRIPPerc\n" if $pvalUp < 0.05;
print "$input (Original < Shuffled) = $pvalDown, Shuffled_Intersect = $mean, Original_Intersect = $NT2_DRIPPerc\n" if $pvalDown < 0.05;
print "$input NOT SIGNIFICANT (pVal Orig > Shuff = $pvalUp, pval Orig < Shuff = $pvalDown), Shuffled_Intersect = $mean, Original_Intersect = $NT2_DRIPPerc\n" if $pvalDown > 0.05 and $pvalUp > 0.05;
__END__
my $Rscript = "

$resG
print(length(rpa))
print(length(NT2_DRIP))
df = data.frame(Original_RPA=rep($rpaPerc,length(rpa)),Shuffled_RPA=rpa,Original_NT2_DRIP=rep($NT2_DRIPPerc,length(NT2_DRIP)),Shuffled_NT2_DRIP=NT2_DRIP)
library(ggplot2)
library(reshape2)
dfm =melt(df)
print(dfm)
pdf(\"$fileName\.pdf\")
ggplot(dfm,aes(variable,value)) + geom_boxplot(aes(fill=variable))
dev.off()

";

R_toolbox::execute_Rscript($Rscript);
