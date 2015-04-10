#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox;

my ($input) = @ARGV;
die "$0 <input>\n" unless @ARGV == 1;
my ($folder, $fileName) = mitochy::getFilename($input, "folder");

my $g2hax = `bedtools intersect -wao -a $input -b g2hax.bed`;
my $rpa = `bedtools intersect -wao -a $input -b rpa.bed`;
my ($g2haxOrig)  = `bedtools intersect -u -a ../$fileName.name -b g2hax.bed | wc -l` =~ /^(\d+)/;
my ($g2haxTotal) = `wc -l ../$fileName.name` =~ /^(\d+)/;
my $g2haxPerc = $g2haxOrig / $g2haxTotal;
my ($rpaOrig)  = `bedtools intersect -u -a ../$fileName.name -b rpa.bed | wc -l` =~ /^(\d+)/;
my ($rpaTotal) = `wc -l ../$fileName.name` =~ /^(\d+)/;
my $rpaPerc = $rpaOrig / $rpaTotal;

my @rpa = split("\n", $rpa);
my @g2hax = split("\n", $g2hax);

my %g2hax;
my %temp;
my @resG;
foreach my $line (@g2hax) {
	my ($chr, $start, $end, $name, $val, $strand, $names, $dot, $num, $num2, $num3) = split("\t", $line);
	my $count = (keys %{$temp{$name}});
	$temp{$name}{$count} = $num3 == 0 ? 0 : 1;
	$g2hax{$count}{val} = 0 if not defined($g2hax{$count}{val});
	$g2hax{$count}{val} ++ if $num3 != 0;
	$g2hax{$count}{tot} ++;
}

my ($nle, $nge, $nrun) = (0,0,0);
foreach my $count (sort {$a <=> $b} keys %g2hax) {
	my $value = $g2hax{$count}{val} / $g2hax{$count}{tot};
	$nge ++ if $value >= $g2haxPerc;
	$nle ++ if $value <= $g2haxPerc;
	$nrun++;
	push(@resG, $value);
	print "Value = $value g2hax = $g2haxOrig / $g2haxOrig = $g2haxPerc\n";
}

my $pvalUp = ($nge + 1) / ($nrun + 1);
my $pvalDown = ($nle + 1) / ($nrun + 1);
print "G2HAX Pval Up = $pvalUp\nG2HAX Pval Down = $pvalDown Nrun = $nrun\n";

my %rpa;
my @resR;
%temp = ();
foreach my $line (@rpa) {
	my ($chr, $start, $end, $name, $val, $strand, $names, $dot, $num, $num2, $num3) = split("\t", $line);
	my $count = (keys %{$temp{$name}});
	$temp{$name}{$count} = $num3 == 0 ? 0 : 1;
	$rpa{$count}{val} = 0 if not defined($rpa{$count}{val});
	$rpa{$count}{val} ++ if $num3 != 0;
	$rpa{$count}{tot} ++;
}

($nle, $nge, $nrun) = (0,0,0);
foreach my $count (sort {$a <=> $b} keys %rpa) {
	my $value = $rpa{$count}{val} / $rpa{$count}{tot};
	$nge ++ if $value >= $rpaPerc;
	$nle ++ if $value <= $rpaPerc;
	$nrun++;
	print "Value = $value rpa = $rpaOrig / $rpaOrig = $rpaPerc\n";
	push(@resR, $value);
}

$pvalUp = ($nge + 1) / ($nrun + 1);
$pvalDown = ($nle + 1) / ($nrun + 1);
print "RPAPval Up = $pvalUp\nRPAPval Down = $pvalDown Nrun = $nrun\n";

my $resR = R_toolbox::newRArray(\@resR, "rpa");
my $resG = R_toolbox::newRArray(\@resG, "g2hax");

my $Rscript = "

$resR
$resG
print(length(rpa))
print(length(g2hax))
g2hax = g2hax[1:length(rpa)]
df = data.frame(Original_RPA=rep($rpaPerc,length(rpa)),Shuffled_RPA=rpa,Original_G2HAX=rep($g2haxPerc,length(g2hax)),Shuffled_G2HAX=g2hax)
library(ggplot2)
library(reshape2)
dfm =melt(df)
print(dfm)
pdf(\"$fileName\.pdf\")
ggplot(dfm,aes(variable,value)) + geom_boxplot(aes(fill=variable))
dev.off()

";

R_toolbox::execute_Rscript($Rscript);
