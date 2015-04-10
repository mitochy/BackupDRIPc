#!/usr/bin/perl

use strict; use warnings; use mitochy; use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

# original bed files
my (@orig) = <./*_*_orig.bed>;

#my @check;
my @check = qw(Rad21 CTCF Znf143);
#my @check = qw(Hdac2 Hdac6 H3K9ac Sin3a);
my %used;
#open (my $outFinal, ">", "Result_round2.txt") or die "Cannot open Result_round2.txt: $!\n";
#print $outFinal "Feature\tChIP1\tChIP2\tEnrichment\tOrigEnrich\tShufEnrich\tOrigObs\tShufObs\tOrigExp\tShufExp\n";
#print $outFinal "$feature1\t$chip1\t$chip2\tEnrich\t$OrigEnrichment\t$ShufEnrichment\t$OrigChIP12Total\t$ShufChIP12Total\t$OrigExp\t$ShufExp\n";#$OrigChIPTotal1\t$OrigChIPTotal2\t$FeatOrigTotal\t
foreach my $Orig1 (@orig) {
	my ($fileName1) = mitochy::getFilename($Orig1);
	my ($feature1, $chip1) = $fileName1 =~ /^(\w+)_(\w+)\_orig/; die "Died at $Orig1 undef feature1 $feature1 chip1 $chip1\n" if not defined($feature1) or not defined($chip1);
	die "Feature is not promoter/genebody/terminal ($feature1)\n" if $feature1 ne "promoter" and $feature1 ne "terminal" and $feature1 ne "genebody";

	# Next if chip1 is not in check (which we want because we don't want 150 x 150 (22k) files)
	next if not grep(/^$chip1$/i, @check) and @check != 0;

	# Next if chip1 is already used before
	if (defined($used{$feature1})) {
		next if defined($used{$feature1}{$chip1});
	}
	
	# Check if shuffle file exist
	#print GREEN "Processing $Orig1: $feature1 $chip1\n";
	my $Shuf1 = "$feature1\_$chip1\_shuf.bed"; check($Shuf1);
	my $FeatOrig = "$feature1\_orig.bed"; check($FeatOrig);
	my $FeatShuf = "$feature1\_shuf.bed"; check($FeatShuf);
	
	# Get number of line of feature with or w/o dripc total
	my ($FeatOrigTotal) = `wc -l < $FeatOrig`; chomp($FeatOrigTotal);
	my ($FeatShufTotal) = `wc -l < $FeatShuf`; chomp($FeatShufTotal);

	# Get number of line of feature with or w/o dripc that has chip1
	my ($OrigChIPTotal1) = `wc -l < $Orig1`; chomp($OrigChIPTotal1);
	my ($ShufChIPTotal1) = `wc -l < $Shuf1`; chomp($ShufChIPTotal1);

	# Mark chip1 as used
	$used{$feature1}{$chip1} = 1;

	# Now get mark #2 and intersect	(only use whatever not in @used)
	foreach my $Orig2 (@orig) {
		my ($fileName2) = mitochy::getFilename($Orig2);
		my ($feature2, $chip2) = $fileName2 =~ /^(\w+)_(\w+)\_orig/; die "Died at $Orig2 undef feature2 $feature2 chip2 $chip2\n" if not defined($feature2) or not defined($chip2);
		die "Feature is not promoter/genebody/terminal ($feature2)\n" if $feature2 ne "promoter" and $feature2 ne "terminal" and $feature2 ne "genebody";
	
		# Next if feature 2 is not feature 1 (only compare promoter with promoter)
		next if $feature1 ne $feature2;
		
		# Next if chip2 is not in check (which we want because we don't want 150 x 150 (22k) files)
		next if not grep(/^$chip2$/i, @check) and @check != 0;

		# Next if chip2 is already used before
		if (defined($used{$feature2})) {
			next if defined($used{$feature2}{"$chip2"});
			next if defined($used{$feature2}{"$chip1\_$chip2"});
			next if defined($used{$feature2}{"$chip2\_$chip1"});
		}
		
		# Check if shuffle file exist
		#print YELLOW "\tProcessing $feature1 $chip1 with $feature1 $chip2\n";
		my $Shuf2 = "$feature2\_$chip2\_shuf.bed"; check($Shuf2);
		
		# Get number of line of feature /w/o dripc that has chip2
		my ($OrigChIPTotal2) = `wc -l < $Orig2`; chomp($OrigChIPTotal2);
		my ($ShufChIPTotal2) = `wc -l < $Shuf2`; chomp($ShufChIPTotal2);

		# Mark chip2 as used
		$used{$feature2}{"$chip1\_$chip2"} = 1;

		# Now do bedtools intersect
		# Expected (at original) would be promoter_Bach1_orig.bed / promoter_orig.bed * promoter_Hdac2_orig.bed
		# Expected = OrigChipTotal1 * OrigChipTotal2 / FeatOrigTotal. Plus 1 in case denominator is 0
		my $OrigExp = int(($OrigChIPTotal1+1) * ($OrigChIPTotal2+1) / ($FeatOrigTotal+1)+1);
		my $ShufExp = int(($ShufChIPTotal1+1) * ($ShufChIPTotal2+1) / ($FeatShufTotal+1)+1);

		# Observed is bedtools intersect Orig1 and Orig2 or Shuf1 and Shuf2
		my $OrigChIP12 = "$feature1\_$chip1\_$chip2\_orig.bed2";
		my $ShufChIP12 = "$feature1\_$chip1\_$chip2\_shuf.bed2";
		`bedtools intersect -u -a $Orig1 -b $Orig2 > $OrigChIP12` if not -e $OrigChIP12 or -z $OrigChIP12 eq 1;
		`bedtools intersect -u -a $Shuf1 -b $Shuf2 > $ShufChIP12` if not -e $ShufChIP12 or -z $ShufChIP12 eq 1;
		my ($OrigChIP12Total) = `wc -l < $OrigChIP12`; chomp($OrigChIP12Total);
		my ($ShufChIP12Total) = `wc -l < $ShufChIP12`; chomp($ShufChIP12Total);

		my $OrigEnrichment = int($OrigChIP12Total / $OrigExp * 100)/100;
		my $ShufEnrichment = int($ShufChIP12Total / $ShufExp * 100)/100;
		my $Enrichment = int(($OrigEnrichment+1) / ($ShufEnrichment+1)*100)/100;
		#DEBUG#
		print GREEN "ORIG\n";
		print YELLOW "\tEXP: $feature1 $chip1 $chip2:\n\t= $Orig1 * $Orig2 / $FeatOrig\n\t= $OrigChIPTotal1 * $OrigChIPTotal2 / $FeatOrigTotal = $OrigExp\n";
		print YELLOW "\tOBS: Obs ($OrigChIP12) $OrigChIP12Total Exp = $OrigExp Enrichment = $OrigEnrichment\n";
		print GREEN "SHUF\n";
		print YELLOW "\tEXP: $feature1 $chip1 $chip2:\n\t= $Shuf1 * $Shuf2 / $FeatShuf\n\t= $ShufChIPTotal1 * $ShufChIPTotal2 / $FeatShufTotal = $ShufExp\n";
		print YELLOW "\tOBS Shuf: Obs ($ShufChIP12) $ShufChIP12Total Exp = $ShufExp Enrichment = $ShufEnrichment\n";

		if ($OrigEnrichment > 1.2 and $OrigEnrichment > $ShufEnrichment) {
			print RED "$feature1 $chip1 $chip2: ORIG $OrigEnrichment SHUF $ShufEnrichment\n";
		}
		elsif ($OrigEnrichment > $ShufEnrichment) {
			print YELLOW "$feature1 $chip1 $chip2: ORIG $OrigEnrichment SHUF $ShufEnrichment\n";
		}
		elsif ($OrigEnrichment < 0.8 and $OrigEnrichment < $ShufEnrichment) {
			print BLUE "$feature1 $chip1 $chip2: ORIG $OrigEnrichment SHUF $ShufEnrichment\n";
		}
		elsif ($OrigEnrichment < $ShufEnrichment) {
			print CYAN "$feature1 $chip1 $chip2: ORIG $OrigEnrichment SHUF $ShufEnrichment\n";
		}
		else {
			print WHITE "$feature1 $chip1 $chip2: ORIG $OrigEnrichment SHUF $ShufEnrichment\n";
		}
		print WHITE "\n\n";

		#print $outFinal "Feature=$feature1 ChIP1=$chip1 ChIP2=$chip2 FoldEnrichment=$Enrichment Orig_Enrich=$OrigEnrichment Shuf_Enrich=$ShufEnrichment Obs_Orig=$OrigChIP12Total Obs_Shuf=$ShufChIP12Total Exp_Orig=$OrigExp Exp_Shuf=$ShufExp TotalOrig_with_ChIP1=$OrigChIPTotal1 TotalOrig_with_ChIP2=$OrigChIPTotal2 TotalOrig=$FeatOrigTotal TotalShuf_with_ChIP1=$ShufChIPTotal1 TotalShuf_with_ChIP2=$ShufChIPTotal2 Total_Shuf=$FeatShufTotal\n";
	}
}

sub check {
	my ($file) = @_;
	die "$file doesn't exist\n" if not -e $file;
}

__END__
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	$line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}
close $in1;

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

