#!/usr/bin/perl

use strict; use warnings; use mitochy; use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

# original bed files
my (@orig) = <./*orig.bed>;

my @check = qw(Rad21 CTCF Znf143);
my %used;
#my @check = qw(Hdac2 Hdac6 H3K9ac Sin3a);
foreach my $Orig1 (@orig) {
	my ($fileName1) = mitochy::getFilename($Orig1);
	my ($feature1, $chip1) = $fileName1 =~ /^(\w+)_(\w+)\_orig/; die "Died at $Orig1\n" if not defined($feature1) or not defined($chip1);
	die "Feature is not promoter/genebody/terminal ($feature1)\n" if $feature1 ne "promoter" and $feature1 ne "terminal" and $feature1 ne "genebody";

	# Next if chip1 is not in check (which we want because we don't want 150 x 150 (22k) files)
	next if not grep(/^$chip1$/i, @check);

	# Next if chip1 is already used before
	if (defined($used{$feature1})) {
		next if defined($used{$feature1}{$chip1});
	}
	
	# Check if shuffle file exist
	print GREEN "Processing $chip1\n";
	my $Shuf1 = "$feature1\_$chip1\_shuf.bed"; check($Shuf1);
	my $FeatOrig = "$feature1\_orig.bed"; check($FeatOrig);
	my $FeatShuf = "$feature1\_shuf.bed"; check($FeatShuf);
	
	# Get number of line of feature with or w/o dripc total
	my ($FeatOrigTotal) = `wc -l < $FeatOrig` =~ /^(\d+) /;
	# Get number of line of feature with or w/o dripc that has chip1
	my ($OrigChIPTotal1) = `wc -l $Orig1` =~ /^(\d+) /;
	my ($ShufChIPTotal1) = `wc -l $Shuf1` =~ /^(\d+) /;

	# Mark chip1 as used

	# Now get mark #2 and intersect	(only use whatever not in @used)
	foreach my $Orig2 (@orig) {
		my ($fileName2) = mitochy::getFilename($Orig2);
		my ($feature2, $chip2) = $fileName2 =~ /^(\w+)_(\w+)\_orig/; die "Died at $Orig2\n" if not defined($feature2) or not defined($chip2);
		die "Feature is not promoter/genebody/terminal ($feature2)\n" if $feature2 ne "promoter" and $feature2 ne "terminal" and $feature2 ne "genebody";
	
		# Next if feature 2 is not feature 1 (only compare promoter with promoter)
		next if $feature1 ne $feature2;
		
		# Next if chip2 is not in check (which we want because we don't want 150 x 150 (22k) files)
		next if not grep(/^$chip2$/i, @check);

		# Next if chip2 is already used before
		if (defined($used{$feature2})) {
			next if defined($used{$feature2}{$chip2});
		}
		
		# Check if shuffle file exist
		print YELLOW "\tProcessing $chip1 with $chip2\n";
		my $Shuf2 = "$feature2\_$chip2\_shuf.bed"; check($Shuf2);
		
		# Get number of line of feature /w/o dripc that has chip2
		my ($OrigChIPTotal2) = `wc -l $Orig2` =~ /^(\d+) /;
		my ($ShufChIPTotal2) = `wc -l $Shuf2` =~ /^(\d+) /;

		# Mark chip2 as used
		$used{$feature2}{$chip2} = 1;

		# Now do bedtools intersect
		# Expected (at original) would be promoter_Bach1_orig.bed / promoter_orig.bed * promoter_Hdac2_orig.bed / promoter_orig.bed
		# Expected = OrigChipTotal1 * OrigChipTotal2 / (FeatOrigTotal)**2. Plus 1 in case denominator is 0
		my $Expected = ($OrigChIPTotal1+1) * ($OrigChIPTotal2+1) / ($FeatOrigTotal+1)**2
		
		# Observed is
	}

	# Mark chip2 as used
	$used{$feature1}{$chip1} = 1;

}

sub check {
	my ($file) = @_;
	die "$file doesn't exist\n" if not -e $file;
}

__END__
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}
close $in1;

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;
