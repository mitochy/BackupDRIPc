#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my @chip;
my %data;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	next if $line =~ /^Feature/;
	my ($feature, $chip1, $chip2, $value) = split("\t", $line);
	push(@chip, $chip1) if not grep(/^$chip1$/, @chip);
	push(@chip, $chip2) if not grep(/^$chip2$/, @chip);
	$data{$feature}{$chip1}{$chip2} = $value;
}
close $in1;

foreach my $feature (sort keys %data) {
	my $check = 0;
	my $Rscript = "library(GMD)\nlibrary(RColorBrewer)\n";
	foreach my $chip1 (sort @chip) {
		my @values;
		foreach my $chip2 (sort @chip) {
			my $value = $chip1 eq $chip2 ? 1 : $data{$feature}{$chip1}{$chip2};
			$value = $data{$feature}{$chip2}{$chip1} if (not defined($value));
			die "Died at $feature $chip1 $chip2\n" if not defined($value);
			push(@values, $value);
		}
		my $valueR = "c(" . join(",", @values) . ")";
		$Rscript .= "df = data.frame($chip1 = $valueR)\n" if $check == 0;
		$Rscript .= "df\$$chip1 = $valueR\n" if $check == 1;
		$check = 1;
	}
	@chip = sort(@chip);
	my $nameR = R_toolbox::newRArray(\@chip, "chip", "with_quote");
	$Rscript .= "
	$nameR
	rownames(df) = chip
	write.table(df,file=\"$feature.Rdf\")
	pdf(\"$fileName1\_$feature.pdf\")
	heatmap.3(as.matrix(df),breaks=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2),color.FUN=function(x) {rev(brewer.pal(8,\"RdBu\"))},cexRow=0.5,cexCol=0.5,main=\"$feature\")#,cellnote=as.matrix(df),notecol=\"black\")
	dev.off()
	";
	R_toolbox::execute_Rscript($Rscript);
}
