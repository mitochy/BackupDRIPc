#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <FOLD.fold>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

#terminal	Tbp	6.06282400723294e-05	110	125,0,0	7291	7525	213	256	1.14096916299559

my %data;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($feature, $chip, $pval, $rna, $info, $notused, $total, $used, $notused2, $fold) = split("\t", $line);
	next if $rna == 0;
	push(@{$data{$feature}{$chip}{fold}}, $fold);
	push(@{$data{$feature}{$chip}{pval}}, $pval);
}
close $in1;

mkdir "Rdata2" if not -d "Rdata2";
foreach my $feature (sort keys %data) {
	foreach my $chip (sort keys %{$data{$feature}}) {
		my $chipmed = median(@{$data{$feature}{$chip}{fold}});
		$data{$feature}{$chip}{median} = int(100*$chipmed)/100;
	}
}

my @alphabet = qw(A1 A2 A3 A4 A5 A6 A7 A8 B C D E F G H I J K L M N O P Q R S T U V W X Y Z Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9);
foreach my $feature (sort keys %data) {
	my $Rscript .= "library(ggplot2);df = numeric(0)\n";
	my $count = 0;
	foreach my $chip (sort {$data{$feature}{$b}{median} <=> $data{$feature}{$a}{median}} keys %{$data{$feature}}) {
		my $median = $data{$feature}{$chip}{median};
		my $pval = median(@{$data{$feature}{$chip}{pval}});
		my $alp = $alphabet[$count]; die "Died at $feature $chip $count\n" if not defined($alp);
		$count ++;
		open (my $out1, ">", "Rdata2/$feature\_$chip.tsv") or die "Cannot write to Rdata2/$feature\_$chip.tsv: $!\n";
		foreach my $fold (@{$data{$feature}{$chip}{fold}}) {
			print $out1 "$fold\n";
		}
		$Rscript .= "
		temp = read.table(\"$feature\_$chip.tsv\")
		df = rbind(df,data.frame(id=\"$alp\_$chip\_$median\_$pval\", val=temp\$V1))
		";
		close $out1;
	}
	$Rscript .= "
	pdf(\"$feature.pdf\")
	ggplot(df,aes(id,val)) + geom_boxplot(aes(fill=id)) + theme(legend.text = element_text(size=3))
	dev.off()
	";
	open (my $out2, ">", "Rdata2/$feature.R") or die;
	print $out2 "$Rscript\n";
	close $out2;
}

sub median {
	my @arr = @_;
	@arr = sort {$a <=> $b} @arr;
	return($arr[int(@arr/2)]);
}
