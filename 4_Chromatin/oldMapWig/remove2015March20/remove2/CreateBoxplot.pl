#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my (@fold) = <Rdata/*_fold.tsv>;
my $promR = "library(ggplot2);df = as.numeric(0)\n";
my $termR = "library(ggplot2);df = as.numeric(0)\n";
my $gbodR = "library(ggplot2);df = as.numeric(0)\n";

foreach my $fold (@fold) {

	my ($folder1, $fileName1) = mitochy::getFilename($fold,"folder");
	#my $orig = $fold; $orig =~ s/fold/orig/; die if not -e $orig;
	#my $shuf = $fold; $shuf =~ s/fold/shuf/; die if not -e $shuf;
	#print "$fold\t$orig\t$shuf\n";
	my ($chip, $feature) = $fileName1 =~ /^(\w+)_(\w+)_fold/;
	if ($feature eq "promoter") {
		$promR .= "df = rbind(df,data.frame(id=\"$chip\",value=read.table(\"$fold\")\$V1))\n";
	}
	elsif ($feature eq "terminal") {
		$termR .= "df = rbind(df,data.frame(id=\"$chip\",value=read.table(\"$fold\")\$V1))\n";
	}
	else {
		$gbodR .= "df = rbind(df,data.frame(id=\"$chip\",value=read.table(\"$fold\")\$V1))\n";
	}
	
}

$promR .= "pdf(\"Promoter.pdf\"); ggplot(df,aes(id,value)) + geom_boxplot(aes(fill=id),outlier.shape=NA); dev.off()\n";
$termR .= "pdf(\"Terminal.pdf\"); ggplot(df,aes(id,value)) + geom_boxplot(aes(fill=id),outlier.shape=NA); dev.off()\n";
$gbodR .= "pdf(\"Genebody.pdf\"); ggplot(df,aes(id,value)) + geom_boxplot(aes(fill=id),outlier.shape=NA); dev.off()\n";

R_toolbox::execute_Rscript($promR);
R_toolbox::execute_Rscript($termR);
R_toolbox::execute_Rscript($gbodR);

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
