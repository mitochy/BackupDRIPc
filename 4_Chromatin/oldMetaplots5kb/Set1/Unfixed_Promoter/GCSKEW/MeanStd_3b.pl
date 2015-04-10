#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my (@densFile) = <./dripcNODRIP_promoter_orig_zero.out.fa.dens.tsv>;
my (@contFile) = <./dripcNODRIP_promoter_orig_zero.out.fa.cont.tsv>;
my (@skewFile) = <./dripcNODRIP_promoter_orig_zero.out.fa.skew.tsv>;
process(\@densFile, \@contFile, \@skewFile, "Promoter_Orig_zero");
(@densFile) = <./dripcNODRIP_promoter_orig_low.out.fa.dens.tsv>;
(@contFile) = <./dripcNODRIP_promoter_orig_low.out.fa.cont.tsv>;
(@skewFile) = <./dripcNODRIP_promoter_orig_low.out.fa.skew.tsv>;
process(\@densFile, \@contFile, \@skewFile, "Promoter_Orig_low");
(@densFile) = <./dripcNODRIP_promoter_orig_med.out.fa.dens.tsv>;
(@contFile) = <./dripcNODRIP_promoter_orig_med.out.fa.cont.tsv>;
(@skewFile) = <./dripcNODRIP_promoter_orig_med.out.fa.skew.tsv>;
process(\@densFile, \@contFile, \@skewFile, "Promoter_Orig_med");
(@densFile) = <./dripcNODRIP_promoter_orig_high.out.fa.dens.tsv>;
(@contFile) = <./dripcNODRIP_promoter_orig_high.out.fa.cont.tsv>;
(@skewFile) = <./dripcNODRIP_promoter_orig_high.out.fa.skew.tsv>;
process(\@densFile, \@contFile, \@skewFile, "Promoter_Orig_High");
(@densFile) = <./dripcNODRIP_promoter_orig_super.out.fa.dens.tsv>;
(@contFile) = <./dripcNODRIP_promoter_orig_super.out.fa.cont.tsv>;
(@skewFile) = <./dripcNODRIP_promoter_orig_super.out.fa.skew.tsv>;
process(\@densFile, \@contFile, \@skewFile, "Promoter_Orig_super");

sub process {
	my ($densFile, $contFile, $skewFile, $name) = @_;
	my @densFile = @{$densFile};
	my @contFile = @{$contFile};
	my @skewFile = @{$skewFile};
	
	my (@dens, @cont, @skew, @pos);

	my @densS;
	my @densM;
	my @contS;
	my @contM;
	my @skewS;
	my @skewM;
	foreach my $input (@densFile) {
		my $count = 0;
		open (my $in, "<", $input) or die;
		while (my $line = <$in>) {
			chomp($line);
			my ($pos, $val) = split("\t", $line);
			push(@{$dens[$count]}, $val);
			$pos[$count] = $pos;
			$count++;
		}
		close $in;
	}
	foreach my $input (@contFile) {
		my $count = 0;
		open (my $in, "<", $input) or die;
		while (my $line = <$in>) {
			chomp($line);
			my ($pos, $val) = split("\t", $line);
			push(@{$cont[$count]}, $val);
			$count++;
		}
		close $in;
	}
	foreach my $input (@skewFile) {
		my $count = 0;
		open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
		while (my $line = <$in>) {
			chomp($line);
			my ($pos, $val) = split("\t", $line);
			push(@{$skew[$count]}, $val);
			$count++;
		}
		close $in;
	}

	
	for (my $i = 0; $i < @dens; $i++) {
		($densM[$i], $densS[$i]) = meanStd(@{$dens[$i]});
		($contM[$i], $contS[$i]) = meanStd(@{$cont[$i]});
		($skewM[$i], $skewS[$i]) = meanStd(@{$skew[$i]});
	}
	my $densM = R_toolbox::newRArray(\@densM, "dens");	
	my $contM = R_toolbox::newRArray(\@contM, "cont");	
	my $skewM = R_toolbox::newRArray(\@skewM, "skew");	
	my $densS = R_toolbox::newRArray(\@densS, "densSD");	
	my $contS = R_toolbox::newRArray(\@contS, "contSD");	
	my $skewS = R_toolbox::newRArray(\@skewS, "skewSD");	
	my $pos   = R_toolbox::newRArray(\@pos, "x");

	my $Rscript = "
	library(ggplot2)
	library(reshape2)
	$densM
	$contM
	$skewM
	$densS
	$contS
	$skewS
	$pos

df = data.frame(dens=(dens-0)/(1.2-0)*100, 
                cont=(cont-0.3)/(0.7-0.3)*100, 
                skew=(skew+0.05)/(0.11+0.05)*100)
densMax = ((dens+densSD-0)/(1.2-0)*100)
densMin = ((dens-densSD-0)/(1.2-0)*100)
contMax = ((cont+contSD-0.3)/(0.7-0.3)*100)
contMin = ((cont-contSD-0.3)/(0.7-0.3)*100)
skewMax = ((skew+skewSD+0.05)/(0.11+0.05)*100)
skewMin = ((skew-skewSD+0.05)/(0.11+0.05)*100)
dfm=melt(df)
dfm\$x = c(x,x,x)
dfm\$ymin = c(densMin,contMin,skewMin)
dfm\$ymax = c(densMax,contMax,skewMax)
pdf(\"$name.pdf\");
#ggplot(dfm,aes(dfm\$x,value)) + stat_smooth(span=0.5,aes(color=variable),n=50,level=0.99) + scale_color_manual(values=c(\"blue4\",\"green3\",\"red3\")) + scale_fill_manual(values=c(\"blue4\",\"green3\",\"red3\")) +
ggplot(dfm,aes(dfm\$x,value)) + geom_line(aes(color=variable)) + scale_color_manual(values=c(\"blue4\",\"green3\",\"red3\")) + scale_fill_manual(values=c(\"blue4\",\"green3\",\"red3\")) +
 # theme(axis.ticks.x = element_blank(),
 #       axis.ticks.y = element_blank(),
 #       axis.text.x = element_blank(),
 #       axis.text.y = element_blank(),
 #       axis.title.x = element_blank(),
 #       axis.title.y = element_blank(),
theme(legend.position=\"none\") +
  coord_cartesian(ylim=c(0,100)) +

  geom_ribbon(aes(x=dfm\$x,ymin=ymin,ymax=ymax,fill=variable),alpha=0.25)
";

R_toolbox::execute_Rscript($Rscript);
}


sub meanStd {
	my (@data) = @_;
	my $mean = 0;
	my $sd = 0;
	for (my $i = 0; $i < @data; $i++) {
		$mean += $data[$i] / @data;
	}
	for (my $i = 0; $i < @data; $i++) {
		$sd += ($mean - $data[$i])**2/@data;
	}
	return($mean, sqrt($sd));
}
