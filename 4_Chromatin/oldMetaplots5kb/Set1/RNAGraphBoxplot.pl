#!/usr/bin/perl
# VERSION 3 March 2015
# dripc_promoter.shuffled is in /data/mitochi/Work/Project/DRIPc/4_Chromatin/1_Shuffle/Result_125/dripc_promoter.shuffled
use strict; use warnings; use mitochy; use Getopt::Std; use R_toolbox; use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use vars qw($opt_r $opt_i $opt_a $opt_b $opt_g $opt_n);
getopts("r:i:a:b:g:n");

print RED "\n0. Sanity Check\n";
preCheck();

my ($rnaFile, $origFile, $shufFile) = ($opt_r, $opt_a, $opt_b);
my ($folder1, $fileNameOrig) = mitochy::getFilename($origFile, "folder");
my ($folder2, $fileNameShuf) = mitochy::getFilename($shufFile, "folder");

print "1. Parsing RNA file $rnaFile\n";
my %rna = %{parse_rna($rnaFile)};

my %orig;
my %shuf;
open (my $in, "<", $origFile) or die "Cannot read from $origFile: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name) = split("\t", $line);
	my $rna = defined($rna{$name}) ? $rna{$name} : 0;
	$line =~ s/\tNA/\t0/g if not $opt_n;
	push(@{$orig{super}}, $rna) if $rna >= 200;
	push(@{$orig{high}}, $rna) if $rna >= 100 and $rna < 200;
	push(@{$orig{med}}, $rna) if $rna >= 50 and $rna < 100;
	push(@{$orig{low}}, $rna) if $rna >= 10 and $rna < 50;
}
close $in;

open (my $in2, "<", $shufFile) or die "Cannot read from $shufFile: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name) = split("\t", $line);
	my $rna = defined($rna{$name}) ? $rna{$name} : 0;
	$line =~ s/\tNA/\t0/g if not $opt_n;
	push(@{$shuf{super}}, $rna) if $rna >= 200;
	push(@{$shuf{high}}, $rna) if $rna >= 100 and $rna < 200;
	push(@{$shuf{med}}, $rna) if $rna >= 50 and $rna < 100;
	push(@{$shuf{low}}, $rna) if $rna >= 10 and $rna < 50;
}
close $in2;

my $origsuper = R_toolbox::newRArray(\@{$orig{super}}, "origsuper");
my $orighigh = R_toolbox::newRArray(\@{$orig{high}}, "orighigh");
my $origmed = R_toolbox::newRArray(\@{$orig{med}}, "origmed");
my $origlow = R_toolbox::newRArray(\@{$orig{low}}, "origlow");
my $shufsuper = R_toolbox::newRArray(\@{$shuf{super}}, "shufsuper");
my $shufhigh = R_toolbox::newRArray(\@{$shuf{high}}, "shufhigh");
my $shufmed = R_toolbox::newRArray(\@{$shuf{med}}, "shufmed");
my $shuflow = R_toolbox::newRArray(\@{$shuf{low}}, "shuflow");

print "4. Running Rscript\n";
my $Rscript = "
library(ggplot2)

new_theme_empty <- theme_bw()
#new_theme_empty\$line <- element_blank()
#new_theme_empty\$rect <- element_blank()
new_theme_empty\$panel.grid= element_blank()
new_theme_empty\$axis.line = element_blank()
new_theme_empty\$axis.ticks.x = element_blank()
new_theme_empty\$axis.ticks.y = element_line(size=1)
new_theme_empty\$strip.text <- element_blank()
new_theme_empty\$axis.text.x <- element_blank()
new_theme_empty\$plot.title <- element_blank()
new_theme_empty\$axis.title <- element_blank()
new_theme_empty\$legend.position = \"none\"
#new_theme_empty\$plot.margin <- structure(c(0, 0, -1, -1), unit = \"lines\", valid.unit = 3L, class = \"unit\")

$origsuper
$orighigh
$origmed
$origlow
$shufsuper
$shufhigh
$shufmed
$shuflow
pval = wilcox.test(origsuper,shufsuper)\$p.value;print(pval)
pval = wilcox.test(origsuper,shufsuper)\$p.value;print(pval)
pval = wilcox.test(origsuper,shufsuper)\$p.value;print(pval)
pval = wilcox.test(origsuper,shufsuper)\$p.value;print(pval)

df = data.frame(id=\"0_origsuper\",val=log(origsuper))
df = rbind(df, data.frame(id=\"1_shufsuper\",val=log(shufsuper)))
df = rbind(df, data.frame(id=\"2_orighigh\",val=log(orighigh)))
df = rbind(df, data.frame(id=\"3_shufhigh\",val=log(shufhigh)))
df = rbind(df, data.frame(id=\"4_origmed\",val=log(origmed)))
df = rbind(df, data.frame(id=\"5_shufmed\",val=log(shufmed)))
df = rbind(df, data.frame(id=\"6_origlow\",val=log(origlow)))
df = rbind(df, data.frame(id=\"7_shuflow\",val=log(shuflow)))

pdf(\"RNA.pdf\")
ggplot(df,aes(id,val)) + geom_boxplot(aes(fill=id),outlier.shape=NA) + coord_cartesian(ylim=c(2,7.5)) + new_theme_empty
dev.off()
";

R_toolbox::execute_Rscript($Rscript);

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

sub parse_peak {
	my ($input) = @_;
	my %data;
	my %used;
	# Parse the bedfile
	open (my $in1, "<", $input) or die "Cannot read from $input: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
		my ($orig, $shuf) = $info =~ /ORIG=(\w+\.\d+),chr.+TWIN=(\w+\.\d+),chr/;
		die "Undefined twin gene at line: $line\n" if not defined($shuf);
		push(@{$data{orig}{$orig}{temp}}, $shuf);
		next if defined($used{$name}) and $used{$name} == 1;
		$data{orig}{$orig}{count} ++;
		$data{name}{$name} = $orig;
		$used{$name} = 1;
	}
	
	close $in1;

	# randomize array and take first 10
	foreach my $name (keys %{$data{orig}}) {
		my @value = shuffle(@{$data{orig}{$name}{temp}});
		@{$data{orig}{$name}{temp}} = ();
		for (my $i = 0; $i < @{$data{orig}{$name}{temp}}; $i++) {
			$data{shuf}{$value[$i]}{$name} ++;
		}
		@value = ();
	}
	return(\%data);
}

sub shuffle {
        my (@value) = @_;
	my $shuffleTimes = @value < 1000 ? 1000 : @value;
        for (my $i = 0; $i < $shuffleTimes; $i++) {
                my $rand1 = int(rand(@value));
                my $rand2 = int(rand(@value));
                my $val1 = $value[$rand1];
                my $val2 = $value[$rand2];
                $value[$rand1] = $val2;
                $value[$rand2] = $val1;
        }
        return(@value);
}

sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}

sub preCheck {
	my $usage = "
usage: $0 [Options] -r <RNAseq> -a <orig.tsv e.g. MCF_groseq_promoter_orig.tsv> -b <shuf.tsv>
	
E.g. $0 -r NT2.rpkm  -a dripcNODRIP_promoter_orig.input -b dripcNODRIP_promoter_shuf.input
";
	if (not defined($opt_r) or not defined($opt_a) or not defined($opt_b)) {
		print "$usage\n";
		#print GREEN "$rnaWarning";
		#print YELLOW "$fileWarning";
		die "\n";
	}
	else {
		print RED "\#################################\nSanity check success!!\n";
	}
}

__END__

		#print $outzero1  "$line\n" if $rna == 0;
		#print $outlow1   "$line\n" if $rna > 10 and $rna <= 50;
		#print $outmed1   "$line\n" if $rna > 50 and $rna <= 100;
		#print $outhigh1  "$line\n" if $rna > 100 and $rna <= 200;
		#print $outsuper1 "$line\n" if $rna > 200;


./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -i ../../1_Shuffle/Result_125/dripc_promoter.shuffled -a MCF_groseq_promoter_orig.tsv -b MCF_groseq_promoter_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -i ../../1_Shuffle/Result_125/dripc_terminal.shuffled -a MCF_groseq_terminal_orig.tsv -b MCF_groseq_terminal_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -i ../../1_Shuffle/Result_125/dripc_genebody.shuffled -a MCF_groseq_genebody_orig.tsv -b MCF_groseq_genebody_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -i ../../1_Shuffle/Result_125/dripc_antisense.shuffled -a MCF_groseq_antisense_orig.tsv -b MCF_groseq_antisense_shuf.tsv

