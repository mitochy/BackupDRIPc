#!/usr/bin/perl
# VERSION 3 March 2015
# dripc_promoter.shuffled is in /data/mitochi/Work/Project/DRIPc/4_Chromatin/1_Shuffle/Result_125/dripc_promoter.shuffled
use strict; use warnings; use mitochy; use Getopt::Std; use R_toolbox; use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use vars qw($opt_r $opt_g $opt_n $opt_m $opt_f);
getopts("r:g:nmf");

print RED "\n###################### SANITY CHECK ######################\n";
preCheck();
my ($origFile) = @ARGV; die "ORIGINAL FILE undefined or does not exist!\n" if not -e $origFile;
my $rnaFile = defined($opt_r) ? $opt_r : "../../data/NT2.rpkm"; die RED "RNA FILE $rnaFile does not exist!\n" unless -e $rnaFile;
my ($shufFile) = $origFile; $shufFile =~ s/orig/shuf/; die "SHUFFLED FILE $shufFile does not exist!\n" if not -e $shufFile;
#my $graphScript = defined($opt_g) ? $opt_g : "RGraph_4Exp_Separate.pl";
#print RED "Graph Script File $graphScript does not exist!\n" if not -e $graphScript; die "\n" if not -e $graphScript;
my ($pdf) = $origFile =~ /^(.+)_orig/; die if not defined($pdf);

#print RED "\n$pdf.pdf already exist, will not overwrite\n" and die "\n" if -e "$pdf.pdf" and not $opt_f;

my ($folder1, $fileNameOrig) = mitochy::getFilename($origFile, "folder");
my ($folder2, $fileNameShuf) = mitochy::getFilename($shufFile, "folder");
my $Rscript = "

library(ggplot2)
library(reshape2)
origlow1=read.table(\"$fileNameOrig\_low.temp\")[,seq(31,71,by=1)];colnames(origlow1) = seq(-2000,2000,by=100);origlow=melt(origlow1);colnames(origlow)=c(\"x\",\"meth\");origlow\$id=\"orig\";origlow\$rna=\"0_low\"
shuflow1=read.table(\"$fileNameShuf\_low.temp\")[,seq(31,71,by=1)];colnames(shuflow1) = seq(-2000,2000,by=100);shuflow=melt(shuflow1);colnames(shuflow)=c(\"x\",\"meth\");shuflow\$id=\"shuf\";shuflow\$rna=\"0_low\"
origmed1=read.table(\"$fileNameOrig\_med.temp\")[,seq(31,71,by=1)];colnames(origmed1) = seq(-2000,2000,by=100);origmed=melt(origmed1);colnames(origmed)=c(\"x\",\"meth\");origmed\$id=\"orig\";origmed\$rna=\"1_med\"
shufmed1=read.table(\"$fileNameShuf\_med.temp\")[,seq(31,71,by=1)];colnames(shufmed1) = seq(-2000,2000,by=100);shufmed=melt(shufmed1);colnames(shufmed)=c(\"x\",\"meth\");shufmed\$id=\"shuf\";shufmed\$rna=\"1_med\"
orighigh1=read.table(\"$fileNameOrig\_high.temp\")[,seq(31,71,by=1)];colnames(orighigh1) = seq(-2000,2000,by=100);orighigh=melt(orighigh1);colnames(orighigh)=c(\"x\",\"meth\");orighigh\$id=\"orig\";orighigh\$rna=\"2_high\"
shufhigh1=read.table(\"$fileNameShuf\_high.temp\")[,seq(31,71,by=1)];colnames(shufhigh1) = seq(-2000,2000,by=100);shufhigh=melt(shufhigh1);colnames(shufhigh)=c(\"x\",\"meth\");shufhigh\$id=\"shuf\";shufhigh\$rna=\"2_high\"
origsuper1=read.table(\"$fileNameOrig\_super.temp\")[,seq(31,71,by=1)];colnames(origsuper1) = seq(-2000,2000,by=100);origsuper=melt(origsuper1);colnames(origsuper)=c(\"x\",\"meth\");origsuper\$id=\"orig\";origsuper\$rna=\"3_super\"
shufsuper1=read.table(\"$fileNameShuf\_super.temp\")[,seq(31,71,by=1)];colnames(shufsuper1) = seq(-2000,2000,by=100);shufsuper=melt(shufsuper1);colnames(shufsuper)=c(\"x\",\"meth\");shufsuper\$id=\"shuf\";shufsuper\$rna=\"3_super\"

plow = as.numeric(0)
pmed = as.numeric(0)
phigh = as.numeric(0)
psuper = as.numeric(0)
for (i in 1:dim(origlow1)[2]) {
	pval = wilcox.test(origlow1[i][which(is.na(origlow1[1]) != TRUE),],shuflow1[i][which(is.na(shuflow1[1]) != TRUE),],alternative=\"less\",na.rm=TRUE)\$p.value
	if (pval > 0.05) {
		plow = cbind(plow,\"\")
	} else if (pval > 1E-5) {
		plow = cbind(plow,\"*\")
	} else if (pval > 1E-10) {
		plow = cbind(plow,\"**\")
	} else if (pval > 1E-20) {
		plow = cbind(plow,\"***\")
	} else {
		plow = cbind(plow,\"****\")
	}

	pval = wilcox.test(origmed1[i][which(is.na(origmed1[1]) != TRUE),],shufmed1[i][which(is.na(shufmed1[1]) != TRUE),],alternative=\"less\",na.rm=TRUE)\$p.value
	if (pval > 0.05) {
		pmed = cbind(pmed,\"\")
	} else if (pval > 1E-5) {
		pmed = cbind(pmed,\"*\")
	} else if (pval > 1E-10) {
		pmed = cbind(pmed,\"**\")
	} else if (pval > 1E-20) {
		pmed = cbind(pmed,\"***\")
	} else {
		pmed = cbind(pmed,\"****\")
	}

	pval = wilcox.test(orighigh1[i][which(is.na(orighigh1[1]) != TRUE),],shufhigh1[i][which(is.na(shufhigh1[1]) != TRUE),],alternative=\"less\",na.rm=TRUE)\$p.value
	if (pval > 0.05) {
		phigh = cbind(phigh,\"\")
	} else if (pval > 1E-5) {
		phigh = cbind(phigh,\"*\")
	} else if (pval > 1E-10) {
		phigh = cbind(phigh,\"**\")
	} else if (pval > 1E-20) {
		phigh = cbind(phigh,\"***\")
	} else {
		phigh = cbind(phigh,\"****\")
	}

	pval = wilcox.test(origsuper1[i][which(is.na(origsuper1[1]) != TRUE),],shufsuper1[i][which(is.na(shufsuper1[1]) != TRUE),],alternative=\"less\",na.rm=TRUE)\$p.value
	if (pval > 0.05) {
		psuper = cbind(psuper,\"\")
	} else if (pval > 1E-5) {
		psuper = cbind(psuper,\"*\")
	} else if (pval > 1E-10) {
		psuper = cbind(psuper,\"**\")
	} else if (pval > 1E-20) {
		psuper = cbind(psuper,\"***\")
	} else {
		psuper = cbind(psuper,\"****\")
	}
}
plow = plow[,-1];
pmed = pmed[,-1];
phigh = phigh[,-1];
psuper = psuper[,-1];

#plow = melt(plow);plow\$rna=\"0_low\"
#pmed = melt(pmed);pmed\$rna=\"1_med\"
#phigh = melt(phigh);phigh\$rna=\"2_high\"
#psuper = melt(psuper);psuper\$rna=\"3_super\"
#pvalues = rbind(plow,pmed,phigh,psuper)
#colnames(pvalues)=c(\"x\",\"pval\",\"rna\")
#origlow\$pval = plow\$value
#origmed\$pval = pmed\$value
#orighigh\$pval = phigh\$value
#origsuper\$pval = psuper\$value
#shuflow\$pval = \"\"
#shufmed\$pval = \"\"
#shufhigh\$pval = \"\"
#shufsuper\$pval = \"\"

pvalues = data.frame(x=seq(-2000,2000,by=100),pval=plow,rna=\"0_low\",id=\"orig\")
pvalues = rbind(pvalues,data.frame(x=seq(-2000,2000,by=100),pval=pmed,rna=\"1_med\",id=\"orig\"))
pvalues = rbind(pvalues,data.frame(x=seq(-2000,2000,by=100),pval=phigh,rna=\"2_high\",id=\"orig\"))
pvalues = rbind(pvalues,data.frame(x=seq(-2000,2000,by=100),pval=psuper,rna=\"3_super\",id=\"orig\"))
#pvalues\$pval = log(pvalues\$pval,base=10)

library(extrafont)
theme_blank <- theme_bw()
library(grid)
theme_blank\$axis.text=element_blank()
theme_blank\$axis.ticks=element_blank()
theme_blank\$panel.grid=element_blank()
theme_blank\$axis.title.x=element_blank()
theme_blank\$axis.title.y=element_blank()
theme_blank\$legend.position=\"none\"
theme_blank\$axis.ticks.length = unit(0,\"lines\")
theme_blank\$axis.ticks.margin = unit(0,\"lines\")
theme_blank\$panel.margin = unit(0,\"lines\")
theme_blank\$plot.margin = unit(c(0,0,0,0),\"lines\")
theme_blank\$panel.background = element_blank()
theme_blank\$strip.text = element_blank()
theme_blank\$strip.background = element_blank()
#theme_blank\$panel.border=element_blank()

df = rbind(origlow,shuflow,origmed,shufmed,orighigh,shufhigh,origsuper,shufsuper)
#df\$x = as.numeric(df\$x)
#df = df[which(df\$x > -2000 & df\$x < 2000),]
pdf(\"MethBoxplot.pdf\",width=3.75,height=0.75)
ggplot(df,aes(x,meth)) + geom_boxplot(aes(fill=id),outlier.shape=NA,lwd=0,fatten=2) + #geom_text(aes(x=x,y=100,label=pval,fill=id)) + 
scale_fill_manual(values=c(rgb(241,120,50,maxColorValue=255),rgb(141,191,220,maxColorValue=255))) +
facet_grid(rna~.,scales=\"free\") + theme_blank 
dev.off()

pdf(\"stars.pdf\",width=3.75,height=0.75)
ggplot(pvalues,aes(x,pval)) + geom_text(aes(y=100,label=pval)) + coord_cartesian(ylim=c(99,101)) + facet_grid(rna~.) #+ theme_blank
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
Usage: $0 [Options] <orig.tsv>

THIS SCRIPT REMOVE NA AND USE MEDIAN


Options:
-m: Use MEAN instead of MEDIAN
-n: Don't convert NAs into 0
-g: Use this Graph.pl file (Default: RGraph_4Exp_Separate.pl)

E.g. :
	$0 -r K562.rpkm -a H3K4me3_dripc_original.tsv -b H3K4me3_dripc_shuffled.tsv
	$0 -r NT2.rpkm  -a MCF_groseq_antisense_orig.tsv -b MCF_groseq_antisense_shuf.tsv
";

	my $rnaWarning = "
RNA seq file format:
<string name1>	<float value>
<string name2>	<float value>

E.g. 
ENST00000001	259.9
ENST00000002	500

*Your string name has to be identical as the names in -a or -b names
**Any gene name in -a or -b files not found in the RNAseq rpkm file will be assigned 0
";

	my $fileWarning = "
-a and -b file format:
<string nameA>	<float valueA1>	<float valueA2>	<float valueA3> ....
<string nameB>	<float valueB1>	<float valueB2>	<float valueB3> ....

E.g.
ENST000000001

*NAs will be treated as 0
**If number of values are not the same at each row R will be angry

";
	if (not defined($ARGV[0])) {
		print CYAN "\###################### USAGE STATEMENT ######################";
		print YELLOW "$usage\n\n";
		print CYAN "\###################### RNA FILE FORMAT ######################\n";
		print YELLOW "$rnaWarning\n";
		print CYAN "\###################### TSV FILE FORMAT ######################\n";
		print YELLOW "$fileWarning\n";
		die "\n";
	}
	else {
		print CYAN "\#################################\nSanity check success!!\n";
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

