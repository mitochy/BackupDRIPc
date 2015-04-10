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
my $rnaFile = defined($opt_r) ? $opt_r : "../../../data/NT2.rpkm"; die RED "RNA FILE $rnaFile does not exist!\n" unless -e $rnaFile;
my ($shufFile) = $origFile; $shufFile =~ s/orig/shuf/; die "SHUFFLED FILE $shufFile does not exist!\n" if not -e $shufFile;
#my $graphScript = defined($opt_g) ? $opt_g : "RGraph_4Exp_Separate.pl";
#print RED "Graph Script File $graphScript does not exist!\n" if not -e $graphScript; die "\n" if not -e $graphScript;
my ($pdf) = $origFile =~ /^(.+)_orig/; die if not defined($pdf);

print RED "\n$pdf.pdf already exist, will not overwrite\n" and die "\n" if -e "$pdf.pdf" and not $opt_f;

my ($folder1, $fileNameOrig) = mitochy::getFilename($origFile, "folder");
my ($folder2, $fileNameShuf) = mitochy::getFilename($shufFile, "folder");

print "1. Parsing RNA file $rnaFile\n";
my %rna = %{parse_rna($rnaFile)};

print "2. Dividing Original tsv $origFile into 4 expression quartiles (e.g. $fileNameOrig\_super.temp)\n";
my %meth;
open (my $outsuper1,  ">", "$fileNameOrig\_super.temp")	or die "Cannot write to $fileNameOrig\_super.temp: $!\n";
open (my $outhigh1,   ">", "$fileNameOrig\_high.temp")  or die "Cannot write to $fileNameOrig\_high.temp: $!\n";
open (my $outmed1,    ">", "$fileNameOrig\_med.temp")	or die "Cannot write to $fileNameOrig\_med.temp: $!\n";
open (my $outlow1,    ">", "$fileNameOrig\_low.temp")	or die "Cannot write to $fileNameOrig\_low.temp: $!\n";
open (my $in, "<", $origFile) or die "Cannot read from $origFile: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $vals, $name, $value2, $strand, $info) = split("\t", $line); my $name2 = $name;
	my ($dripname) = $info =~ /DRIP=(.+);ORIG=/;
        ($name2) = $name =~ /^(ENS\w+\.\d+)_/ if $name =~ /_/;
	die "Died at $origFile undefined name2\n" if not defined($name2);
        my $rna = defined($rna{$name2}) ? $rna{$name2} : 0;
	$meth{super}{$dripname}{orig} = $vals if $rna >= 200;
	$meth{high}{$dripname}{orig}  = $vals if $rna >= 100 and $rna < 200;
	$meth{med}{$dripname}{orig}   = $vals if $rna >= 50 and $rna < 100;
	$meth{low}{$dripname}{orig}   = $vals if $rna >= 10 and $rna < 50;
}
close $in;

print "3. Dividing Shuffled tsv $shufFile into 4 expression quartiles (e.g. $fileNameShuf\_super.temp)\n";
open (my $outsuper2,  ">", "$fileNameShuf\_super.temp")  or die "Cannot write to $fileNameShuf\_super.temp: $!\n";
open (my $outhigh2,   ">", "$fileNameShuf\_high.temp")  or die "Cannot write to $fileNameShuf\_high.temp: $!\n";
open (my $outmed2,    ">", "$fileNameShuf\_med.temp")  or die "Cannot write to $fileNameShuf\_med.temp: $!\n";
open (my $outlow2,    ">", "$fileNameShuf\_low.temp")  or die "Cannot write to $fileNameShuf\_low.temp: $!\n";
open (my $in2, "<", $shufFile) or die "Cannot read from $shufFile: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $vals, $nameshuf, $value, $strand, $info) = split("\t", $line); 
	my ($dripname) = $info =~ /DRIP=(.+);ORIG=/;
	my ($name) = $info =~ /ORIG=(ENS\w+\.\d+),/;
	my $name2 = $name;
        ($name2) = $name =~ /^(ENS\w+\.\d+)_/ if $name =~ /_/;
	die "Died at $origFile undefined name2\n" if not defined($name2);
        my $rna = defined($rna{$name2}) ? $rna{$name2} : 0;
	#$line =~ s/\tNA/\t0/g if $opt_n;
	$meth{super}{$dripname}{shuf}{count} ++ if $rna >= 200 and $vals ne "NA";
	$meth{high}{$dripname}{shuf}{count}  ++ if $rna >= 100 and $rna < 200 and $vals ne "NA";
	$meth{med}{$dripname}{shuf}{count}   ++ if $rna >= 50 and $rna < 100 and $vals ne "NA";
	$meth{low}{$dripname}{shuf}{count}   ++ if $rna >= 10 and $rna < 50 and $vals ne "NA";
	$vals = 0 if $vals eq "NA";
	$meth{super}{$dripname}{shuf}{value} += $vals if $rna >= 200;
	$meth{high}{$dripname}{shuf}{value}  += $vals if $rna >= 100 and $rna < 200;
	$meth{med}{$dripname}{shuf}{value}   += $vals if $rna >= 50 and $rna < 100;
	$meth{low}{$dripname}{shuf}{value}   += $vals if $rna >= 10 and $rna < 50;
}
close $in2;

foreach my $levels (keys %meth) {
	foreach my $dripname (keys %{$meth{$levels}}) {
		my $orig = $meth{$levels}{$dripname}{orig}; next if $orig eq "NA";
		next if not defined($meth{$levels}{$dripname}{shuf}) or not defined($meth{$levels}{$dripname}{shuf}{count});
		my $shuf = $meth{$levels}{$dripname}{shuf}{value} / $meth{$levels}{$dripname}{shuf}{count};
		print $outsuper1 "$orig\n" if $levels eq "super";
		print $outsuper2 "$shuf\n" if $levels eq "super";
		print $outhigh1 "$orig\n" if $levels eq "high";
		print $outhigh2 "$shuf\n" if $levels eq "high";
		print $outmed1 "$orig\n" if $levels eq "med";
		print $outmed2 "$shuf\n" if $levels eq "med";
		print $outlow1 "$orig\n" if $levels eq "low";
		print $outlow2 "$shuf\n" if $levels eq "low";
	}
}
#print "4. Running $graphScript\n";
my $Rscript = "

library(ggplot2)
origlow=read.table(\"$fileNameOrig\_low.temp\")\$V1
shuflow=read.table(\"$fileNameShuf\_low.temp\")\$V1
origmed=read.table(\"$fileNameOrig\_med.temp\")\$V1
shufmed=read.table(\"$fileNameShuf\_med.temp\")\$V1
orighigh=read.table(\"$fileNameOrig\_high.temp\")\$V1
shufhigh=read.table(\"$fileNameShuf\_high.temp\")\$V1
origsuper=read.table(\"$fileNameOrig\_super.temp\")\$V1
shufsuper=read.table(\"$fileNameShuf\_super.temp\")\$V1

low=wilcox.test(origlow,shuflow,alternative=\"less\")\$p.value
med=wilcox.test(origmed,shufmed,alternative=\"less\")\$p.value
high=wilcox.test(orighigh,shufhigh,alternative=\"less\")\$p.value
super=wilcox.test(origsuper,shufsuper,alternative=\"less\")\$p.value

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
theme_blank\$panel.border=element_blank()

df = data.frame(id=\"orig\",rna=\"0_low\",meth=origlow)
df = rbind(df,data.frame(id=\"shuf\",rna=\"0_low\",meth=shuflow))
df = rbind(df,data.frame(id=\"orig\",rna=\"1_med\",meth=origmed))
df = rbind(df,data.frame(id=\"shuf\",rna=\"1_med\",meth=shufmed))
df = rbind(df,data.frame(id=\"orig\",rna=\"2_high\",meth=orighigh))
df = rbind(df,data.frame(id=\"shuf\",rna=\"2_high\",meth=shufhigh))
df = rbind(df,data.frame(id=\"orig\",rna=\"3_super\",meth=origsuper))
df = rbind(df,data.frame(id=\"shuf\",rna=\"3_super\",meth=shufsuper))

print(low)
print(med)
print(high)
print(super)
pdf(\"MethBoxplot.pdf\")
ggplot(df,aes(id,meth)) + geom_boxplot(aes(fill=id)) + facet_grid(.~rna) + theme_blank +
stat_summary(fun.y=mean, colour=\"darkred\", geom=\"point\", 
               shape=18, size=3,show_guide = FALSE)
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

