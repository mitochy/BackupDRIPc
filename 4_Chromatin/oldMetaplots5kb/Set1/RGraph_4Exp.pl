#!/usr/bin/perl
# THIS SCRIPT IS VERY NOT USER FRIENDLY!

use strict; use warnings; use mitochy; use R_toolbox; use Getopt::Std;
use vars qw($opt_m);
getopts("m");

my ($input) = @ARGV;
die "usage: $0 [-m for mean | DEFAULT MEDIAN] <input>\n" unless @ARGV;

my ($folder, $filename) = mitochy::getFilename($input, "folder");
print "$filename\n";
my $currfeature;
my ($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+_\w+)_(drip\w*)_(\w+\_\w+)_(\w+)_(\w+)$/;
($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+)_(drip\w*)_(\w+\_\w+)_(\w+)_(\w+)$/ if not defined($driptype);
($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+_\w+)_(drip\w*)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+)_(drip\w*)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+_\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+)_(\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
$currfeature = $feature;
die "Died at pre-processing $input\n" if not defined($driptype);
print "Processing histone $histone feature $feature DRIPTYPE $driptype\n";
my @files = <./$histone*$tag\_*$feature*.temp>;
print "FILES = @files\n";

# ALL NA REMOVED
my $Rscript = "
library(ggplot2)
library(reshape2)
dat = numeric(0)
colMedian <- function(x, na.rm = TRUE) {apply(x, 2, median, na.rm =TRUE)}
";

my ($ymin, $ymax) = (0,0);
for (my $i = 0; $i < @files; $i++) {
	my $file = $files[$i];
	my ($folder2, $filename2) = mitochy::getFilename($file, "folder");
	my ($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+_\w+)_(drip\w*)_(\w+\_\w+)_(\w+)_(\w+)$/;
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+)_(drip\w*)_(\w+\_\w+)_(\w+)_(\w+)$/ if not defined($driptype);
	#die "\n\n$filename2\nHISTONE $histone TAG $tag FEAT $feature DRIP $driptype LEVEL $level\n\n" if $filename2 =~ /other/;
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+_\w+)_(drip\w*)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+)_(drip\w*)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+_\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+)_(\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
	print "$currfeature not $feature\n" and next if $currfeature ne $feature;
	next if ($level ne "zero" and $input =~ /intergenic/);
	my $id;
	$id = 0 if $level eq "super" 	and $driptype =~ /orig/i;
	$id = 1 if $level eq "high" 	and $driptype =~ /orig/i;
	$id = 2 if $level eq "med" 	and $driptype =~ /orig/i;
	$id = 3 if $level eq "low" 	and $driptype =~ /orig/i;
	$id = 4 if $level eq "super" 	and $driptype =~ /shuf/i;
	$id = 5 if $level eq "high" 	and $driptype =~ /shuf/i;
	$id = 6 if $level eq "med" 	and $driptype =~ /shuf/i;
	$id = 7 if $level eq "low" 	and $driptype =~ /shuf/i;
	print "\t$i\tProcessing file $file driptype $driptype\n";

	$Rscript .= "
	df = read.table(\"$file\"); df = df[,-1];
	length = dim(df)[2]
	";

	# USING MEAN AND SD
	#df.mean = colMeans(df)#,na.rm=TRUE)
	#e = sqrt(length(df.mean))
	#df.pos  = apply(df,2,sd)#/sqrt(length(df.mean)))
	#df.pos  = df.mean + (df.pos / e)
	#df.neg  = apply(df,2,sd)#/sqrt(length(df.mean)))
	#df.neg  = df.mean - (df.neg / e)

	# USING MEDIAN AND CONF INTERVAL
	$Rscript .= "df.mean = colMedian(df,na.rm=TRUE)\n" if not $opt_m;
	$Rscript .= "df.mean = colMeans(df,na.rm=TRUE)\n" if $opt_m;
	$Rscript .= "df.pos=numeric(0);df.neg=numeric(0);\nfor (i in 1:dim(df)[2]) {\n";
	$Rscript .= "\t\ta = apply(matrix(sample(df[,i], rep=TRUE, 10*length(df[,i])), nrow=20), 1, median,na.rm=TRUE)\n" if not $opt_m;
	$Rscript .= "\t\ta = apply(matrix(sample(df[,i], rep=TRUE, 10*length(df[,i])), nrow=20), 1, mean,na.rm=TRUE)\n" if $opt_m;
	$Rscript .= "
		df.pos[i]=as.numeric(quantile(a,0.95));
		df.neg[i]=as.numeric(quantile(a,0.05))
	}

	x = seq(-4900,4900,by=length)
	
	dm = melt(df.mean)
	colnames(dm) = c(\"Val\")
	dm\$Pos = x
	dm\$ID = \"$id\_$driptype\_$level\"
	dm\$SE.pos = df.pos
	dm\$SE.neg = df.neg
	dat = rbind(dat, dm)
	";
}

my $color = color();
$ymax = $input =~ /meth/ ? 100 : $input =~ /dens/ ? 0.8 : $input =~ /cont/ ? 0.8 : $input =~ /skew/ ?  0.15 : "max(dat\$Val,na.rm=TRUE)";
$ymin = $input =~ /meth/ ?   0 : $input =~ /dens/ ? 0.4 : $input =~ /cont/ ? 0.4 : $input =~ /skew/ ? -0.05 : "min(dat\$Val,na.rm=TRUE)";
$Rscript .= "

pdf(\"$histone\_$tag\_$feature.pdf\");
ymax = $ymax
ymin = $ymin
ymax2 = as.integer((ymax+(ymax - ymin) * 0.1) * 100)/100
ymin2 = as.integer((ymin-(ymax - ymin) * 0.1) * 100)/100
labelmax=as.integer((ymax-(ymax-ymin)*0.1)*100)/100
labelmin=as.integer((ymin+(ymax-ymin)*0.1)*100)/100
theme_blank <- theme_bw()
library(grid)
theme_blank\$axis.text=element_blank()
theme_blank\$axis.ticks=element_blank()
theme_blank\$panel.grid=element_blank()
theme_blank\$axis.title.x=element_blank()
theme_blank\$axis.title.y=element_blank()
theme_blank\$legend.position=\"none\"
theme_blank\$axis.ticks.length = unit(0,\"null\")
theme_blank\$axis.ticks.margin = unit(0,\"null\")
theme_blank\$panel.margin = unit(0,\"null\")
theme_blank\$plot.margin = rep(unit(0,\"null\"),4)
theme_blank\$panel.background = element_blank()

ggplot(dat,aes(Pos,Val)) +
geom_ribbon(aes(ymax=SE.pos,ymin=SE.neg,alpha=ID,fill=ID)) + 
geom_line(aes(fill=ID,color=ID)) +
$color
scale_alpha_manual(values=c(0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25)) +
scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
coord_cartesian(ylim=c(ymin2,ymax2)) +
theme_blank +
annotate(geom = \"segment\",x=-2000,xend=-2100,y=labelmax,yend=labelmax) +
annotate(geom = \"text\",x=-2000,y=labelmax,label=labelmax) +
annotate(geom = \"segment\",x=-2000,xend=-2100,y=labelmin,yend=labelmin)+
annotate(geom = \"text\",x=-2000,y=labelmin,label=labelmin)
dev.off()
";

open (my $out, ">", "$histone\_$tag\_$feature.R") or die "Cannot write to $histone\_$tag\_$feature.R: $!\n";
print $out $Rscript;
close $out;

system("run_Rscript.pl $histone\_$tag\_$feature.R > $histone\_$tag\_$feature.Rlog 2>&1 && tail $histone\_$tag\_$feature.Rlog");

sub color {
	my $color;
	if ($input !~ /intergenic/) {
		$color = "
		scale_color_manual(values=c(
		# Line color for 1st sample (from highest exp to lowest) - red
		rgb(140,45,4,maxColorValue=255),
		rgb(217,72,1,maxColorValue=255),
		rgb(241,105,19,maxColorValue=255),
		rgb(253,141,60,maxColorValue=255),
		#rgb(253,174,107,maxColorValue=255),
		# Line color for 2nd sample (from highest exp to lowest) - blue
		rgb(8,69,148,maxColorValue=255),
		rgb(33,113,181,maxColorValue=255),
		rgb(66,146,198,maxColorValue=255),
		#rgb(107,174,214,maxColorValue=255),
		rgb(158,202,225,maxColorValue=255))) +
		
		scale_fill_manual(values=c(
		# Fill for 1st sample (from highest exp to lowest) - red
		rgb(140,45,4,maxColorValue=255),
		rgb(217,72,1,maxColorValue=255),
		rgb(241,105,19,maxColorValue=255),
		rgb(253,141,60,maxColorValue=255),
		#rgb(253,174,107,maxColorValue=255),
		# Fill for 2nd sample (from highest exp to lowest) - blue
		rgb(8,69,148,maxColorValue=255),
		rgb(33,113,181,maxColorValue=255),
		rgb(66,146,198,maxColorValue=255),
		#rgb(107,174,214,maxColorValue=255),
		rgb(158,202,225,maxColorValue=255))) +
		";
	}
	else {
		$color = "
		scale_color_manual(values=c(
		# Line color for 1st sample (from highest exp to lowest) - red
		rgb(140,45,4,maxColorValue=255),
		# Line color for 2nd sample (from highest exp to lowest) - blue
		rgb(8,69,148,maxColorValue=255))) +
		
		scale_fill_manual(values=c(
		# Fill for 1st sample (from highest exp to lowest) - red
		rgb(140,45,4,maxColorValue=255),
		# Fill for 2nd sample (from highest exp to lowest) - blue
		rgb(8,69,148,maxColorValue=255))) +
		";
	}
	return($color);
}




__END__
#scale_fill_manual(values=c(rgb(166,54,3,maxColorValue=255),rgb(250,85,13,maxColorValue=255),rgb(8,81,156,maxColorValue=255),rgb(49,130,189,maxColorValue=255))) + 
scale_alpha_manual(values=c(0.25,0.25,0.25,0.25))

my $Rline = "";
my $distR = "";
if ($input =~ /CENTER_promoter/ or $input =~ /CENTER_terminal/ or $input =~ /CENTER_antisense/) {
	my @dist;
	my $location = $input =~ /promoter/ ? "promoter" : $input =~ /terminal/ ? "terminal" : "antisense";
	open (my $ins, "<", "dripcCENTER_$location\_orig.dist") or die "Cannot open dripcCENTER_$location\_orig.dist: $!\n";
	while (my $line = <$ins>) {
		chomp($line);
		next if $line == 0;
		push(@dist, $line);
	}
	$distR = R_toolbox::newRArray(\@dist, "dist");
	$distR .= "
	distFrame = data.frame(grp=1,dist=dist,ymin=-1,ymax=0.5)
	";
	#$Rline = "+ geom_vline(data=distFrame,aes(xintercept=dist),color=\"black\",alpha=0.01)";
	$Rline = "geom_segment(data=distFrame,aes(x=dist,xend=dist,y=ymin,yend=ymax),alpha=0.01) +";
	#$Rline = "+ geom_rect(data=distFrame,aes(xmin=distFrame\$dist,xmax=distFrame\$dist+1,ymin=distFrame\$ymin,ymax=distFrame\$ymax),size=1,alpha=0.01)";
}
