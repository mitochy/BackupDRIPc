#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;
my ($folder, $filename) = mitochy::getFilename($input, "folder");
print "$filename\n";
my ($histone, $feature, $driptype, $level) = $filename =~ /^(\w+_\w+)_(\w+)_(\w+)_(\w+)$/;
($histone, $feature, $driptype, $level) = $filename =~ /^(\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
die "Died at pre-processing $input\n" if not defined($driptype);
print "Processing histone $histone feature $feature DRIPTYPE $driptype\n";

my @files = (<./$histone*$feature*orig*.temp>, <./$histone*$feature*shuf*.temp>);
print "FILES = @files\n";
my $Rscript = "library(ggplot2)\nlibrary(reshape2)\ndat = numeric(0)\nrowMedian <- function(x, na.rm = FALSE)\napply(x, 2, median, na.rm = na.rm)\n";
for (my $i = 0; $i < @files; $i++) {
	my $file = $files[$i];
	my ($folder2, $filename2) = mitochy::getFilename($file, "folder");
	($histone, $feature, $driptype, $level) = $filename2 =~ /^(\w+_\w+)_(\w+)_(\w+)_(\w+)$/;
	($histone, $feature, $driptype, $level) = $filename2 =~ /^(\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
	my $id = 0 if $level eq "super" and $driptype =~ /orig/i;
	$id = 1 if $level eq "high" 	and $driptype =~ /orig/i;
	$id = 2 if $level eq "med" 	and $driptype =~ /orig/i;
	$id = 3 if $level eq "low" 	and $driptype =~ /orig/i;
	$id = 4 if $level eq "zero" 	and $driptype =~ /orig/i;
	$id = 5 if $level eq "super" 	and $driptype =~ /shuf/i;
	$id = 6 if $level eq "high" 	and $driptype =~ /shuf/i;
	$id = 7 if $level eq "med" 	and $driptype =~ /shuf/i;
	$id = 8 if $level eq "low" 	and $driptype =~ /shuf/i;
	$id = 9 if $level eq "zero" 	and $driptype =~ /shuf/i;
	print "\t$i\tProcessing file $file driptype $driptype\n";
	$Rscript .= "
	df = read.table(\"$file\"); df = df[,-1]; 
	length = dim(df)[2]
	#df.mean = colMeans(df)
	df.mean.$driptype.$level = colMeans(df)
	";
	if ($file =~ /shuf/) {
	$Rscript .= "
	dat = rbind(dat,df.mean.orig.$level / df.mean.shuf.$level)
	";
	}
}


my $Rline = "";
my $distR = "";
if ($input =~ /CENTER_promoter/ or $input =~ /CENTER_terminal/) {
	my $location = $input =~ /promoter/ ? "promoter" : "terminal";
	open (my $ins, "<", "dripcCENTER_$location\_orig.dist") or die "Cannot open dripcCENTER_$location\_orig.dist: $!\n";
	my @dist;
	while (my $line = <$ins>) {
		chomp($line);
		push(@dist, $line);
	}
	$distR = R_toolbox::newRArray(\@dist, "dist");
	$distR .= "
	distFrame = data.frame(grp=1,dist=dist,ymin=-1,ymax=0)
	";
	#$Rline = "+ geom_vline(data=dist,aes(xintercept=dist),color=\"azure3\",size=0.5,alpha=0.01)";
	$Rline = "+ geom_segment(data=distFrame,aes(x=dist,xend=dist,y=ymin,yend=ymax),alpha=0.01)";
	#$Rline = "+ geom_rect(data=distFrame,aes(xmin=distFrame\$dist,xmax=distFrame\$dist+1,ymin=distFrame\$ymin,ymax=distFrame\$ymax),size=1,alpha=0.01)";
}

$Rscript .= "
dat.mean = colMeans(dat)
dat.pos  = dat.mean + (apply(dat,2,sd)/sqrt(length(dat.mean)))
dat.neg  = dat.mean - (apply(dat,2,sd)/sqrt(length(dat.mean)))
x = seq(-4900,4900,by=100)
	
dm = melt(dat.mean)
colnames(dm) = c(\"Val\")
dm\$Pos = x
dm\$ID  = 1
dm\$SE.pos = dat.pos
dm\$SE.neg = dat.neg

$distR
pdf(\"$histone\_$feature.pdf\");
ggplot(dm,aes(Pos,Val)) +
geom_ribbon(aes(ymax=SE.pos,ymin=SE.neg,fill=rgb(217,72,1,maxColorValue=255),alpha=0.25)) + 
geom_line(aes(color=rgb(217,72,1,maxColorValue=255))) + theme(legend.position=\"none\")
#scale_color_manual(values=c(rgb(140,45,4,maxColorValue=255),rgb(217,72,1,maxColorValue=255))) +
#scale_fill_manual(values=c(rgb(140,45,4,maxColorValue=255),rgb(217,72,1,maxColorValue=255))) #$Rline
dev.off()
";

open (my $out, ">", "$histone\_$feature.R") or die "Cannot write to $histone\_$feature.R: $!\n";
print $out $Rscript;
close $out;
__END__
#scale_fill_manual(values=c(rgb(166,54,3,maxColorValue=255),rgb(250,85,13,maxColorValue=255),rgb(8,81,156,maxColorValue=255),rgb(49,130,189,maxColorValue=255))) + 
scale_alpha_manual(values=c(0.25,0.25,0.25,0.25))

