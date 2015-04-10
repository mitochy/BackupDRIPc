#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;
my ($folder, $filename) = mitochy::getFilename($input, "folder");
print "$filename\n";
my ($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+_\w+)_(drip\w*)_(\w+)_(\w+)_(\w+)$/;
($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+)_(drip\w*)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+_\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
($histone, $tag, $feature, $driptype, $level) = $filename =~ /^(\w+)_(\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
die "Died at pre-processing $input\n" if not defined($driptype);
print "Processing histone $histone feature $feature DRIPTYPE $driptype\n";
my @files = <./$histone*$tag\_*$feature*.temp>;
print "FILES = @files\n";
my $Rscript = "library(ggplot2)\nlibrary(reshape2)\ndat = numeric(0)\nrowMedian <- function(x, na.rm = FALSE)\napply(x, 2, median, na.rm = na.rm)\n";
for (my $i = 0; $i < @files; $i++) {
	my $file = $files[$i];
	next if $file !~ /super/;
	my ($folder2, $filename2) = mitochy::getFilename($file, "folder");
	my ($histone, $feature, $driptype, $level) = $filename =~ /^(\w+_\w+)_drip\w*_(\w+)_(\w+)_(\w+)$/;
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+_\w+)_(drip\w*)_(\w+)_(\w+)_(\w+)$/;
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+)_(drip\w*)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+_\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
	($histone, $tag, $feature, $driptype, $level) = $filename2 =~ /^(\w+)_(\w+)_(\w+)_(\w+)_(\w+)_(\w+)$/ if not defined($driptype);
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
	df.mean = colMeans(df)#,na.rm=TRUE)
	#df.mean = rowMedian(df,na.rm=TRUE)
	e = sqrt(length(df.mean))
	df.pos  = apply(df,2,sd)#/sqrt(length(df.mean)))
	df.pos  = df.mean + (df.pos / e)
	df.neg  = apply(df,2,sd)#/sqrt(length(df.mean)))
	df.neg  = df.mean - (df.neg / e)
	x = seq(-4900,4900,by=100)
	
	dm = melt(df.mean)
	colnames(dm) = c(\"Val\")
	dm\$Pos = x
	dm\$ID = \"$id\_$driptype\_$level\"
	dm\$SE.pos = df.pos
	dm\$SE.neg = df.neg
	#$filename2.dat = colMeans(df)
	#$filename2.dat = rowMedian(df)
	dat = rbind(dat, dm)
	";
}

my $Rline = "";
my $distR = "";
if ($input =~ /CENTER_promoter/ or $input =~ /CENTER_terminal/) {
	my @dist;
	my $location = $input =~ /promoter/ ? "promoter" : "terminal";
	open (my $ins, "<", "dripcCENTER_$location\_orig.dist") or die "Cannot open dripcCENTER_$location\_orig.dist: $!\n";
	my @dist;
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
$Rscript .= "
$distR
pdf(\"$histone\_$tag\_$feature.pdf\");
ggplot(dat,aes(Pos,Val)) + #$Rline
geom_ribbon(aes(ymax=SE.pos,ymin=SE.neg,alpha=ID,fill=ID)) + 
geom_line(aes(fill=ID,color=ID),size=0.1) +
scale_color_manual(values=c(
rgb(140,45,4,maxColorValue=255),
rgb(8,69,148,maxColorValue=255))) +
scale_fill_manual(values=c(
rgb(140,45,4,maxColorValue=255),
rgb(8,69,148,maxColorValue=255))) +
scale_alpha_manual(values=c(0.25,0.25)) +
xlab(\"bp from $feature DRIPc Peak\") + ylab(\"$histone Signal\") + 
theme(legend.position=\"none\") + ggtitle(\"$histone at $feature\")# + coord_cartesian(xlim=c(-2000,2000),ylim=c(0,100))
dev.off()
";

open (my $out, ">", "$histone\_$tag\_$feature.R") or die "Cannot write to $histone\_$tag\_$feature.R: $!\n";
print $out $Rscript;
close $out;

system("run_Rscript.pl $histone\_$tag\_$feature.R > $histone\_$tag\_$feature.Rlog 2>&1 && tail $histone\_$tag\_$feature.Rlog");
__END__
#scale_fill_manual(values=c(rgb(166,54,3,maxColorValue=255),rgb(250,85,13,maxColorValue=255),rgb(8,81,156,maxColorValue=255),rgb(49,130,189,maxColorValue=255))) + 
scale_alpha_manual(values=c(0.25,0.25,0.25,0.25))

