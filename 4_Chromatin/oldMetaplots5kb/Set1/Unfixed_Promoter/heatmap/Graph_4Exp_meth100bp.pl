#!/usr/bin/perl
# THIS SCRIPT IS VERY NOT USER FRIENDLY!

use strict; use warnings; use mitochy; use R_toolbox;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;
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
my $Rscript = "library(ggplot2)\nlibrary(reshape2)\ndat = numeric(0)\ncolMedian <- function(x, na.rm = FALSE)\napply(x, 2, median, na.rm = na.rm)\n";
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
	my $id = 0 if $level eq "super" and $driptype =~ /orig/i;
	$id = 1 if $level eq "high" 	and $driptype =~ /orig/i;
	$id = 2 if $level eq "med" 	and $driptype =~ /orig/i;
	$id = 3 if $level eq "low" 	and $driptype =~ /orig/i;
	$id = 4 if $level eq "super" 	and $driptype =~ /shuf/i;
	$id = 5 if $level eq "high" 	and $driptype =~ /shuf/i;
	$id = 6 if $level eq "med" 	and $driptype =~ /shuf/i;
	$id = 7 if $level eq "low" 	and $driptype =~ /shuf/i;
	print "\t$i\tProcessing file $file driptype $driptype\n";
	my $startXseq = "-4900, 4900";
	$startXseq = "-4850, 4850" if $input =~ /dens/i or $input =~ /cont/i or $input =~ /skew/i;
	$Rscript .= "
	df = read.table(\"$file\"); df = df[,-1];
	length = dim(df)[2]
	
	# USING SD
	#df.m = colMeans(df,na.rm=TRUE)
	#e = sqrt(length(df.m))
	#df.pos  = apply(df,2,sd,na.rm=TRUE)/sqrt(length(df.m)))
	#df.pos  = df.m + (df.pos / e)
	#df.neg  = apply(df,2,sd,na.rm=TRUE)/sqrt(length(df.m)))
	#df.neg  = df.m - (df.neg / e)

	# USING CONF INTERVAL
	df.m = colMedian(df,na.rm=TRUE); # Median
	#df.m = colMeans(df,na.rm=TRUE); # Mean
	df.pos=numeric(0);
	df.neg=numeric(0);
	for (i in 1:dim(df)[2]) {
		a = apply(matrix(sample(df[,i], rep=TRUE, 10*length(df[,i])), nrow=20), 1, median,na.rm=TRUE) # Median
		#a = apply(matrix(sample(df[,i], rep=TRUE, 10*length(df[,i])), nrow=20), 1, mean,na.rm=TRUE) # Mean
		df.pos[i]=as.numeric(quantile(a,0.95))
		df.neg[i]=as.numeric(quantile(a,0.05))
	}
	x = seq($startXseq,by=100)
	
	dm = melt(df.m)
	colnames(dm) = c(\"Val\")
	dm\$Pos = x
	dm\$ID = \"$id\_$driptype\_$level\"
	dm\$SE.pos = df.pos
	dm\$SE.neg = df.neg
	#$filename2.dat = colMeans(df)
	#$filename2.dat = colMedian(df)
	dat = rbind(dat, dm)
	";
}

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
	$distR .= "distFrame = data.frame(grp=1,dist=dist,ymin=-1,ymax=0.5)\n";
	$Rline = "geom_segment(data=distFrame,aes(x=dist,xend=dist,y=ymin,yend=ymax),alpha=0.01) +";
}

# COLOR
my $color = color();
$ymax = $input =~ /dens/i ? 0.9 : $input =~ /cont/i ? 0.7   : $input =~ /skew/i ? 0.12 : $input =~ /meth/i ? 100 : "max(dat\$Val)";
$ymin = $input =~ /dens/i ? 0.1 : $input =~ /cont/i ? 0.4 : $input =~ /skew/i ? -0.02 : $input =~ /meth/i ? 0   : "min(dat\$Val)";
$Rscript .= "
$distR
pdf(\"$histone\_$tag\_$feature.pdf\");
ymax = $ymax
ymin = $ymin

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


summary(dat)
print(ymax)
ggplot(dat,aes(Pos,Val)) + #$Rline
geom_ribbon(aes(ymax=SE.pos,ymin=SE.neg,alpha=ID,fill=ID)) + 
geom_line(aes(fill=ID,color=ID)) +
$color
scale_alpha_manual(values=c(0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25)) +
xlab(\"\") + ylab(\"\") + coord_cartesian(ylim=c(ymin,ymax))# + new_theme_empty
#xlab(\"bp from $feature DRIPc Peak\") + ylab(\"$histone Signal\")
#coord_cartesian(ylim=c(0.7,1.35)) +
#theme(legend.position=\"none\") + ggtitle(\"$histone at $feature\") #+ coord_cartesian(xlim=c(-2000,2000),ylim=c(0,100)) +
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

