#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my @chip;
my %data;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	next if $line =~ /^Feature/;
	my @arr = split("\t", $line);
	my ($feature, $chip, $val1, $val2, $enrich) = ($arr[0],$arr[1],$arr[7], $arr[8], $arr[9]);
	my $value = 0.5 * (($val1+1) + ($val2+1)); #M
	$enrich = ($val1+1) / ($val2+1); #A
	push(@chip, $chip) if not grep(/^$chip$/, @chip);
	$data{$feature}{$chip}{M} = $value; 
	$data{$feature}{$chip}{A} = $enrich;
}
close $in1;
@chip = sort(@chip);
my $nameR = R_toolbox::newRArray(\@chip, "type", "with_quote");
my $Rscript = "library(GMD)\nlibrary(RColorBrewer)\nlibrary(ggplot2)\n$nameR\n";
my $check = 0;
foreach my $feature (sort keys %data) {
	my @values;
	my @enrich;
	foreach my $chip (sort @chip) {
		my $value = $data{$feature}{$chip}{M};
		my $enrich = $data{$feature}{$chip}{A};
		push(@values, $value);
		push(@enrich, $enrich);
	}
	my $valueR = "c(" . join(",", @values) . ")";
	my $enrichR = "c(" . join(",", @enrich) . ")";
	$Rscript .= "df = data.frame(M = $valueR, A = $enrichR, type=type, feature=\"$feature\")\n" if $check == 0;
	$Rscript .= "df = rbind(df, data.frame(M = $valueR, A = $enrichR, type=type, feature=\"$feature\"))\n" if $check == 1;
	$check = 1;
}
$Rscript .= "
#rownames(df) = c(type,type,type)
print(df)
pdf(\"$fileName1\_PeakEnrichment.pdf\")
ggplot(df, aes(M,A)) + geom_point(aes(color=type)) + geom_text(aes(y=A * 1.1, color=type,label=type),size=1) + facet_grid(feature~.) 
dev.off()
";
R_toolbox::execute_Rscript($Rscript);

__END__
#heatmap.3(as.matrix(df),breaks=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2),color.FUN=function(x) {rev(brewer.pal(8,\"RdBu\"))},cexRow=0.5,cexCol=0.5,main=\"$feature\")#,cellnote=as.matrix(df),notecol=\"black\")

