#!/usr/bin/perl
# MEDIAN DISTANCE of TTS with end of DRIPC IS 3500bp therefore this is what I'm going to use
use strict; use warnings; use mitochy; use Statistics::Basic qw(:all); use R_toolbox;

my $bedFile   = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
my $dripcFile = "/data/mitochi/Work/Project/DRIPc/bed/dripc.bed";

# First, find how long 5' DRIPc usually is by getting all DRIPc length of genes with no neighbor 10k downstream of termination site
# Merge bed file
runbash("bedtools merge -nms -s -i $bedFile \| perl -pi -e \'s\/\(.+\)\\t\(.+\)\\t\(.+\)\\t\(.+\)\\t\(.+\)\$\/\$1\\t\$2\\t\$3\\t\$4\\t0\\t\$5\/\' > merged.bed");
# Get 5' and change to 20k distance
runbash("bedtools_bed_change.pl -a -x -20000 -y -1 -i merged.bed -o promoter20k.bed");
# Intersect with gene
runbash("bedtools intersect -v -a promoter20k.bed -b $bedFile > promoter20k_nogene.bed");
# Get genes with at least 10kb length
runbash("bedLength.pl promoter20k_nogene.bed 10000 > promoter20k_nogene_10kb.bed");
# Change 5' to -1500 + 5000 distance
runbash("bedtools_bed_change.pl -x 18500 -y 5000 -i promoter20k_nogene_10kb.bed -o promoter_nogene.bed && rm promoter20k_nogene_10kb.bed && rm promoter20k.bed && rm promoter20k_nogene.bed");
# intersect with DRIPc bed file
runbash("bedtools intersect -s -wao -a promoter_nogene.bed -b $dripcFile | grep -v -P \"\\t-1\\t\" > promoter_nogene_dripc.bed");
# Calculate distance from TTS, size of peak, etc
my %data;
open (my $in2, "<", "promoter_nogene_dripc.bed") or die "Cannot read from promoter_nogene_dripc.bed: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	my ($chr, $start, $end, $name, $zero, $strand, $chr2, $start2, $end2, $name2, $zero2, $strand2, $length) = split("\t", $line);
	$name = "$chr $start $end $name";
	# Peak count
	$data{$name}{chr} = $chr;
	$data{$name}{start} = $start;
	$data{$name}{count} ++;
	
	# Peak length
	push(@{$data{$name}{length}}, abs($end2 - $start2));

	# Peak position
	my $startDist = $strand eq "+" ? $start2 - $start : $end - $end2;
	my $endDist   = $strand eq "+" ? $end2 - $start :   $end - $start2;
	push(@{$data{$name}{startdist}}, $startDist);
	push(@{$data{$name}{enddist}}, $endDist);
}
close $in2;

open (my $out, ">", "Result.tsv") or die "Cannot write to Result.tsv: $!\n";
print $out "Name\tCount\tLength\tStartDist\tEndDist\n";
foreach my $name (sort {$data{$a}{chr} cmp $data{$b}{chr} || $data{$a}{start} <=> $data{$b}{start}} keys %data) {
	my $count = $data{$name}{count};
	for (my $i = 0; $i < $count; $i++) {
		my $length = $data{$name}{length}[$i];
		my $startdist = $data{$name}{startdist}[$i];
		my $enddist = $data{$name}{enddist}[$i];
		print $out "$name\t$count\t$length\t$startdist\t$enddist\n";
	}
}
close $out;
print "Output: Result.tsv\n";

my $Rscript = "

library(ggplot2)
df = read.table(\"Result.tsv\",header=TRUE, sep=\"\t\")
df\$Group = df\$Count
pdf(\"PromoterDRIPc.pdf\")
ylim1 = boxplot.stats(df\$Length)\$stats[c(1,5)]
ggplot(df,aes(y=Length)) + geom_boxplot(aes(x=as.factor(Group), fill=Group)) + coord_cartesian(ylim = ylim1*1.05)
ylim2 = boxplot.stats(df\$StartDist)\$stats[c(1,5)]
ggplot(df,aes(y=StartDist)) + geom_boxplot(aes(x=as.factor(Group), fill=Group)) + coord_cartesian(ylim = ylim2*1.05)
ylim3 = boxplot.stats(df\$EndDist)\$stats[c(1,5)]
ggplot(df,aes(y=EndDist)) + geom_boxplot(aes(x=as.factor(Group), fill=Group)) + coord_cartesian(ylim = ylim3*1.05)
ggplot(df,aes(y=Count)) + geom_boxplot(aes(x=1))
dev.off()
";

R_toolbox::execute_Rscript($Rscript);


sub intersect {
	my ($start1, $end, $start2, $end2, $strand) = @_;
	if ($strand eq "+") {
		
	}
	else {
	}

}


sub runbash {
	my ($cmd) = @_;
	system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}















__END__





# Process Bed File
my ($input) = @_;
my %bed;
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	$bed{$chr}{$start}{$end}
}
close $in;
	
	



# Parse DRIPC data
my %dripc;
open (my $in1, "<", "terminal10k_nogene_dripc.bed") or die "Cannot read from terminal10k_nogene_dripc.bed: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	$dripc{$name}{chr} = $chr;
	$dripc{$name}{start} = $start;
	$dripc{$name}{end} = $end;
	$dripc{$name}{strand} = $strand;
}
close $in1;

print "Count\tMeanLength\tSDLength\tMeanDist\tSDDist\n";
foreach my $count (sort keys %stat) {
	for (my $i = 0; $i < $count; $i++) {
		my $meanLength = mean(@{$stat{$count}{length}});
		my $sdLength = stddev(@{$stat{$count}{length}});
		my $meanDist = mean(@{$stat{$count}{dist}});
		my $sdDist = stddev(@{$stat{$count}{dist}});
		print "$count\t$meanLength\t$sdLength\t$meanDist\t$sdDist\n";
	}
}
print "min Peak Dist = $minPeakDist\n";
print "max Peak Dist = $maxPeakDist\n";


