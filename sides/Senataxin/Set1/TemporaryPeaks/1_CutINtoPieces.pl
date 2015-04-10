#!/usr/bin/perl
# Previously on OldPeaks, sort and merge
use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <OldPeaks.bed>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end) = split("\t", $line);
	if ($end - $start > 1000) {
		for (my $i = $start; $i < $end - 1; $i+= 1000) {
			my $currend = $i + 999;
			print $out1 "$chr\t$i\t$currend\n";
			print "$i\t$currend\n" if ($start == 110453232);
		}
	}
	else {
		print $out1 "$chr\t$start\t$end\n";
	}
}
close $in1;
close $out1;

system("bedtools intersect -wao -a $fileName1.out -b hg19_refseq_MERGE.bed > OldPeaks_cut.bed");
system("perl -pi -e 's/^\(.+\)\t\(.+\)\t\(.+\)\t\(.+\)\t\(.+\)\t\(.+\)\t\(.+\)\t\(.+\)\t\(.+\)\t\(.+\)\$/\$1\t\$2\t\$3\t\$7\t\$8\t\$9\t\$4\t\$5\t\$6/' OldPeaks_cut.bed");

print "Output: OldPeaks_cut.bed

screen map_wig_to_bed_BIG.pl -w ../Wig/C.wig -m OldPeaks_cut.bed
screen map_wig_to_bed_BIG.pl -w ../Wig/D.wig -m OldPeaks_cut.bed
screen map_wig_to_bed_BIG.pl -w ../Wig/E.wig -m OldPeaks_cut.bed
screen map_wig_to_bed_BIG.pl -w ../Wig/F.wig -m OldPeaks_cut.bed

";
