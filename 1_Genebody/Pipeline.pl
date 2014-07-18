#!/usr/bin/perl
# MEDIAN DISTANCE of TTS with end of DRIPC IS 3500bp therefore this is what I'm going to use
use strict; use warnings; use mitochy; use Statistics::Basic qw(:all); use R_toolbox;

my $bedFile   = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
my $dripcFile = "/data/mitochi/Work/Project/DRIPc/bed/dripc.bed";

#5' limit: 912 to 3712 of TSS
#3' limit: 729 to 3539 of TTS
# Get 5' of all genes
runbash("bedtools_bed_change.pl -a -x -500 -y 5916 -i $bedFile -o promoter.bed");
# Get 3' of all genes
runbash("bedtools_bed_change.pl -b -x -2161 -y 6877 -i $bedFile -o terminal.bed");
# Get antisense of all genes
runbash("bedtools_bed_change.pl -a -x -7095 -y 0 -i $bedFile -o antisense.bed");
# intersect with DRIPc bed file to get promoter
runbash("echo \"track name=dripc_promoter itemRgb=On\" > dripc_promoter.bed");
runbash("echo \"track name=dripc_terminal itemRgb=On\" > dripc_terminal.bed");
runbash("echo \"track name=dripc_both itemRgb=On\" > dripc_both.bed");
runbash("echo \"track name=dripc_genebody itemRgb=On\" > dripc_genebody.bed");
runbash("echo \"track name=dripc_intergenic itemRgb=On\" > dripc_intergenic.bed");
runbash("echo \"track name=dripc_antisense itemRgb=On\" > dripc_antisense.bed");

runbash("bedtools intersect -u -s -a $dripcFile -b promoter.bed > dripc_promoter_temp.bed");
runbash("bedtools intersect -u -s -a $dripcFile -b terminal.bed >> dripc_terminal_temp.bed");
runbash("bedtools intersect -u -s -a dripc_promoter_temp.bed -b dripc_terminal_temp.bed >> dripc_both.bed");
runbash("bedtools intersect -v -s -a $dripcFile -b promoter.bed > temp");
runbash("bedtools intersect -v -s -a temp -b terminal.bed > temp2");
runbash("bedtools intersect -v -s -a temp2 -b dripc_both.bed > temp");
runbash("bedtools intersect -u -s -a temp -b $bedFile >> dripc_genebody.bed");
runbash("bedtools intersect -v -s -a dripc_promoter_temp.bed -b dripc_both.bed > temp && cat temp >> dripc_promoter.bed && rm dripc_promoter_temp.bed");
runbash("bedtools intersect -v -s -a dripc_terminal_temp.bed -b dripc_both.bed > temp && cat temp >> dripc_terminal.bed && rm dripc_terminal_temp.bed");
runbash("bedtools intersect -v -s -a $dripcFile -b dripc_both.bed > temp");
runbash("bedtools intersect -v -s -a temp -b dripc_promoter.bed > temp2");
runbash("bedtools intersect -v -s -a temp2 -b dripc_terminal.bed > temp");
runbash("bedtools intersect -v -s -a temp -b dripc_genebody.bed > temp2");
runbash("bedtools intersect -u -S -a temp2 -b antisense.bed >> dripc_antisense.bed");
runbash("bedtools intersect -v -s -a temp2 -b dripc_antisense.bed >> dripc_intergenic.bed");
runbash("rm temp");

runbash("perl -pi -e 's\/^\(.+\)\\t\(.+\)\\t\(.+\)\\t\(.+\)\\t\(.+\)\\t\(.+\)\$\/\$1\\t\$2\\t\$3\\t\$4\\t\$5\\t\$6\\t\$2\\t\$3\\tCOLOR\/' dripc*.bed");
runbash("perl -pi -e 's/(POS.+)\tCOLOR/\$1\t255,0,0/' dripc*.bed");
runbash("perl -pi -e 's/(NEG.+)\tCOLOR/\$1\t0,0,255/' dripc*.bed");
# Resolve conflicting dripc peak assignment

sub runbash {
	my ($cmd) = @_;
	system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}
