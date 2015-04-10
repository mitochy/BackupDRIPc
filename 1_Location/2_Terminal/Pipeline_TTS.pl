#!/usr/bin/perl
# This script find clean genes, find DRIPc peak at these clean genes, 
# and find the average distance of these peaks from TTS and their average length
# This is assuming that DRIPc in general will behave similarly therefore DRIPc 
# result from clean genes will be representative of DRIPc from all genes

use strict; use warnings; use mitochy; use Statistics::Basic qw(:all); use R_toolbox; use Cache::FileCache;

my $remove = 1;  # Remove this to not remove all TEMP filesuse strict; use warnings; use mitochy; use Statistics::Basic qw(:all); use R_toolbox; use Cache::FileCache;
my $bedFile   = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
my $dripcFile = "/data/mitochi/Work/Project/DRIPc/bed/dripc.bed";
my $fileType  = "terminal";
my $cache     = new Cache::FileCache;
$cache->set_cache_root("/data/mitochi/Work/Cache/");
my $x_factor  = -2000;

########################################
# 1. Find how long 5' DRIPc usually is #
# by getting all DRIPc length of genes #
# with no neighbor 5k downstream of   #
# termination site       	       #
########################################
#check_gid("TEMP_$fileType\_nogene.bed");
# a. Get 3'-end of each gene and change to 20k distance upstream and 5k downstream
runbash("bedtools_bed_change.pl -b -x -10000 -y -1 -i $bedFile -o TEMP_$fileType\_up5k.bed");
runbash("bedtools_bed_change.pl -b -x 1 -y 20000 -i $bedFile -o TEMP_$fileType\_down20k.bed");

# b. Intersect with gene to get genes that doesn't have anything within -20kb of TSS and +5kb of TSS
# Downstream is easy; just -v intersect with bedfile
runbash("bedtools intersect -v -a TEMP_$fileType\_down20k.bed -b $bedFile > TEMP_$fileType\_down20k_nogene.bed");
# Upstream cannot even have a 3'-end of gene isoform
runbash("cp TEMP_$fileType\_up5k.bed TEMP_$fileType\_up5k2.bed"); # Copy and intersect with itself
# Run bedools intersect with itself. If only intersect with itself then pass, otherwise fail
# This step (fix_bedtools) also only get genes that exist at upstream and downstream nogene.bed
runbash("bedtools intersect -wao -a TEMP_$fileType\_up5k.bed -b TEMP_$fileType\_up5k2.bed > TEMP_$fileType\_up5k_nogene.bed");
runbash("rm TEMP_$fileType\_up5k2.bed");
fix_bedtools("TEMP_$fileType\_down20k_nogene.bed", "TEMP_$fileType\_up5k_nogene.bed");
# Output is TEMP_$fileType\_nogene.bed

# c. Get genes with at least 5kb length
runbash("bedLength.pl TEMP_$fileType\_nogene.bed 10000 > TEMP_$fileType\_nogene_5kb.bed");

# Change 3'' to +/- 2kb distance (start has to be within first 3kb)
runbash("bedtools_bed_change.pl -b -x $x_factor -y 2000 -i TEMP_$fileType\_nogene_5kb.bed -o TEMP_$fileType\_nogene.bed && rm TEMP_$fileType\_nogene_5kb.bed");

# intersect with DRIPc bed file
runbash("bedtools intersect -s -wao -a TEMP_$fileType\_nogene.bed -b $dripcFile | grep -v -P \"\\t-1\\t\" > TEMP_$fileType\_nogene_dripc.bed");
# Calculate distance from TTS, size of peak, etc
my %data;
open (my $in2, "<", "TEMP_$fileType\_nogene_dripc.bed") or die "Cannot read from TEMP_$fileType\_nogene_dripc.bed: $!\n";
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

	# Say X factor is 1500
	# S---TTS---E

	# + strand
	#25000-26500--30000---35000
	#S-----TTS-----[S2----E2]
	# Start Dist: S2 - S = 30000 - 25000 = 5000. S2 to TTS is 3500 not 5000.      Therefore TTS = Dist + (-1500)
	# End Dist:   E2 - S = 35000-25000 = 10000.  E2 to TTS is 35000-26500 = 8500. Therefore TTS = Dist + (-1500)

	# - strand  
	# 20k--25000---28500--30000
	# [S2---E2]----TTS---E----------
	# Start Dist: E - E2 = 30000 - 25000 = 5000.  E2 to TTS is 3500 not 5000.    Therefore TTS = E2 + (-1500)
	# End Dist:   E - S2 = 30000 - 20000 = 10000. TTS to S2 is 28.5k-20k = 8.5k. Therefore TTS = S2 + (-1500)
	
	#if ($strand eq "+" and $start2 < $end) {
	#	$start2 = $end;
	#}
	#if ($strand eq "-" and $end2 > $start) {
	#	$end2 = $start;
	#}
	#die "start $start end $end start2 $start2 end2 $end2\n" if $line =~ /ENST00000377376.4/;
	my $startDist = $strand eq "+" ? $start2 - $start  + ($x_factor) : $end - $end2   + ($x_factor);
	my $endDist   = $strand eq "+" ? $end2   - $start  + ($x_factor) : $end - $start2 + ($x_factor);
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
summary = summary(df[which(df\$Count < 2),-1])
summary
dim(df[which(df\$Count < 2),])
dim(df)
df\$Group = df\$Count
pdf(\"TEMP_$fileType\_DRIPc.pdf\")
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

if ($remove == 1) {
        runbash("rm TEMP_*");
}

###############
# SUBROUTINES #
###############

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

sub fix_bedtools {

	# Only get genes if they also exist in downstream ($data{$name}{exist} == 1)
	# Only get genes if they don't intersect with different genes ($name1 ne $name2)

	#fix_bedtools("TEMP_$fileType\_down20k_nogene.bed", "TEMP_$fileType\_up5k_nogene.bed");
	my ($input1, $input2) = @_;
	my %data;

	my %temp;
	# Parse downstream nogene bed file
	open (my $in, "<", $input1) or die "Cannot read from $input1: $!\n";
	while (my $line = <$in>) {
		next if $line =~ /^\#/;
		my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
		$temp{$name}{exist} = 1;
	}
	close $in;

	# Parse upstream nogene bed file
	# No gene will have a "-1" coz they all will intersect with itself
	my ($folder, $filename) = mitochy::getFilename($input2, "folder");
	open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
	while (my $line = <$in2>) {
		next if $line =~ /^\#/;
		my ($chr1, $start1, $end1, $name1, $zero1, $strand1, $chr2, $start2, $end2, $name2, $zero2, $strand2, $length) = split("\t", $line);
	
		# Next if gene is not defined in downstream
		next if (not defined($temp{$name1}{exist}));
	
		# Get genes that intersect with diff gene name ($data{$name1}{num} will be 1)
		if ($name1 ne $name2) {
			$data{$name1}{num}  = 0;
			$data{$name1}{chr} = "$chr1";
			$data{$name1}{start} = "$start1";
			$data{$name1}{end} = "$end1";
			$data{$name1}{strand} = "$strand1";
		}
		else {
			if (not defined($data{$name1})) {
				$data{$name1}{chr} = "$chr1";
				$data{$name1}{start} = "$start1";
				$data{$name1}{end} = "$end1";
				$data{$name1}{strand} = "$strand1";
				$data{$name1}{num}  = 1;
			}
		}
	}
	close $in2;
	open (my $out, ">", "$filename\_TEMP") or die "Cannot write to $filename\_TEMP: $!\n";
	my %allgene = %{$cache->get("gencode_v19_annotation")};
	foreach my $name (sort {$data{$a}{chr} cmp $data{$b}{chr} || $data{$a}{start} <=> $data{$b}{start}} keys %data) {
		next if $data{$name}{num} == 0;

		# Get protein coding genes only 
		my ($gid, $type) = ("$name\_UNK", "UNKNOWN");
		if (defined($allgene{by_tid}{$name})) {
			$gid  = $allgene{by_tid}{$name};
			$type = $allgene{by_gene}{$gid}{type};
		}
		next if $type ne "protein_coding" and $type ne "lincRNA" and $type ne "antisense";
		my $chr    = $data{$name}{chr};
		my $start  = $data{$name}{start};
		my $end    = $data{$name}{end};
		my $strand = $data{$name}{strand};

		print $out "$chr\t$start\t$end\t$name\t0\t$strand\n";
	}
	close $out;

	# Convert coordinate into its real start and end
	runbash("getGene.pl $filename\_TEMP $bedFile 3");
	runbash("rm $filename\_TEMP");
	runbash("mv $filename\_TEMP.out TEMP_$fileType\_nogene.bed");
}

__END__
