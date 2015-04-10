#!/usr/bin/perl
# MEDIAN DISTANCE of TTS with end of drip IS 3500bp therefore this is what I'm going to use
use strict; use warnings; use mitochy; use Statistics::Basic qw(:all); use R_toolbox;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my $bedFile    = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
my $gtfFile    = "/data/mitochi/Work/Project/DRIPc/gtf/gencode_v19_annotation.gtf";
my $dripFile = "/data/mitochi/Work/Project/DRIPc/bed/drip_neat1.bed";
my ($folder2, $dripName) = mitochy::getFilename($dripFile, "folder");
my $genomeFile = "/data/mitochi/Work/Project/DRIPc/bed/hg19_genome.bed";

# promoter and terminal are +/- 2k bp of TSS/TTS (no antisense)
# promoter_ext limit :  2k to  2.5k bp of TSS
# terminal_ext limit :  2k to  2.5k bp of TTS
print BLUE "A. Getting Genomic Regions\n";
# 1. Get genomic_promoter.bed
print RED "1. Getting genomic_promoter.bed\n";
runbash("bedtools_bed_change.pl -a -x -2000 -y 2000 -i $bedFile -o genomic_promoter.bed");
runbash("bedtools_bed_change.pl -b -x 1 -y 2200 -i genomic_promoter.bed -o genomic_promoter_ext.bed");

# 2. Get genomic_terminal.bed
print RED "2. Getting genomic_terminal.bed\n";
runbash("bedtools_bed_change.pl -b -x -2000 -y 2000 -i $bedFile -o genomic_terminal.bed");
runbash("bedtools_bed_change.pl -b -x 1 -y 1700 -i genomic_terminal.bed -o genomic_terminal_ext.bed");

# 3. Getgenomic_both.bed
print RED "3. Getting genomic_both.bed\n";
runbash("bedtools intersect -s -a genomic_promoter.bed -b genomic_terminal.bed > genomic_both.bed");
runbash("bedtools subtract -s -a genomic_promoter.bed -b genomic_both.bed > temp && mv temp genomic_promoter.bed");
runbash("bedtools subtract -s -a genomic_terminal.bed -b genomic_both.bed > temp && mv temp genomic_terminal.bed");

# Promoter, Terminal, and Both are fixed
runbash("cat genomic_promoter.bed genomic_terminal.bed genomic_both.bed > Fixed.bed");

# 4. Get promoter_extended.bed
print RED "4. Getting genomic_promoter_ext.bed\n";
runbash("bedtools subtract -s -a genomic_promoter_ext.bed -b Fixed.bed > temp && mv temp genomic_promoter_ext.bed");

# 5. Get terminal_extended.bed
print RED "5. Getting genomic_terminal_ext.bed\n";
runbash("bedtools subtract -s -a genomic_terminal_ext.bed -b Fixed.bed > temp && mv temp genomic_terminal_ext.bed");

# 6. Get genomic_both_ext.bed (0 currently)
print RED "6. Getting genomic_both_ext.bed\n";
runbash("bedtools intersect -s -a genomic_promoter_ext.bed -b genomic_terminal_ext.bed > genomic_both_ext.bed");
runbash("bedtools subtract -s -a genomic_promoter_ext.bed -b genomic_both_ext.bed > temp && mv temp genomic_promoter_ext.bed");
runbash("bedtools subtract -s -a genomic_terminal_ext.bed -b genomic_both_ext.bed > temp && mv temp genomic_terminal_ext.bed");
# Promoter_ext, Terminal_ext, and Both_ext are fixed
runbash("cat genomic_promoter.bed genomic_promoter_ext.bed genomic_terminal.bed genomic_terminal_ext.bed genomic_both.bed genomic_both_ext.bed > Fixed.bed");

# 7. Get genomic_genebody.bed
print RED "7. Getting genomic_genebody.bed\n";
runbash("bedtools subtract -s -a $bedFile -b Fixed.bed > genomic_genebody.bed");

# 8. Get genomic_intergenic.bed
print RED "8. Getting genomic_intergenic.bed\n";
runbash("cat genomic_promoter.bed genomic_promoter_ext.bed genomic_terminal.bed genomic_terminal_ext.bed genomic_both.bed genomic_both_ext.bed genomic_genebody.bed > Fixed.bed");
runbash("bedtools subtract -a $genomeFile -b Fixed.bed > genomic_intergenic.bed");
runbash("rm Fixed.bed");

# 9. intersect with DRIP bed file to get promoter, terminal, and drip
print BLUE "B. Intersecting with drip Files\n";
# Intersecting Scripts
# Get Promoter
print RED "1. Getting drip_promoter.bed\n";
runbash("bedtools intersect -a $dripFile -b genomic_promoter.bed > drip_promoter.bed");
runbash("bedtools intersect -a $dripFile -b genomic_promoter_ext.bed  > drip_promoter_ext.bed");

# Get Terminal
print RED "2. Getting drip_terminal.bed\n";
runbash("bedtools intersect -a $dripFile -b genomic_terminal.bed  > drip_terminal.bed");
runbash("cat drip_promoter.bed drip_promoter_ext.bed drip_terminal.bed genomic_both.bed genomic_both_ext.bed genomic_genebody.bed > Fixeddrip.bed");
runbash("bedtools subtract -a $dripFile -b Fixeddrip.bed > DRIP_temp");
runbash("bedtools intersect -a $dripFile -b genomic_terminal_ext.bed > drip_terminal_ext_temp.bed");
runbash("bedtools intersect -u -a DRIP_temp -b genomic_terminal_ext.bed > drip_terminal_ext.bed");
# Put anything in drip_terminal_ext.bed that wasn't in genomic terminal ext.bed in there
runbash("bedtools subtract -a drip_terminal_ext.bed -b drip_terminal_ext_temp.bed >> genomic_terminal_ext.bed");
runbash("rm drip_terminal_ext_temp.bed");
# Remove any antisense* peak that overlap with the new genomic terminal
runbash("bedtools subtract -a genomic_intergenic.bed -b genomic_terminal_ext.bed > temp && mv temp genomic_intergenic.bed");
# Get Both
print RED "3. Getting drip_both.bed\n";
runbash("bedtools intersect -a $dripFile -b genomic_both.bed > drip_both.bed");
runbash("bedtools intersect -a $dripFile -b genomic_both_ext.bed > drip_both_ext.bed");
# Get Genebody
print RED "4. Getting drip_genebody.bed\n";
runbash("bedtools intersect -a $dripFile -b genomic_genebody.bed > drip_genebody.bed");
# Get Intergenic
print RED "5. Getting drip_intergenic.bed\n";
runbash("cat drip_promoter.bed drip_promoter_ext.bed drip_terminal.bed drip_terminal_ext.bed drip_both.bed drip_both_ext.bed drip_genebody.bed > Fixeddrip.bed");
runbash("bedtools subtract -a $dripFile -b Fixeddrip.bed > DRIP_temp");
runbash("bedtools intersect -a DRIP_temp -b genomic_intergenic.bed > drip_intergenic.bed");
runbash("rm Fixeddrip.bed");
# Add the region to intergenic.bed (this is treated as genomic)
runbash("rm DRIP_temp");

# Merge Genomics and sort
print BLUE "C. Merging Genomic Peaks\n";
my @genomic = qw(genomic_promoter.bed genomic_promoter_ext.bed genomic_terminal.bed genomic_terminal_ext.bed genomic_both.bed genomic_both_ext.bed genomic_genebody.bed genomic_intergenic.bed);
foreach my $genomicfile (@genomic) {
	runbash("cat $genomicfile | sort -k1,1 -k2,2n > TEMP && mv TEMP $genomicfile");
	mergeBed("$genomicfile");
}

# Merge drip, remove peaks less than 100bp, rename peaks, and add header
print BLUE "D. Merging drip Peaks\n";
my @drip = <./drip_*.bed>;
foreach my $dripFile (@drip) {
	runbash("cat $dripFile | sort -k1,1 -k2,2n > TEMP && mv TEMP $dripFile");
	mergeBed($dripFile);
 	runbash("bedLength.pl $dripFile 100 > $dripFile\_TEMP");
	renameBed("$dripFile\_TEMP");
}

# Adding track and colors
print BLUE "E. Adding tracks and colors\n";
runbash("echo track name=$dripName\_drip_promoter itemRgb=On > drip_promoter.bed");
runbash("echo track name=$dripName\_drip_promoter_ext itemRgb=On > drip_promoter_ext.bed");
runbash("echo track name=$dripName\_drip_terminal itemRgb=On > drip_terminal.bed");
runbash("echo track name=$dripName\_drip_terminal_ext itemRgb=On > drip_terminal_ext.bed");
runbash("echo track name=$dripName\_drip_both itemRgb=On > drip_both.bed");
runbash("echo track name=$dripName\_drip_both_ext itemRgb=On > drip_both_ext.bed");
runbash("echo track name=$dripName\_drip_genebody itemRgb=On > drip_genebody.bed");
runbash("echo track name=$dripName\_drip_intergenic itemRgb=On > drip_intergenic.bed");

foreach my $dripFile (@drip) {
	runbash("cat $dripFile\_TEMP >> $dripFile");
	runbash("rm $dripFile\_TEMP");
}

# Add color to peaks
#Using YoongWearn's Color Scheme (Microsoft Office's Paper)
#Promoter: 233,162,60
#Terminal: 221,182,11
#Promoter and Terminal: 161,178,147
#Gene bod: 158,133,188
#Antisense: 205, 144,163
#Intergenic: 120, 149, 179
# EXT is just add 20 to every value

runbash("perl -pi -e \'s\/^\(.+\)\\t\(.+\)\\t\(.+\)\\t\(.+\)\\t\(.+\)\\t\(.+\)\$\/\$1\\t\$2\\t\$3\\t\$4\\t\$5\\t\$6\\t\$2\\t\$3\\tCOLOR\/\' *.bed");
runbash("perl -pi -e \'s\/^\(.+\)\\t\(.+\)\\t\(.+\)\$/\$1\\t\$2\\t\$3\\tINTERGENIC\\t0\\t\+\\t\$2\\t\$3\\tCOLOR\/\' genomic_intergenic.bed");
runbash("perl -pi -e 's/\tCOLOR/\t233,162,60/' drip_promoter.bed");
runbash("perl -pi -e 's/\tCOLOR/\t255,182,80/' drip_promoter_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t221,182,11/' drip_terminal.bed");
runbash("perl -pi -e 's/\tCOLOR/\t241,202,31/' drip_terminal_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t161,178,147/' drip_both.bed");
runbash("perl -pi -e 's/\tCOLOR/\t181,198,167/' drip_both_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t158,133,188/' drip_genebody.bed");
runbash("perl -pi -e 's/\tCOLOR/\t120,149,179/' drip_intergenic.bed");

runbash("perl -pi -e 's/\tCOLOR/\t233,162,60/' genomic_promoter.bed");
runbash("perl -pi -e 's/\tCOLOR/\t255,182,80/' genomic_promoter_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t221,182,11/' genomic_terminal.bed");
runbash("perl -pi -e 's/\tCOLOR/\t241,202,31/' genomic_terminal_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t161,178,147/' genomic_both.bed");
runbash("perl -pi -e 's/\tCOLOR/\t181,198,167/' genomic_both_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t158,133,188/' genomic_genebody.bed");
runbash("perl -pi -e 's/\tCOLOR/\t120,149,179/' genomic_intergenic.bed");

# Move 
mkdir "drip"   if not -d "drip";
mkdir "Genomic" if not -d "Genomic";
runbash("mv drip_*.bed drip");
runbash("cat *.bed | grep -v track > AllGenomic.bed");
runbash("mv *.bed Genomic");
#runbash("cp Genomic/*.bed ../");
#
###############
# SUBROUTINES #
###############

sub runbash {
	my ($cmd) = @_;
	print "\t$cmd\n";
	system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}

sub mergeBed {
	my ($input) = @_;
	if ($input !~ /^genomic_intergenic.bed$/) {
		runbash("bedtools merge -nms -i $input > $input.temp");
		runbash("perl -pi -e \'s/\\t([\+\-])\$/\\t0\\t\$1/\' $input.temp");
		runbash("mv $input.temp $input");
	}
	else {
		print "INPUT $input is here\n";
		runbash("bedtools merge -i $input > $input.temp");
		runbash("mv $input.temp $input");
	}

}

sub renameBed {
	my ($input) = @_;
	runbash("cat $input | sort -k1,1 -k2,2n > TEMP && mv TEMP $input");
	my $feature = get_feature($input);
	
	open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
	open (my $out, ">", "TEMP") or die "Cannot write to TEMP: $!\n";
	my $linecount = 0;
	while (my $line = <$in>) {
		chomp($line);
		print $out "$line\n" and next if $line =~ /track/;
		$linecount++;
		my ($chr, $start, $end) = split("\t", $line);
		print $out "$chr\t$start\t$end\t$feature$linecount\t0\t+\n";
	}
	close $in;
	close $out;
	runbash("mv TEMP $input");
}

sub get_feature {
	my ($input) = @_;
	return("PROMOTER") if $input =~ /drip_promoter.bed/;
	return("PROMOTER_EXT") if $input =~ /drip_promoter_ext.bed/;
	return("TERMINAL") if $input =~ /drip_terminal.bed/;
	return("TERMINAL_EXT") if $input =~ /drip_terminal_ext.bed/;
	return("BOTH") if $input =~ /drip_both.bed/;
	return("BOTH_EXT") if $input =~ /drip_both_ext.bed/;
	return("GENEBODY") if $input =~ /drip_genebody.bed/;
	return("ANTISENSE_PROM") if $input =~ /drip_antisense.bed/;
	return("ANTISENSE_GENIC") if $input =~ /drip_antisense_other.bed/;
	return("INTERGENIC") if $input =~ /drip_intergenic.bed/;
	return("ANTISENSE_EXT") if $input =~ /drip_antisense_ext.bed/;
	die "Cannot find feature for $input\n";
}
