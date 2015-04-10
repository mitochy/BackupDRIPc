#!/usr/bin/perl
# MEDIAN DISTANCE of TTS with end of DRIPC IS 3500bp therefore this is what I'm going to use
use strict; use warnings; use mitochy; use Statistics::Basic qw(:all); use R_toolbox;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my $bedFile    = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
my $gtfFile    = "/data/mitochi/Work/Project/DRIPc/gtf/gencode_v19_annotation.gtf";
my $dripcFile  = "/data/mitochi/Work/Project/DRIPc/bed/dripc.bed";
my $genomeFile = "/data/mitochi/Work/Project/DRIPc/bed/hg19_genome.bed";

# promoter, terminal, and antisense limit are +/- 2k bp of TSS/TTS (antisense is TSS but reversed)
# promoter_ext limit :  2k to  5.2k bp of TSS
# terminal_ext limit :  2k to  3.3k bp of TTS
# antisense_ext limit: -2k to -4.5k bp of TSS
print BLUE "A. Getting Genomic Regions\n";

# 1. Get genomic_promoter.bed
print RED "1. Getting genomic_promoter.bed\n";
#runbash("bedtools_bed_change.pl -a -x 0 -y 2000 -i $bedFile -o genomic_promoter_forintersect.bed"); # For intersect with dripc purpose only
runbash("bedtools_bed_change.pl -a -x -2000 -y 2000 -i $bedFile -o genomic_promoter.bed");
#runbash("bedtools_bed_change.pl -a -x 0 -y 4000 -i $bedFile -o genomic_promoter.bed");
runbash("bedtools_bed_change.pl -b -x 1 -y 3200 -i genomic_promoter.bed -o genomic_promoter_ext.bed");

# 2. Get genomic_terminal.bed
print RED "2. Getting genomic_terminal.bed\n";
runbash("bedtools_bed_change.pl -b -x -2000 -y 2000 -i $bedFile -o genomic_terminal.bed");
runbash("bedtools_bed_change.pl -b -x 1 -y 1300 -i genomic_terminal.bed -o genomic_terminal_ext.bed");

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

# 8. Get genomic_antisense.bed
print RED "8. Getting genomic_antisense.bed\n";
runbash("cat genomic_promoter.bed genomic_promoter_ext.bed genomic_terminal.bed genomic_terminal_ext.bed genomic_both.bed genomic_both_ext.bed genomic_genebody.bed > Fixed.bed");
runbash("bedtools_bed_change.pl -a -x -2000 -y 2000 -i $bedFile -o genomic_antisense.bed");
runbash("bedtools subtract -S -a genomic_antisense.bed -b Fixed.bed > temp && mv temp genomic_antisense.bed");

# 9. Get genomic_antisense_extended.bed
print RED "9. Getting genomic_antisense_extended.bed\n";
runbash("bedtools_bed_change.pl -a -x -4500 -y -2001 -i $bedFile -o genomic_antisense_ext.bed");
runbash("bedtools subtract -S -a genomic_antisense_ext.bed -b Fixed.bed > temp");
runbash("bedtools subtract -s -a temp -b genomic_antisense.bed > genomic_antisense_ext.bed");
runbash("rm temp");
runbash("perl -pi -e \'s/\\+/NEGATIVES/\' genomic_antisense*.bed && perl -pi -e \'s/\\-/\\+/\' genomic_antisense*.bed && perl -pi -e \'s/NEGATIVES/\\-/\' genomic_antisense*.bed");

# 10. Getting genomic_antisense_other
runbash("cat genomic_terminal.bed genomic_terminal_ext.bed genomic_promoter_ext.bed genomic_genebody.bed > Fixed.bed");
runbash("bedtools merge -s -nms -i Fixed.bed | sort -k1,1 -k2,2n | perl -pi -e 's/\t([\+\-])/\t0\t$1/' > SADF");
runbash("bedtools subtract -S -a SADF -b SADF | sort -k1,1 -k2,2n |perl -pi -e 's/\+/NEG/' | perl -pi -e 's/\-/\+/' | perl -pi -e 's/\tNEG/\t\-/' > genomic_antisense_other.bed && rm SADF");
runbash("bedtools subtract -a genomic_antisense_other.bed -b genomic_antisense.bed > TEMP");
runbash("bedtools subtract -a TEMP -b genomic_antisense_ext.bed > genomic_antisense_other.bed && rm TEMP"); 

# 11. Get genomic_intergenic.bed
print RED "10. Getting genomic_intergenic.bed\n";
runbash("cat genomic_promoter.bed genomic_promoter_ext.bed genomic_terminal.bed genomic_terminal_ext.bed genomic_both.bed genomic_both_ext.bed genomic_genebody.bed genomic_antisense.bed genomic_antisense_ext.bed genomic_antisense_other.bed > Fixed.bed");
runbash("bedtools subtract -a $genomeFile -b Fixed.bed > genomic_intergenic.bed");
runbash("cat $genomeFile | perl -pi -e \'s\/\$\/\\tGENOME\\t0\\t\+\/\' > genomic_intergenic_pos.bed");
runbash("cat $genomeFile | perl -pi -e \'s\/\$\/\\tGENOME\\t0\\t\-\/\' > genomic_intergenic_neg.bed");
runbash("rm Fixed.bed");

# 11. intersect with DRIPc bed file to get promoter, terminal, and dripc
print BLUE "B. Intersecting with DRIPc Files\n";
# Intersecting Scripts
# Get Promoter
print RED "1. Getting dripc_promoter.bed\n";
runbash("bedtools intersect -s -a $dripcFile -b genomic_promoter.bed > dripc_promoter.bed");
runbash("bedtools intersect -s -a $dripcFile -b genomic_promoter_ext.bed > dripc_promoter_ext.bed");
# Get Terminal
print RED "2. Getting dripc_terminal.bed\n";
runbash("bedtools intersect -s -a $dripcFile -b genomic_terminal.bed > dripc_terminal.bed");
runbash("cat dripc_promoter.bed dripc_promoter_ext.bed dripc_terminal.bed genomic_both.bed genomic_both_ext.bed genomic_genebody.bed > FixedDripc.bed");
runbash("bedtools subtract -s -a $dripcFile -b FixedDripc.bed > dripc_temp");
runbash("bedtools intersect -s -a $dripcFile -b genomic_terminal_ext.bed > dripc_terminal_ext_temp.bed");
runbash("bedtools intersect -u -s -a dripc_temp -b genomic_terminal_ext.bed > dripc_terminal_ext.bed");
# Put anything in dripc_terminal_ext.bed that wasn't in genomic terminal ext.bed in there
runbash("bedtools subtract -s -a dripc_terminal_ext.bed -b dripc_terminal_ext_temp.bed >> genomic_terminal_ext.bed");
runbash("rm dripc_terminal_ext_temp.bed");
# Remove any antisense* peak that overlap with the new genomic terminal
runbash("bedtools subtract -s -a genomic_antisense.bed -b genomic_terminal_ext.bed > temp && mv temp genomic_antisense.bed");
runbash("bedtools subtract -s -a genomic_antisense_ext.bed -b genomic_terminal_ext.bed > temp && mv temp genomic_antisense_ext.bed");
runbash("bedtools subtract -s -a genomic_intergenic_pos.bed -b genomic_terminal_ext.bed > temp && mv temp genomic_intergenic_pos.bed");
runbash("bedtools subtract -s -a genomic_intergenic_neg.bed -b genomic_terminal_ext.bed > temp && mv temp genomic_intergenic_neg.bed");
# Get Both
print RED "3. Getting dripc_both.bed\n";
runbash("bedtools intersect -s -a $dripcFile -b genomic_both.bed > dripc_both.bed");
runbash("bedtools intersect -s -a $dripcFile -b genomic_both_ext.bed > dripc_both_ext.bed");
# Get Genebody
print RED "4. Getting dripc_genebody.bed\n";
runbash("bedtools intersect -s -a $dripcFile -b genomic_genebody.bed > dripc_genebody.bed");
# Get Antisense
print RED "5. Getting dripc_antisense.bed\n";
runbash("cat dripc_promoter.bed dripc_promoter_ext.bed dripc_terminal.bed dripc_terminal_ext.bed dripc_genebody.bed dripc_both.bed dripc_both_ext.bed> FixedDripc.bed");
runbash("bedtools subtract -s -a $dripcFile -b FixedDripc.bed > dripc_temp");
runbash("bedtools intersect -s -a dripc_temp -b genomic_antisense.bed > dripc_antisense.bed");
# Get Antisense
print RED "6. Getting dripc_antisense_ext.bed\n";
runbash("cat dripc_promoter.bed dripc_promoter_ext.bed dripc_terminal.bed dripc_terminal_ext.bed dripc_genebody.bed dripc_both.bed dripc_both_ext.bed dripc_antisense.bed > FixedDripc.bed");
runbash("bedtools subtract -s -a $dripcFile -b FixedDripc.bed > dripc_temp");
runbash("bedtools intersect -s -u -a dripc_temp -b genomic_antisense_ext.bed > dripc_antisense_ext.bed");
# Get Genic Antisense: Intersect antisense with gene region that has minus strand
print RED "7. Getting dripc_antisense_other.bed\n";
runbash("cat dripc_promoter.bed dripc_promoter_ext.bed dripc_terminal.bed dripc_terminal_ext.bed dripc_both.bed dripc_both_ext.bed dripc_genebody.bed dripc_antisense.bed dripc_antisense_ext.bed > FixedDripc.bed");
runbash("bedtools subtract -s -a $dripcFile -b FixedDripc.bed > dripc_temp");
runbash("bedtools intersect -s -a dripc_temp -b genomic_antisense_other.bed > dripc_antisense_other.bed");
# Add the region to antisense.bed (this is treated as genomic)
#runbash("bedtools intersect -S -wb -a dripc_temp -b $bedFile | cut -f1-6,10 | perl -pi -e \'s\/(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$\/\$1\\t\$2\\\t\$3\\t\$7\\t\$5\\t\$6/\' > genomic_antisense_other.bed");
# Get Intergenic
print RED "8. Getting dripc_intergenic.bed\n";
runbash("cat dripc_promoter.bed dripc_promoter_ext.bed dripc_terminal.bed dripc_terminal_ext.bed dripc_both.bed dripc_both_ext.bed dripc_genebody.bed dripc_antisense.bed dripc_antisense_ext.bed dripc_antisense_other.bed > FixedDripc.bed");
runbash("bedtools subtract -s -a $dripcFile -b FixedDripc.bed > dripc_temp");
runbash("bedtools intersect -a dripc_temp -b genomic_intergenic.bed > dripc_intergenic.bed");
runbash("bedtools intersect -s -a dripc_temp -b genomic_intergenic_pos.bed >> dripc_intergenic.bed");
runbash("bedtools intersect -s -a dripc_temp -b genomic_intergenic_neg.bed >> dripc_intergenic.bed");
runbash("rm FixedDripc.bed");
# Add the region to intergenic.bed (this is treated as genomic)
runbash("bedtools intersect -s -a dripc_temp -b genomic_intergenic_pos.bed | cut -f1-3 >> genomic_intergenic.bed");
runbash("bedtools intersect -s -a dripc_temp -b genomic_intergenic_neg.bed | cut -f1-3 >> genomic_intergenic.bed");
runbash("rm dripc_temp");
runbash("rm genomic_intergenic_*.bed");
# Merge Genomics and sort
my @genomic = qw(genomic_promoter.bed genomic_promoter_ext.bed genomic_terminal.bed genomic_terminal_ext.bed genomic_both.bed genomic_both_ext.bed genomic_genebody.bed genomic_antisense.bed genomic_antisense_ext.bed genomic_antisense_other.bed genomic_intergenic.bed);
print BLUE "C. Merging Genomic Peaks\n";
for (my $i = 0; $i < @genomic; $i++) {
	my $genomicfile = $genomic[$i];
	print RED "$i. Merging $genomicfile\n";
	runbash("cat $genomicfile | sort -k1,1 -k2,2n > TEMP && mv TEMP $genomicfile");
	mergeBed($genomicfile);
	renameBed($genomicfile) if $genomicfile =~ /intergenic/;
}
# Merge DRIPc, remove peaks less than 100bp, rename peaks, and add header
print BLUE "D. Merging DRIPc Peaks\n";
my @dripc = <./dripc_*.bed>;
for (my $i = 0; $i < @dripc; $i++) {
	my $dripcfile = $dripc[$i];
	print RED "$i. Merging $dripcfile\n";
	runbash("cat $dripcfile | sort -k1,1 -k2,2n > TEMP && mv TEMP $dripcfile");
	mergeBed($dripcfile);
 	runbash("bedLength.pl $dripcfile 100 > $dripcfile\_TEMP");
	renameBed("$dripcfile\_TEMP");
}
# Adding track and colors
print BLUE "E. Adding tracks and colors\n";
runbash("echo track name=dripc_promoter itemRgb=On > dripc_promoter.bed");
runbash("echo track name=dripc_promoter_ext itemRgb=On > dripc_promoter_ext.bed");
runbash("echo track name=dripc_terminal itemRgb=On > dripc_terminal.bed");
runbash("echo track name=dripc_terminal_ext itemRgb=On > dripc_terminal_ext.bed");
runbash("echo track name=dripc_both itemRgb=On > dripc_both.bed");
runbash("echo track name=dripc_both_ext itemRgb=On > dripc_both_ext.bed");
runbash("echo track name=dripc_genebody itemRgb=On > dripc_genebody.bed");
runbash("echo track name=dripc_antisense itemRgb=On > dripc_antisense.bed");
runbash("echo track name=dripc_antisense_other itemRgb=On > dripc_antisense_other.bed");
runbash("echo track name=dripc_antisense_ext itemRgb=On > dripc_antisense_ext.bed");
runbash("echo track name=dripc_intergenic itemRgb=On > dripc_intergenic.bed");

foreach my $dripcfile (@dripc) {
	runbash("cat $dripcfile\_TEMP >> $dripcfile");
	runbash("rm $dripcfile\_TEMP");
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
#runbash("perl -pi -e 's/^(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)$/$1\t$2\t$3\t$4\t$5\t$6\t$2\t$3\tCOLOR/' *.bed");
#runbash("perl -pi -e \'s\/^\(.+\)\\t\(.+\)\\t\(.+\)\$/\$1\\t\$2\\t\$3\\tINTERGENIC\\t0\\t\+\\t\$2\\t\$3\\tCOLOR\/\' genomic_intergenic.bed");
runbash("perl -pi -e 's/\tCOLOR/\t233,162,60/' dripc_promoter.bed");
runbash("perl -pi -e 's/\tCOLOR/\t255,182,80/' dripc_promoter_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t221,182,11/' dripc_terminal.bed");
runbash("perl -pi -e 's/\tCOLOR/\t241,202,31/' dripc_terminal_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t161,178,147/' dripc_both.bed");
runbash("perl -pi -e 's/\tCOLOR/\t181,198,167/' dripc_both_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t158,133,188/' dripc_genebody.bed");
runbash("perl -pi -e 's/\tCOLOR/\t205,144,163/' dripc_antisense.bed");
runbash("perl -pi -e 's/\tCOLOR/\t205,144,163/' dripc_antisense_other.bed");
runbash("perl -pi -e 's/\tCOLOR/\t225,146,183/' dripc_antisense_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t120,149,179/' dripc_intergenic.bed");

runbash("perl -pi -e 's/\tCOLOR/\t233,162,60/' genomic_promoter.bed");
runbash("perl -pi -e 's/\tCOLOR/\t255,182,80/' genomic_promoter_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t221,182,11/' genomic_terminal.bed");
runbash("perl -pi -e 's/\tCOLOR/\t241,202,31/' genomic_terminal_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t161,178,147/' genomic_both.bed");
runbash("perl -pi -e 's/\tCOLOR/\t181,198,167/' genomic_both_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t158,133,188/' genomic_genebody.bed");
runbash("perl -pi -e 's/\tCOLOR/\t205,144,163/' genomic_antisense.bed");
runbash("perl -pi -e 's/\tCOLOR/\t205,144,163/' genomic_antisense_other.bed");
runbash("perl -pi -e 's/\tCOLOR/\t225,146,183/' genomic_antisense_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t120,149,179/' genomic_intergenic.bed");

foreach my $genomicfile (@genomic) {
	fix_names($genomicfile);
}
foreach my $dripc (@dripc) {
	fix_names($dripc);
}

# Move 
mkdir "DRIPc"   if not -d "DRIPc";
mkdir "Genomic" if not -d "Genomic";
runbash("mv dripc_*.bed DRIPc");
runbash("cat genomic*.bed | grep -v track > AllGenomic.bed");
runbash("mv *.bed Genomic");
runbash("cp Genomic/*.bed ../");

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
		runbash("bedtools merge -s -nms -i $input > $input.temp");
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
		my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
		$strand = "+" if $input =~ /genomic_intergenic/;
		print $out "$chr\t$start\t$end\t$feature$linecount\t0\t$strand\n";
	}
	close $in;
	close $out;
	runbash("mv TEMP $input");
}

sub get_feature {
	my ($input) = @_;
	return("PROMOTER") if $input =~ /promoter.bed/;
	return("PROMOTER_EXT") if $input =~ /promoter_ext.bed/;
	return("TERMINAL") if $input =~ /terminal.bed/;
	return("TERMINAL_EXT") if $input =~ /terminal_ext.bed/;
	return("BOTH") if $input =~ /both.bed/;
	return("BOTH_EXT") if $input =~ /both_ext.bed/;
	return("GENEBODY") if $input =~ /genebody.bed/;
	return("ANTISENSE_PROM") if $input =~ /antisense.bed/;
	return("ANTISENSE_GENIC") if $input =~ /antisense_other.bed/;
	return("INTERGENIC") if $input =~ /intergenic.bed/;
	return("ANTISENSE_EXT") if $input =~ /antisense_ext.bed/;
	die "Cannot find feature for $input\n";
}

sub fix_names {
	# This fix multiple same names on the name column (e.g. ENST00000379694.4 dripc_both has 51x repeats)
	my ($input) = @_;
	open (my $in, "<", $input) or die "fix_names: Failed to open $input: $!\n";
	open (my $out, ">", "$input\_FixedNames.temp") or die "fix_names: Failed to open $input\_FixedNames.temp: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		if ($line =~ /^#/ or $line =~ /^track/) {
			print $out "$line\n";
		}
		else {
			my ($chr, $start, $end, $names, @others) = split("\t", $line);
			die "fix_names: Died due to not defined names at line $line\n" unless defined($names);
			my $others = join("\t", @others); #just in case others is not BED6
	
			# Fix multiple same names
			my @names = split(";", $names);
			my @fixed;
			for (my $i = 0; $i < @names; $i++) {
				next if $names[$i] =~ /^POS\d+$/;
				next if $names[$i] =~ /^NEG\d+$/;
				push(@fixed, $names[$i]) if not grep(/^$names[$i]$/, @fixed);
			}
			my $fixed = join(";", @fixed);
	
			# Print out
			print $out "$chr\t$start\t$end\t$fixed\t$others\n";
		}
	}
	close $in;
	close $out;
	runbash("mv $input\_FixedNames.temp $input");	
}
