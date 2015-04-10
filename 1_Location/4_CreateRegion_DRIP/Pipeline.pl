#!/usr/bin/perl
# version = 1.0.0: Can be used for DRIP
# MEDIAN DISTANCE of TTS with end of DRIPC IS 3500bp therefore this is what I'm going to use
use strict; use warnings; use mitochy; use Statistics::Basic qw(:all); use R_toolbox;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my ($regionsFile)  = @ARGV;
die "Usage: $0 <Regions File>\n" unless @ARGV == 1 and defined($regionsFile);
die "\nERROR: $regionsFile does not exist!\n\n" if not -e ($regionsFile);

my ($folder, $fileName) = mitochy::getFilename($regionsFile, "folder");

# Human
my $bedFile    = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
my $genomeFile = "/data/mitochi/Work/Project/DRIPc/bed/hg19_genome.bed";

print "\n\nGENE FILE: $bedFile\nGENOME FILE: $genomeFile\n\n------------------------------\n\n";
# Mouse
#my $bedFile    = "/data/mitochi/Work/Project/DRIPc/bed/mm10_gencode.bed";
#my $genomeFile = "/data/mitochi/Work/Project/DRIPc/bed/mm10_genome.bed";
# promoter, terminal, and antisense limit are +/- 2k bp of TSS/TTS (antisense is TSS but reversed)
# promoter_ext limit :  2k to  5.2k bp of TSS
# terminal_ext limit :  2k to  3.3k bp of TTS
# antisense_ext limit: -2k to -4.5k bp of TSS
print BLUE "A. Getting Genomic $fileName\n";
# 1. Get $fileName\_genomic_promoter.bed
print RED "1. Getting $fileName\_genomic_promoter.bed\n";
runbash("bedtools_bed_change.pl -a -x -2000 -y 2000 -i $bedFile -o $fileName\_genomic_promoter.bed");
runbash("bedtools_bed_change.pl -b -x 1 -y 3200 -i $fileName\_genomic_promoter.bed -o $fileName\_genomic_promoter_ext.bed");

# 2. Get $fileName\_genomic_terminal.bed
print RED "2. Getting $fileName\_genomic_terminal.bed\n";
runbash("bedtools_bed_change.pl -b -x -2000 -y 2000 -i $bedFile -o $fileName\_genomic_terminal.bed");
runbash("bedtools_bed_change.pl -b -x 1 -y 1300 -i $fileName\_genomic_terminal.bed -o $fileName\_genomic_terminal_ext.bed");

# 3. Get$fileName\_genomic_both.bed
print RED "3. Getting $fileName\_genomic_both.bed\n";
runbash("bedtools intersect  -a $fileName\_genomic_promoter.bed -b $fileName\_genomic_terminal.bed > $fileName\_genomic_both.bed");
runbash("bedtools subtract -a $fileName\_genomic_promoter.bed -b $fileName\_genomic_both.bed > $fileName\_genomic_TEMP && mv $fileName\_genomic_TEMP $fileName\_genomic_promoter.bed");
runbash("bedtools subtract -a $fileName\_genomic_terminal.bed -b $fileName\_genomic_both.bed > $fileName\_genomic_TEMP && mv $fileName\_genomic_TEMP $fileName\_genomic_terminal.bed");

# Promoter, Terminal, and Both are fixed
runbash("cat $fileName\_genomic_promoter.bed $fileName\_genomic_terminal.bed $fileName\_genomic_both.bed > $fileName\_Fixed_genomic.bed");

# 4. Get promoter_extended.bed
print RED "4. Getting $fileName\_genomic_promoter_ext.bed\n";
runbash("bedtools subtract -a $fileName\_genomic_promoter_ext.bed -b $fileName\_Fixed_genomic.bed > $fileName\_genomic_TEMP && mv $fileName\_genomic_TEMP $fileName\_genomic_promoter_ext.bed");

# 5. Get terminal_extended.bed
print RED "5. Getting $fileName\_genomic_terminal_ext.bed\n";
runbash("bedtools subtract -a $fileName\_genomic_terminal_ext.bed -b $fileName\_Fixed_genomic.bed > $fileName\_genomic_TEMP && mv $fileName\_genomic_TEMP $fileName\_genomic_terminal_ext.bed");

# 6. Get $fileName\_genomic_both_ext.bed (0 currently)
print RED "6. Getting $fileName\_genomic_both_ext.bed\n";
runbash("bedtools intersect -a $fileName\_genomic_promoter_ext.bed -b $fileName\_genomic_terminal_ext.bed > $fileName\_genomic_both_ext.bed");
runbash("bedtools subtract -a $fileName\_genomic_promoter_ext.bed -b $fileName\_genomic_both_ext.bed > $fileName\_genomic_TEMP && mv $fileName\_genomic_TEMP $fileName\_genomic_promoter_ext.bed");
runbash("bedtools subtract -a $fileName\_genomic_terminal_ext.bed -b $fileName\_genomic_both_ext.bed > $fileName\_genomic_TEMP && mv $fileName\_genomic_TEMP $fileName\_genomic_terminal_ext.bed");
# Promoter_ext, Terminal_ext, and Both_ext are fixed
runbash("cat $fileName\_genomic_promoter.bed $fileName\_genomic_promoter_ext.bed $fileName\_genomic_terminal.bed $fileName\_genomic_terminal_ext.bed $fileName\_genomic_both.bed $fileName\_genomic_both_ext.bed > $fileName\_Fixed_genomic.bed");

# 7. Get $fileName\_genomic_genebody.bed
print RED "7. Getting $fileName\_genomic_genebody.bed\n";
runbash("bedtools subtract -a $bedFile -b $fileName\_Fixed_genomic.bed > $fileName\_genomic_genebody.bed");

# 10. Get $fileName\_genomic_intergenic.bed
print RED "10. Getting $fileName\_genomic_intergenic.bed\n";
runbash("cat $fileName\_genomic_promoter.bed $fileName\_genomic_promoter_ext.bed $fileName\_genomic_terminal.bed $fileName\_genomic_terminal_ext.bed $fileName\_genomic_both.bed $fileName\_genomic_both_ext.bed $fileName\_genomic_genebody.bed> $fileName\_Fixed_genomic.bed");
runbash("bedtools subtract -a $genomeFile -b $fileName\_Fixed_genomic.bed > $fileName\_genomic_intergenic.bed");
runbash("rm $fileName\_Fixed_genomic.bed");

# 11. intersect with Region bed file to get promoter, terminal, and dripc
print BLUE "B. Intersecting with Region Files\n";
# Intersecting Scripts

# Get Promoter
print RED "1. Getting $fileName\_region_promoter.bed\n";
runbash("bedtools intersect -a $regionsFile -b $fileName\_genomic_promoter.bed > $fileName\_region_promoter.bed");
runbash("bedtools intersect -a $regionsFile -b $fileName\_genomic_promoter_ext.bed > $fileName\_region_promoter_ext.bed");

# Get Terminal
print RED "2. Getting $fileName\_region_terminal.bed\n";
runbash("bedtools intersect -a $regionsFile -b $fileName\_genomic_terminal.bed > $fileName\_region_terminal.bed");
runbash("bedtools intersect -a $regionsFile -b $fileName\_genomic_terminal_ext.bed > $fileName\_region_terminal_ext.bed");
## NOT SURE THIS WORK IN DRIP:For terminal, sometimes ext extend way downstream of terminal, so we take the whole peak as ext (see dripc notes)
#runbash("cat $fileName\_region_promoter.bed $fileName\_region_promoter_ext.bed $fileName\_region_terminal.bed $fileName\_genomic_both.bed $fileName\_genomic_both_ext.bed $fileName\_genomic_genebody.bed > $fileName\_Fixed_genomicGenomic.bed");
#runbash("bedtools subtract -a $regionsFile -b $fileName\_Fixed_genomicGenomic.bed > $fileName\_region_TEMP");
#runbash("bedtools intersect -a $regionsFile -b $fileName\_genomic_terminal_ext.bed > $fileName\_region_terminal_ext_TEMP.bed");
#runbash("bedtools intersect -u -a $fileName\_region_temp -b $fileName\_genomic_terminal_ext.bed > $fileName\_region_terminal_ext.bed");
# Put anything in $fileName\_region_terminal_ext.bed that wasn't in genomic terminal ext.bed in there
#runbash("bedtools subtract -a $fileName\_region_terminal_ext.bed -b $fileName\_region_terminal_ext_temp.bed >> $fileName\_genomic_terminal_ext.bed");
#runbash("rm $fileName\_region_terminal_ext_temp.bed");
# Remove any intergenic genomic regions that overlap with the new genomic terminal_ext
#runbash("bedtools subtract -a $fileName\_genomic_intergenic.bed -b $fileName\_genomic_terminal_ext.bed > $fileName\_genomic_TEMP && mv $fileName\_genomic_TEMP $fileName\_genomic_intergenic.bed");

# Get Both
print RED "3. Getting $fileName\_region_both.bed\n";
runbash("bedtools intersect -a $regionsFile -b $fileName\_genomic_both.bed > $fileName\_region_both.bed");
runbash("bedtools intersect -a $regionsFile -b $fileName\_genomic_both_ext.bed > $fileName\_region_both_ext.bed");

# Get Genebody
print RED "4. Getting $fileName\_region_genebody.bed\n";
runbash("bedtools intersect -a $regionsFile -b $fileName\_genomic_genebody.bed > $fileName\_region_genebody.bed");

# Get Intergenic
print RED "8. Getting $fileName\_region_intergenic.bed\n";
runbash("cat $fileName\_region_promoter.bed $fileName\_region_promoter_ext.bed $fileName\_region_terminal.bed $fileName\_region_terminal_ext.bed $fileName\_region_both.bed $fileName\_region_both_ext.bed $fileName\_region_genebody.bed > $fileName\_Fixed_genomicGenomic.bed");
runbash("bedtools subtract -a $regionsFile -b $fileName\_Fixed_genomicGenomic.bed > $fileName\_region_temp");
runbash("bedtools intersect -a $fileName\_region_temp -b $fileName\_genomic_intergenic.bed > $fileName\_region_intergenic.bed");
runbash("rm $fileName\_Fixed_genomicGenomic.bed");
runbash("rm $fileName\_region_temp");

# Check if any files intersect with each other
check_regions_overlap();

# Merge Genomics and sort
print BLUE "C. Merging Genomic Peaks\n";
my @genomic = ("$fileName\_genomic_promoter.bed", "$fileName\_genomic_promoter_ext.bed", "$fileName\_genomic_terminal.bed", "$fileName\_genomic_terminal_ext.bed", "$fileName\_genomic_both.bed", "$fileName\_genomic_both_ext.bed", "$fileName\_genomic_genebody.bed", "$fileName\_genomic_intergenic.bed");
for (my $i = 0; $i < @genomic; $i++) {
	my $genomicfile = $genomic[$i];
	print RED "$i. Merging $genomicfile\n";
	runbash("cat $genomicfile | sort -k1,1 -k2,2n > $fileName\_genomic_TEMP && mv $fileName\_genomic_TEMP $genomicfile");
	mergeBed($genomicfile);
	renameBed($genomicfile) if $genomicfile =~ /intergenic/;
}
# Merge Region, remove peaks less than 100bp, rename peaks, and add header
print BLUE "D. Merging Region Peaks\n";
my @dripc = <./$fileName\_region_*.bed>;
for (my $i = 0; $i < @dripc; $i++) {
	my $regionsFile = $dripc[$i];
	print RED "$i. Merging $regionsFile\n";
	runbash("cat $regionsFile | sort -k1,1 -k2,2n > $fileName\_region_TEMP && mv $fileName\_region_TEMP $regionsFile");
	mergeBed($regionsFile);
 	runbash("bedLength.pl $regionsFile 100 > $regionsFile\_TEMP");
	renameBed("$regionsFile\_TEMP");
}

# Adding track and colors
print BLUE "E. Adding tracks and colors\n";
runbash("echo track name=$fileName\_region_promoter itemRgb=On > $fileName\_region_promoter.bed");
runbash("echo track name=$fileName\_region_promoter_ext itemRgb=On > $fileName\_region_promoter_ext.bed");
runbash("echo track name=$fileName\_region_terminal itemRgb=On > $fileName\_region_terminal.bed");
runbash("echo track name=$fileName\_region_terminal_ext itemRgb=On > $fileName\_region_terminal_ext.bed");
runbash("echo track name=$fileName\_region_both itemRgb=On > $fileName\_region_both.bed");
runbash("echo track name=$fileName\_region_both_ext itemRgb=On > $fileName\_region_both_ext.bed");
runbash("echo track name=$fileName\_region_genebody itemRgb=On > $fileName\_region_genebody.bed");
runbash("echo track name=$fileName\_region_intergenic itemRgb=On > $fileName\_region_intergenic.bed");

foreach my $regionsFile (@dripc) {
	runbash("cat $regionsFile\_TEMP >> $regionsFile");
	runbash("rm $regionsFile\_TEMP");
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
#runbash("perl -pi -e \'s\/^\(.+\)\\t\(.+\)\\t\(.+\)\$/\$1\\t\$2\\t\$3\\tINTERGENIC\\t0\\t\+\\t\$2\\t\$3\\tCOLOR\/\' $fileName\_genomic_intergenic.bed");
runbash("perl -pi -e 's/\tCOLOR/\t233,162,60/' $fileName\_region_promoter.bed");
runbash("perl -pi -e 's/\tCOLOR/\t255,182,80/' $fileName\_region_promoter_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t221,182,11/' $fileName\_region_terminal.bed");
runbash("perl -pi -e 's/\tCOLOR/\t241,202,31/' $fileName\_region_terminal_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t161,178,147/' $fileName\_region_both.bed");
runbash("perl -pi -e 's/\tCOLOR/\t181,198,167/' $fileName\_region_both_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t158,133,188/' $fileName\_region_genebody.bed");
runbash("perl -pi -e 's/\tCOLOR/\t120,149,179/' $fileName\_region_intergenic.bed");

runbash("perl -pi -e 's/\tCOLOR/\t233,162,60/' $fileName\_genomic_promoter.bed");
runbash("perl -pi -e 's/\tCOLOR/\t255,182,80/' $fileName\_genomic_promoter_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t221,182,11/' $fileName\_genomic_terminal.bed");
runbash("perl -pi -e 's/\tCOLOR/\t241,202,31/' $fileName\_genomic_terminal_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t161,178,147/' $fileName\_genomic_both.bed");
runbash("perl -pi -e 's/\tCOLOR/\t181,198,167/' $fileName\_genomic_both_ext.bed");
runbash("perl -pi -e 's/\tCOLOR/\t158,133,188/' $fileName\_genomic_genebody.bed");
runbash("perl -pi -e 's/\tCOLOR/\t120,149,179/' $fileName\_genomic_intergenic.bed");

foreach my $genomicfile (@genomic) {
	fix_names($genomicfile);
}
foreach my $dripc (@dripc) {
	fix_names($dripc);
}

# Move 
mkdir "Region"   if not -d "Region";
mkdir "Genomic" if not -d "Genomic";
runbash("mv $fileName\_region_*.bed Region");
runbash("cat $fileName\_genomic*.bed | grep -v track > $fileName\_AllGenomic.bed");
runbash("mv $fileName\_genomic*.bed Genomic");
#runbash("cp Genomic/$fileName\_genomic*.bed ../");

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
	if ($input !~ /^$fileName\_genomic_intergenic.bed$/) {
		runbash("bedtools merge -nms -i $input > $input.temp");
		runbash("perl -pi -e \'s/\$/\\t0\\t\+/\' $input.temp");
		runbash("mv $input.temp $input");
	}
	else {
		runbash("bedtools merge -i $input > $input.temp");
		runbash("mv $input.temp $input");
	}

}

sub renameBed {
	my ($input) = @_;
	runbash("cat $input | sort -k1,1 -k2,2n > $input.TEMP && mv $input.TEMP $input");
	my $feature = get_feature($input);
	
	open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
	open (my $out, ">", "$input.TEMP") or die "Cannot write to $input.TEMP: $!\n";
	my $linecount = 0;
	while (my $line = <$in>) {
		chomp($line);
		print $out "$line\n" and next if $line =~ /track/;
		$linecount++;
		my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
		if ($input =~ /$fileName\_genomic_intergenic/) {
			$strand = "+";
		}
		print $out "$chr\t$start\t$end\t$feature$linecount\t0\t$strand\n";
	}
	close $in;
	close $out;
	runbash("mv $input.TEMP $input");
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
	print "WARNING: Cannot find feature for $input\n";
	return("NONAME");
}

sub fix_names {
	# This fix multiple same names on the name column (e.g. ENST00000379694.4 $fileName\_region_both has 51x repeats)
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

sub check_regions_overlap {
	my @regions_check = <./$fileName\_region_*.bed>;
	my @genomic_check = <./$fileName\_genomic_*.bed>;
	foreach my $regions_check (@regions_check) {
		foreach my $regions_check2 (@regions_check) {
			next if $regions_check eq $regions_check2;
			my ($number) = `bedtools intersect -u -a $regions_check -b $regions_check2 | wc -l` =~ /^(\d+)/;
			die "$regions_check and $regions_check2 has $number of overlap!\n" if $number != 0;
			
		}
	}
	foreach my $genomic_check (@genomic_check) {
		foreach my $genomic_check2 (@genomic_check) {
			next if $genomic_check eq $genomic_check2;
			my ($number) = `bedtools intersect -u -a $genomic_check -b $genomic_check2 | wc -l` =~ /^(\d+)/;
			die "$genomic_check and $genomic_check2 has $number of overlap!\n" if $number != 0;
		}
	}
}
