#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox; use Cache::FileCache; use Getopt::Std;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use vars qw($opt_w $opt_i);
getopts("w:i:");
# - Map Wig files to each (average mapping -m)
# - Combine DRIP results replacing 0 value with DRIPc values
# - Create boxplot of signal

die "usage: $0 -w WIG -i <region Folder>\n" unless defined($opt_w) and defined($opt_i);
my ($wigFile, $regionFolder) = ($opt_w, $opt_i);

# Check if NT2 wig files cache exist
print BLUE "1. Checking NT2 Wig Cache\n";
my ($wigFolder, $wigFilename) = mitochy::getFilename($wigFile, "folder");
my $cache = new Cache::FileCache;
$cache->set_cache_root("/data/mitochi/Work/Cache/");
my $wigCache = $cache->get("$wigFilename\_BIG");
runbash("echo \"chr1\t500\t600\tTEST\t0\t\+\" > TEMP_MAPWIG_FILE.bed");
runbash("map_wig_to_bed_BIG.pl -w $wigFile -m TEMP_MAPWIG_FILE.bed") if not defined($wigCache);
runbash("rm *TEMP_MAPWIG_FILE*");

# Process DRIPc bed files from 4_CreateRegion into positive and negative files
print BLUE "2. Running MapWig and rename files\n";
runbash("run_script_in_paralel2.pl -v \"map_wig_to_bed_BIG.pl -w $wigFile -o $regionFolder -m FILENAME\" $regionFolder bed 12");
my @files = <$regionFolder/*.txt>;
foreach my $file (@files) {
	my ($folder, $fileName) = mitochy::getFilename($file, "folder");
	my ($newName) = $fileName;
	$newName =~ s/$wigFilename\_//;
	runbash("mv $file $folder/$newName.bed");
}

# Further process map_wig_to_bed.pl output
print BLUE "3. Post-processing\n";
@files = <$regionFolder/*.bed>;print "FILES = @files\n";
my @color;
my $Rdata;
my $Rdata2 = "df = rbind(";
foreach my $file (@files) {
	#runbash("cp $file $file.backup");
	my ($folder, $fileName) = mitochy::getFilename($file, "folder");
	my ($newName) = $fileName;
	$newName =~ s/$wigFilename\_//;
	my $color = getColor($file);
	push(@color, $color);
	runbash("sort -k1,1 -k2,2n $file | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > $file\_TEMP && mv $file\_TEMP $file");
	$Rdata .= "X$fileName.dat = data.frame(id=\"$fileName\", value=read.table(\"$file\")\$V5)\n";
	$Rdata2 .= "X$fileName.dat,";
}
$Rdata2 .= ")";
$Rdata2 =~ s/,\)/\)/;
my $color = R_toolbox::newRArray(\@color, "color");

print BLUE "4. Running RScript\n";
# Get distribution
my $Rscript = "

$Rdata
$Rdata2
$color
library(ggplot2)
pdf(\"$regionFolder/Signal.pdf\")
ggplot(df,aes(id,value)) + geom_boxplot(aes(fill=id),outlier.shape=NA) + 
#coord_cartesian(ylim=c(0,30)) +
scale_fill_manual(values=color) +
theme(axis.text.x=element_blank())
dev.off()
";

R_toolbox::execute_Rscript($Rscript);
runbash("echo \"track name=\'$regionFolder Annotation\' itemRgb=On\" > $regionFolder/CombinedRegion.BED");
runbash("cat $regionFolder/*.bed >> $regionFolder/CombinedRegion.BED");
#runbash("rm dump temp.*.R");

print RED "Output:

1. $regionFolder/*.bed: DRIP output files in $regionFolder
2. $regionFolder/CombinedRegion.BED: Combined output files for UCSC track
3. $regionFolder/Signal.pdf: Boxplot of DRIP values 

";
sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}

sub getColor {
	my ($input) = @_;
	return("rgb(233,162,60,maxColorValue=255)") if $input =~ /promoter/;
	return("rgb(255,182,80,maxColorValue=255)") if $input =~ /promoter_ext/;
	return("rgb(221,182,11,maxColorValue=255)") if $input =~ /terminal/;
	return("rgb(241,202,31,maxColorValue=255)") if $input =~ /terminal_ext/;
	return("rgb(161,178,147,maxColorValue=255)") if $input =~ /both/;
	return("rgb(181,198,167,maxColorValue=255)") if $input =~ /both_ext/;
	return("rgb(158,133,188,maxColorValue=255)") if $input =~ /genebody/;
	return("rgb(120,149,179,maxColorValue=255)") if $input =~ /intergenic/;
	return("rgb(0,0,0,maxColorValue=255)");


}
