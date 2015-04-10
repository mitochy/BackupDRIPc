#!/usr/bin/perl

use strict; use warnings; use mitochy; use R_toolbox; use Cache::FileCache;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

# - Process DRIPc bed files from 4_CreateRegion into positive and negative files
# - Map NT2 wig files to each (average mapping -m)
# - Combine DRIPc results replacing 0 value with DRIPc values
# - Create boxplot of signal

# Check if NT2 wig files cache exist
print BLUE "1. Checking NT2 Wig Cache\n";
my $wigPos = "../../wig/NT2_DRIPc_pos.wig";
my $wigNeg = "../../wig/NT2_DRIPc_neg.wig";
my $cache = new Cache::FileCache;
$cache->set_cache_root("/data/mitochi/Work/Cache/");
my $wigPosCache = $cache->get("NT2_DRIPc_pos.wig");
my $wigNegCache = $cache->get("NT2_DRIPc_neg.wig");
runbash("echo \"chr1\t500\t600\tTEST\t0\t\+\" > TEMP_MAPWIG_FILE.bed");
runbash("map_wig_to_bed.pl -w $wigPos -m -r /data/mitochi/Work/Cache/ TEMP_MAPWIG_FILE.bed") if not defined($wigPosCache);
runbash("map_wig_to_bed.pl -w $wigNeg -m -r /data/mitochi/Work/Cache/ TEMP_MAPWIG_FILE.bed") if not defined($wigNegCache);
runbash("rm *TEMP_MAPWIG_FILE*");

# Process DRIPc bed files from 4_CreateRegion into positive and negative files
print BLUE "2a. Creating Pos and Neg file from ../4_CreateRegion/DRIPc/\n";
runbash("run_script_in_paralel2.pl \"grep + FILENAME > FILENAME.pos\" ../4_CreateRegion/DRIPc/ bed 12");
runbash("run_script_in_paralel2.pl \"grep - FILENAME > FILENAME.neg\" ../4_CreateRegion/DRIPc/ bed 12");
mkdir "DRIPc_PosNeg" if not -d "DRIPc_PosNeg";
runbash("mv ../4_CreateRegion/DRIPc/dripc_*.pos ./DRIPc_PosNeg/ && mv ../4_CreateRegion/DRIPc/dripc_*.neg ./DRIPc_PosNeg/");

print BLUE "2b. Running MapWig\n";
runbash("run_script_in_paralel2.pl -v \"map_wig_to_bed.pl -w $wigPos -m -r /data/mitochi/Work/Cache/ FILENAME\" DRIPc_PosNeg pos 12");
runbash("run_script_in_paralel2.pl -v \"map_wig_to_bed.pl -w $wigNeg -m -r /data/mitochi/Work/Cache/ FILENAME\" DRIPc_PosNeg neg 12");

# Further process map_wig_to_bed.pl output
print BLUE "3. Post-processing\n";
runbash("cat NT2_DRIPc_pos_dripc_promoter.txt NT2_DRIPc_neg_dripc_promoter.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_promoter.bed");
runbash("cat NT2_DRIPc_pos_dripc_promoter_ext.txt NT2_DRIPc_neg_dripc_promoter_ext.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_promoter_ext.bed");
runbash("cat NT2_DRIPc_pos_dripc_terminal.txt NT2_DRIPc_neg_dripc_terminal.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_terminal.bed");
runbash("cat NT2_DRIPc_pos_dripc_terminal_ext.txt NT2_DRIPc_neg_dripc_terminal_ext.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_terminal_ext.bed");
runbash("cat NT2_DRIPc_pos_dripc_both.txt NT2_DRIPc_neg_dripc_both.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_both.bed");
runbash("cat NT2_DRIPc_pos_dripc_both_ext.txt NT2_DRIPc_neg_dripc_both_ext.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_both_ext.bed");
runbash("cat NT2_DRIPc_pos_dripc_antisense.txt NT2_DRIPc_neg_dripc_antisense.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_antisense.bed");
runbash("cat NT2_DRIPc_pos_dripc_antisense_other.txt NT2_DRIPc_neg_dripc_antisense_other.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_antisense_other.bed");
runbash("cat NT2_DRIPc_pos_dripc_antisense_ext.txt NT2_DRIPc_neg_dripc_antisense_ext.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_antisense_ext.bed");
runbash("cat NT2_DRIPc_pos_dripc_genebody.txt NT2_DRIPc_neg_dripc_genebody.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_genebody.bed");
runbash("cat NT2_DRIPc_pos_dripc_intergenic.txt NT2_DRIPc_neg_dripc_intergenic.txt | sort -k1,1 -k2,2n | perl -pi -e 's/^(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\\t(.+)\$/\$1\\t\$2\\t\$3\\t\$5\\t\$4\\t\$7\\t\$8\\t\$9\\t\$10/' > dripc_intergenic.bed");
mkdir "WigMappedResult" if not -d "WigMappedResult";
runbash("mv NT2_DRIPc_*.txt WigMappedResult");
print BLUE "4. Running RScript\n";
# Get distribution
my $Rscript = "

prom = rgb(233,162,60,maxColorValue=255)
promE = rgb(255,182,80,maxColorValue=255)
term = rgb(221,182,11,maxColorValue=255)
termE = rgb(241,202,31,maxColorValue=255)
both = rgb(161,178,147,maxColorValue=255)
bothE = rgb(181,198,167,maxColorValue=255)
gbod = rgb(158,133,188,maxColorValue=255)
anti = rgb(205, 144,163,maxColorValue=255)
antiE = rgb(225, 164,183,maxColorValue=255)
inte = rgb(120, 149, 179,maxColorValue=255)

df = data.frame(id=\"Promoter\",value=(read.table(\"dripc_promoter.bed\")\$V5))
df = rbind(df, data.frame(id=\"Promoter (extended)\",value=read.table(\"dripc_promoter_ext.bed\")\$V5))
df = rbind(df, data.frame(id=\"Terminal\",value=read.table(\"dripc_terminal.bed\")\$V5))
df = rbind(df, data.frame(id=\"Terminal (extended)\",value=read.table(\"dripc_terminal_ext.bed\")\$V5))
df = rbind(df, data.frame(id=\"Both\",value=read.table(\"dripc_both.bed\")\$V5))
df = rbind(df, data.frame(id=\"Both (extended)\",value=read.table(\"dripc_both_ext.bed\")\$V5))
df = rbind(df, data.frame(id=\"Genebody\",value=read.table(\"dripc_genebody.bed\")\$V5))
df = rbind(df, data.frame(id=\"Antisense\",value=read.table(\"dripc_antisense.bed\")\$V5))
df = rbind(df, data.frame(id=\"Antisense_Other\",value=read.table(\"dripc_antisense_other.bed\")\$V5))
df = rbind(df, data.frame(id=\"Antisense (extended)\",value=read.table(\"dripc_antisense_ext.bed\")\$V5))
df = rbind(df, data.frame(id=\"Intergenic\",value=read.table(\"dripc_intergenic.bed\")\$V5))

library(ggplot2)
pdf(\"FigSup_DRIPc_Signal.pdf\")
ggplot(df,aes(id,value)) + geom_boxplot(aes(fill=id),outlier.shape=NA) + coord_cartesian(ylim=c(0,30)) +
scale_fill_manual(values=c(prom,promE,term,termE,both,bothE,gbod,anti,anti,antiE,inte))
dev.off()
";

R_toolbox::execute_Rscript($Rscript);
runbash("cat dripc_*.bed > AllDRIPc.bed && mv dripc_*.bed AllDRIPc.bed ../");
runbash("rm dump temp.*.R");

sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}
