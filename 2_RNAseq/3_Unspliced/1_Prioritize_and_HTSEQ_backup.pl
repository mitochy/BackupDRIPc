#!/usr/bin/perl
# Necessary files:
# hg19Gencode19All.tsv
# hg19Gencode_exon.gtf
# hg19Gencode_intron.gtf
use strict; use warnings; use mitochy; use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my ($RNAsam) = @ARGV;
die "
Usage: $0 <RNA sam>
E.g. $0 NT2_rep1_uniq.sam

" unless @ARGV == 1 and defined($RNAsam);

my ($folder1, $RNAsamName) = mitochy::getFilename($RNAsam, "folder");
mkdir "1_Intron_Further_Process" if not -d "1_Intron_Further_Process";

# 1. Get RP\d+\- gene list
print BOLD RED "\nA. Processing Exon\n";
print GREEN "\t1. Hashing RP[digit]- genes (e.g. RP1-) from hg19Gencode19All.tsv\n\t\tThese are 2nd priority genes\n";
open (my $Gencodein, "<", "hg19Gencode19All.tsv") or die "Cannot read from hg19Gencode19All.tsv: $!\n";
my %good;
while (my $line = <$Gencodein>) {
	chomp($line);
	next if $line =~ /^\#/;
	$good{$line} = 0;
}
close $Gencodein;

# 2. Split exon to those legit and RP only
print GREEN "\t2. Creating exon1.gtf and exon0.gtf: These are made by prioritizing exon, where exon with RP[digit]+ will be second priority (exon0)
\tOutput:
\t\t$RNAsamName\_hg19Gencode_exon1.gtf
\t\t$RNAsamName\_hg19Gencode_exon0.gtf\n";
open (my $exonIn, "<", "hg19Gencode_exon.gtf") or die "Cannot read from hg19Gencode_exon.gtf: $!\n";
open (my $exonOut, ">", "$RNAsamName\_hg19Gencode_exon1.gtf") or die "Cannot write from $RNAsamName\_hg19Gencode_exon1.gtf: $!\n";
open (my $exonOut2, ">", "$RNAsamName\_hg19Gencode_exon0.gtf") or die "Cannot write from $RNAsamName\_hg19Gencode_exon0.gtf: $!\n";
my $group = 0;
while (my $line = <$exonIn>) {
	chomp($line);
	$group++;
	#print "Group $group\n";
	next if $line =~ /^\#/;
	my ($chr, $source, $type, $start, $end, $zero, $strand, $dot, $info) = split("\t", $line);
	my ($name) = $info =~ /gene_id \"(.+)\"; transcript/;

	# Check each name, if there is at least one not defined then we put it in higher priorty (exon1.gtf)
	my @names = split(",", $name);
	my $check = 0;
	foreach my $name (@names) {
		if (defined($good{$name}) and $good{$name} == 0) {
			#print "\tGroup $group\t$name\tDefined\n";
			print $exonOut "$line\n";
			$check = 1;
			last;
		}
	}

	# Else we put in lower priority (exon0.gtf)
	if ($check == 0) {
		print $exonOut2 "$line\n";
	}
}
close $exonIn;

# 3. Do HTseq count for exon1.gtf first
print GREEN "\t3. Running htseq-count with hg19Gencode_exon1.gtf then hg19Gencode_exon0.gtf\n\t\texon1.gtf is first priority, then exon0.gtf\n";
my @files = ("$RNAsamName\_hg19Gencode_exon1.gtf", "$RNAsamName\_hg19Gencode_exon0.gtf");
my $nofeature = $RNAsam;
open (my $exonOut3, ">", "$RNAsamName\_exon_mapping_report.txt") or die "Cannot write to $RNAsamName\_exon_mapping_report.txt: $!\n";
print $exonOut3 "Filename\tUnmap\tMap\n";
for (my $i = 0; $i < @files; $i++) {
	my ($folder2, $filename) = mitochy::getFilename($files[$i], "folder");
	
	# Check linecount of $filename. Next if it that has no gene in it
	#my ($linecount) = `wc -l $filename`;
	#($linecount) = $linecount =~ /^(\d+)/;
	#next if $linecount == 0;

	# Run HTseq-count
	my $command  = "htseq-count -q -s no -m intersection-nonempty -o $RNAsamName.reads.sam $nofeature $files[$i] > $filename.rnacount";
	my $command2 = "grep no_feature $RNAsamName.reads.sam | cut -f1-20 > $RNAsamName.reads_nofeature.sam";
	print "\t\t$i.1 $command\n";
	system($command) == 0 or die "Failed to run $command: $!\n";
	print "\t\t$i.2 $command2\n";
	system($command2) == 0 or die "Failed to run $command2: $!\n";
	#fix_count("$RNAsamName\_$filename.rnacount");
	my ($unmap) = `grep -c no_feature $RNAsamName.reads.sam`;
	my ($map)   = `grep -c -v no_feature $RNAsamName.reads.sam`;
	print $exonOut3 "$filename\t$unmap\t$map\n";
	$nofeature = "$RNAsamName.reads_nofeature.sam";
}

# 4. Combine exon rnacount into one
print GREEN "\t4. Combining exon rnacount into $RNAsamName\_hg19Gencode_exon.rnacount\n";
system("rm $RNAsamName\_hg19Gencode_exon.rnacount") if -e "$RNAsamName\_hg19Gencode_exon.rnacount";
print "\tcat $RNAsamName\_hg19Gencode_exon1.rnacount $RNAsamName\_hg19Gencode_exon0.rnacount > $RNAsamName\_hg19Gencode_exon.rnacount\n";
system("cat $RNAsamName\_hg19Gencode_exon1.rnacount $RNAsamName\_hg19Gencode_exon0.rnacount > $RNAsamName\_hg19Gencode_exon.rnacount");

# 5. Process intron. Get each gene count from final exon count above
print BOLD RED "\nB. Processing Intron\n";
print GREEN "\t1. Getting each gene count from $RNAsamName\_hg19Gencode_exon.rnacount\n";
my $RNAcount = "$RNAsamName\_hg19Gencode_exon.rnacount";
$nofeature = "$RNAsamName.reads_nofeature.sam";

my %rna;
open (my $in, "<", $RNAcount) or die "Cannot read from $RNAcount: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name, $count) = split("\t", $line);
	$rna{$name} = $count
}

close $in;


# 6. Prioritize intron based on the gene's count (highest to lowest):
# >10k count
# 1k-10k count
# 100-1k count
# 50-100 counts
# 10-50 counts
# 1-10 count(s)
# 0 count
# The breaks above is completely arbritary, I did it based on quartile distribution of gene exon counts then by eye approximate so that each group has at least 1-2k genes.
print GREEN "\t2. Splitting main intron gtf hg19Gencode_intron.gtf into 7 based on the intron's gene's count (higher count has higher priority)\n";
# Input file is hg19 Gencode intron.gtf

open (my $in2, "<", "hg19Gencode_intron.gtf");
# Open each output file:
open (my $out1, ">", "$RNAsamName\_hg19Gencode_intron10k.gtf") or die "Cannot write to $RNAsamName\_hg19Gencode_intron10k.gtf: $!\n";
open (my $out2, ">", "$RNAsamName\_hg19Gencode_intron1k.gtf")  or die "Cannot write to $RNAsamName\_hg19Gencode_intron1k.gtf: $!\n";
open (my $out3, ">", "$RNAsamName\_hg19Gencode_intron100.gtf") or die "Cannot write to $RNAsamName\_hg19Gencode_intron100.gtf: $!\n";
open (my $out4, ">", "$RNAsamName\_hg19Gencode_intron50.gtf")  or die "Cannot write to $RNAsamName\_hg19Gencode_intron50.gtf: $!\n";
open (my $out5, ">", "$RNAsamName\_hg19Gencode_intron10.gtf")  or die "Cannot write to $RNAsamName\_hg19Gencode_intron10.gtf: $!\n";
open (my $out6, ">", "$RNAsamName\_hg19Gencode_intron1.gtf")   or die "Cannot write to $RNAsamName\_hg19Gencode_intron1.gtf: $!\n";
open (my $out7, ">", "$RNAsamName\_hg19Gencode_intron0.gtf")   or die "Cannot write to $RNAsamName\_hg19Gencode_intron0.gtf: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	my ($chr, $source, $type, $start, $end, $zero, $strand, $dot, $info) = split("\t", $line);
	my ($name) = $info =~ /gene_id \"(.+)\_\d+\"; transcript/;

	# Bug check. If there is an intron name which doesn't have exon name then die	
	if (not defined($rna{$name})) {
		die "Undefined RNA of name $name\n";
	}

	# Prioritize based on count
	else {
		if ($rna{$name} == 0) {
			print $out7 "$line\n";
		}
		elsif ($rna{$name} < 10) {
			print $out6 "$line\n";
		}
		elsif ($rna{$name} < 50) {
			print $out5 "$line\n";
		}
		elsif ($rna{$name} < 100) {
			print $out4 "$line\n";
		}
		elsif ($rna{$name} < 1000) {
			print $out3 "$line\n";
		}
		elsif ($rna{$name} < 10000) {
			print $out2 "$line\n";
		}
		else {
			print $out1 "$line\n";
		}
	}
}
close $in2;

# 7. Do HTseq for each intron prioritized by 10k until 0. 
print GREEN "\t3. Running htseq-count with $RNAsamName\_hg19Gencode_intron10k.gtf until $RNAsamName\_hg19Gencode_intron0.gtf\n";
@files = ("$RNAsamName\_hg19Gencode_intron10k.gtf","$RNAsamName\_hg19Gencode_intron1k.gtf","$RNAsamName\_hg19Gencode_intron100.gtf","$RNAsamName\_hg19Gencode_intron50.gtf","$RNAsamName\_hg19Gencode_intron10.gtf","$RNAsamName\_hg19Gencode_intron1.gtf","$RNAsamName\_hg19Gencode_intron0.gtf");

open (my $out8, ">", "$RNAsamName\_intron_mapping_report.txt") or die "Cannot write to $RNAsamName\_intron_mapping_report.txt: $!\n";
print $out8 "Filename\tUnmap\tMap\n";
for (my $i = 0; $i < @files; $i++) {
	my ($folder3, $filename) = mitochy::getFilename($files[$i], "folder");
	my $command  = "htseq-count -q -s no -m intersection-nonempty -o $RNAsamName.reads.sam $nofeature $files[$i] > $filename.rnacount";
	my $command2 = "grep no_feature $RNAsamName.reads.sam | cut -f1-20 > $RNAsamName.reads_nofeature.sam";
	print "\t\t$i. $command\n";
	print "\t\t$i.2 $command2\n";
	system($command) == 0 or die "Failed to run $command: $!\n";
	system($command2) == 0 or die "Failed to run $command2: $!\n";
	#fix_count("$RNAsamName\_$filename.rnacount");
	my ($unmap) = `grep -c no_feature $RNAsamName.reads.sam`;
	my ($map)   = `grep -c -v no_feature $RNAsamName.reads.sam`;
	print $out8 "$filename\t$unmap\t$map\n";
	$nofeature = "$RNAsamName.reads_nofeature.sam";
}

# 8. Combine intron rnacount into one
print GREEN "\t4. Combining exon rnacount into $RNAsamName\_hg19Gencode_exon.rnacount\n";
system("rm $RNAsamName\_hg19Gencode_intron.rnacount") if -e "$RNAsamName\_hg19Gencode_intron.rnacount";
print "\tcat $RNAsamName\_hg19Gencode_intron*.rnacount > $RNAsamName\_hg19Gencode_intron.rnacount";
system("cat $RNAsamName\_hg19Gencode_intron*.rnacount > $RNAsamName\_hg19Gencode_intron.rnacount");

# 8. Moving all intermediate files into 1_Intron_Further_Process
print BOLD RED "\nC. Post process: Moving all intermediate files into 1_Intron_Further_Process\n";
system("mv $RNAsamName.reads.sam 1_Intron_Further_Process");
system("mv $RNAsamName.reads_nofeature.sam 1_Intron_Further_Process");
system("mv $RNAsamName\_intron_mapping_report.txt 1_Intron_Further_Process");
system("mv $RNAsamName\_exon_mapping_report.txt 1_Intron_Further_Process");
system("mv $RNAsamName\_*.rnacount 1_Intron_Further_Process");
system("mv $RNAsamName\_hg19Gencode_intron*.gtf 1_Intron_Further_Process");
system("mv $RNAsamName\_hg19Gencode_exon*.gtf 1_Intron_Further_Process");
system("mv 1_Intron_Further_Process/$RNAsamName\_hg19Gencode_exon.rnacount ./");
system("mv 1_Intron_Further_Process/$RNAsamName\_hg19Gencode_intron.rnacount ./");

sub fix_count {
        my ($input) = @_;
        system("perl -pi -e \'s\/^no_featu.\+\\n\$\/\/' $input");
        system("perl -pi -e \'s\/^ambiguou.\+\\n\$\/\/' $input");
        system("perl -pi -e \'s\/^too_low_.\+\\n\$\/\/' $input");
        system("perl -pi -e \'s\/^not_alig.\+\\n\$\/\/' $input");
        system("perl -pi -e \'s\/^alignmen.\+\\n\$\/\/' $input");
}

