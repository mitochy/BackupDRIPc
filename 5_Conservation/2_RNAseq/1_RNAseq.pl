#!/usr/bin/perl
# Necessary files:
use strict; use warnings; use mitochy; use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my ($RNAsam) = @ARGV;
die "
Usage: $0 <RNA sam>
E.g. $0 NT2_rep1_uniq.sam

" unless @ARGV == 1 and defined($RNAsam);

my ($folder1, $RNAsamName) = mitochy::getFilename($RNAsam, "folder");

my $proteinGTF   = "/data/mitochi/Work/Project/DRIPc/5_Conservation/1_RNAseq/mm10_gencode_vM2_exon_proteincoding.gtf";
my $noncodingGTF = "/data/mitochi/Work/Project/DRIPc/5_Conservation/1_RNAseq/mm10_gencode_vM2_exon_noncoding.gtf";
=comment
# 1. Do HTseq count for protein coding.gtf first
print GREEN "\t1. Running htseq-count with $proteinGTF then $noncodingGTF\n";
my @files = ($proteinGTF, $noncodingGTF);
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
	my $command  = "htseq-count -i transcript_id -q -s no -m intersection-nonempty -o $RNAsamName.reads.sam $nofeature $files[$i] > $RNAsamName\_$filename.rnacount";
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
=cut
# 4. Combine exon rnacount into one
print GREEN "\t2. Combining exon rnacount into $RNAsamName\_mm10_gencode_vM2_exon.rnacount\n";
system("rm $RNAsamName\_mm10_gencode_vM2_exon.rnacount") if -e "$RNAsamName\_mm10_gencode_vM2_exon.rnacount";
print "\tcat $RNAsamName\_mm10_gencode_vM2_exon_proteincoding.rnacount $RNAsamName\_mm10_gencode_vM2_exon_noncoding.rnacount > $RNAsamName\_exon.rnacount\n";
system("cat $RNAsamName\_mm10_gencode_vM2_exon_proteincoding.rnacount $RNAsamName\_mm10_gencode_vM2_exon_noncoding.rnacount > $RNAsamName\_mm10_gencode_vM2_exon.rnacount");

fix_count("$RNAsamName\_mm10_gencode_vM2_exon.rnacount");

sub fix_count {
        my ($input) = @_;
        system("perl -pi -e \'s\/^.+no_featu.\+\\n\$\/\/g' $input");
        system("perl -pi -e \'s\/^.+ambiguou.\+\\n\$\/\/g' $input");
        system("perl -pi -e \'s\/^.+too_low_.\+\\n\$\/\/g' $input");
        system("perl -pi -e \'s\/^.+not_alig.\+\\n\$\/\/g' $input");
        system("perl -pi -e \'s\/^.+alignmen.\+\\n\$\/\/g' $input");
}

