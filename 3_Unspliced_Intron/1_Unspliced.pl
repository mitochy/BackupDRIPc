#!/usr/bin/perl
# This script takes NT2 rep1 and 2 for both exon and intron and define unspliced introns
use strict; use warnings; use mitochy;

#my ($exonCount, $intronCount) = @ARGV;
#die "Usage: $0 <exonCount> <intronCount>\n" unless @ARGV;
# Parse hg19 intron and exon to get length
print "0. Parsing hg19 intron and exon\n";
my %data;
open (my $exon_in, "<", "hg19Gencode_exon.gtf") or die "Cannot read from hg19Gencode_exon.gtf: $!\n";
while (my $line = <$exon_in>) {
	chomp($line);
	my ($chr, $source, $type, $start, $end, $zero, $strand, $dot, $name) = split("\t", $line);
	($name) = $name =~ /gene_id \"(.+)\"\; transcript/;
	die "Undefined gtf name at $line\n" unless defined($name);
	$data{exon}{$name} += $end - $start;
	$data{pos}{$name} = "$chr $start $end\t";
}
close $exon_in;
open (my $intron_in, "<", "hg19Gencode_intron.gtf") or die "Cannot read from hg19Gencode_intron.gtf: $!\n";
while (my $line = <$intron_in>) {
	chomp($line);
	my ($chr, $source, $type, $start, $end, $zero, $strand, $dot, $name) = split("\t", $line);
	($name) = $name =~ /gene_id \"(.+)\"\; transcript/;
	die "Undefined gtf name at $line\n" unless defined($name);
	$data{intron}{$name} = $end - $start;
	$data{pos}{$name} = "$chr $start $end";
	$data{line}{$name} = "$chr\t$start\t$end\t$name\t$zero\t$strand";
}
close $intron_in;

# Parse RNAcount from both reps
print "1. Parsing RNA count\n";
my %exon1 = %{parse_count("NT2_rep1_uniq_hg19Gencode_exon.rnacount")};
my %exon2 = %{parse_count("NT2_rep2_uniq_hg19Gencode_exon.rnacount")};
my %intron1 = %{parse_count("NT2_rep1_uniq_hg19Gencode_intron.rnacount")};
my %intron2 = %{parse_count("NT2_rep2_uniq_hg19Gencode_intron.rnacount")};


# Average RNAcount of exon from both reps with these modifications:
print "2. Averaging RNA count from both reps\n";
my %exon;
foreach my $gene (keys %exon1) {
	$exon{$gene}{count} += $exon1{$gene};
}
foreach my $gene (keys %exon2) {
	$exon{$gene}{count} += $exon2{$gene};
}
foreach my $gene (keys %exon) {
	die "Undefined exon name of $gene from reference gtf\n" unless defined($data{exon}{$gene});
	my $total_basepair = $data{exon}{$gene};
	$exon{$gene}{count} = $exon{$gene}{count}*100 / $total_basepair;
#	print "Exon $gene count = $exon{$gene}{count} (divided by $total_basepair)\n";
}
%exon1 = ();
%exon2 = ();

# Average RNAcount of intron from both reps with these modifications:
# - From H1 strand specific wig files: Normalize using table from H1 strand specific intron forward vs reverse polyA-tailed
my %norm;
open (my $normTableIn, "<", "H1_hg19Gencode_intronNORMTABLE.tsv") or die "Cannot read from H1_hg19Gencode_intronNORMTABLE.tsv: $!\n";
while (my $line = <$normTableIn>) {
	chomp($line);
	my ($gene, $factor) = split("\t", $line);
	$norm{$gene} = $factor;	
}
close $normTableIn;
# - From RNA-seq peak at 3' end: Change the count of introns which are within -1973 to 5228 of TTS to 0
my %badintron;
open (my $badIntronIn, "<", "hg19Gencode_BADINTRON.bed") or die "Cannot read from hg19Gencode_BADINTRON.bed: $!\n";
while (my $line = <$badIntronIn>) {
	chomp($line);
	my ($chr, $start, $end, $gene, $type, $strand) = split("\t", $line);
	$badintron{$gene} = 1;	
}
close $badIntronIn;
# Process Intron
my %intron;
foreach my $gene (keys %intron1) {
	# Don't use intron within 0-5228 of TTS
	next if defined($badintron{$gene});
	$intron{$gene}{count} += $intron1{$gene};
}
foreach my $gene (keys %intron2) {
	# Don't use intron within 0-5228 of TTS
	next if defined($badintron{$gene});
	$intron{$gene}{count} += $intron2{$gene};
}
foreach my $gene (keys %intron) {
	
	# Get normalization factor (or 1 if not defined.. since we're using H1 not all intron has expression in both H1 and NT2)
	my $norm = defined($norm{$gene}) ? $norm{$gene} : 1;
	# Process the intron	
	die "Undefined intron name of $gene from reference gtf\n" unless defined($data{intron}{$gene});
	my $total_basepair = $data{intron}{$gene};
	$intron{$gene}{count} = $intron{$gene}{count} * 100 / $total_basepair;
	$intron{$gene}{count} *= $norm;
}
%intron1 = ();
%intron2 = ();

# Get unspliced intron: count is more than 0.2 of the count its gene
print "3. Printing intron vs gene count\n";
my %final;
foreach my $gene (keys %intron) {
	my ($name) = $gene =~ /^(ENS.+\w+\.\d+)\_\d+/;
	my $pos = $data{pos}{$gene};
	my $line = $data{line}{$gene};
	#print "GENE = $name\n";
	die "Died at $gene\n" if (not defined($name));
	if (defined($exon{$name}{count})) {
		my $exon_count   = $exon{$name}{count};
		my $intron_count = $intron{$gene}{count};
		my $ratio = $exon_count == 0 ? $intron_count : $intron_count / $exon_count;
		$final{$gene}{ratio} = $ratio;
		$final{$gene}{line}  = "$gene\_$pos\t$intron_count\t$name\t$exon_count\t$ratio";
	}
	else {
		my $exon_count   = -1;
		my $intron_count = $intron{$gene}{count};
		my $ratio = -1;
		$final{$gene}{ratio} = $intron_count;
		$final{$gene}{line}  = "$gene\_$pos\t$intron_count\t$name\t$exon_count\t$ratio";
	}
}

open (my $out, ">", "NT2_unspliced_all.tsv") or die "Cannot write to NT2_unspliced_all.tsv: $!\n";
print $out "\#name_pos\tintroncount\tgene2\texoncount\tratio_intronexon\n";
foreach my $intron (sort {$final{$b}{ratio} <=> $final{$a}{ratio}} keys %final) {
	print $out "$final{$intron}{line}\n";
}
close $out;

open (my $out2, ">", "NT2_unspliced_filtered.tsv") or die "Cannot write to NT2_unspliced_filtered.tsv: $!\n";
print $out2 "\#chr\tstart\tend\tname\tratio_introncount_exoncount_gene2\tstrand\n";
foreach my $intron (sort {$final{$b}{ratio} <=> $final{$a}{ratio}} keys %final) {
	my $line = $data{line}{$intron};
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	next if $end - $start < 50;
	my ($gene, $intron_count, $gene2, $exon_count, $ratio) = split("\t", $final{$intron}{line});
	$ratio = int($ratio * 100)/100;
	if ($ratio >= 0.2 and $exon_count > 5 and $intron_count > 5) {
		my $norm = defined($norm{$gene}) ? $norm{$gene} : "NA";
		print $out2 "$chr\t$start\t$end\t$name\t$ratio\_$intron_count\_$exon_count\_$gene2\t$strand\n";
	}
	elsif ($ratio >= 0.2 and $intron_count > 5) {
		#print "$chr\t$start\t$end\t$name\t$ratio\_$intron_count\_$exon_count\_$gene2\t$strand\n";
	}
}
close $out2;

print "Output:
1. NT2_unspliced_all.tsv: All intron vs exon with ratio)
2. NT2_unspliced_filtered.tsv: Intron with ratio > 0.2, exon count and intron count > 0.005
";



sub parse_count {
	my ($input) = @_;
	my %count;
	open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
	
	while (my $line = <$in>) {
		chomp($line);
		next if $line =~ /#/;
		next if $line =~ /no_feature/ or $line =~ /ambiguous/ or $line =~ /too_low/ or $line =~ /not_align/ or $line =~ /alignment/;
		my ($genes, $count) = split("\t", $line);
		$count{$genes} = $count;
		#my @genes = split(",", $genes);
		#for (my $i = 0; $i < @genes; $i++) {
		#	if ($input =~ /intron/) {
		#		my ($suffix) = $genes[@genes-1] =~ /ENS.+(\_\d+)/;
		#	}
		#	my $gene = $genes[$i];
		#	$count{$gene} = $count;
		#}
	}
	close $in;
	return(\%count);
}

__END__




#R script
library(ggplot2)
df = read.table("NT2_unspliced_all.tsv",sep="\t")
colnames(df) = c("Name","Intron_count","Name2","Exon_count","Ratio")
df$Exon_count[which(df$Exon_count == -1)] = 0
good = df[which(df$Intron_count > 5 & df$Exon_count > 5 & df$Ratio >= 0.2),]
bad = df[which(!(df$Name %in% good$Name)),]
good$label="unspliced"
bad$label="spliced"
temp = rbind(good,bad)
pdf("NT2_unspliced_all.pdf")
ggplot(temp,aes(sqrt(Intron_count),sqrt(Exon_count))) + geom_point(aes(color=label,alpha=label,size=label)) + 
ylab("Exon read count per bp total exon length (sqrt)") + xlab("Intron read count / bp intron length (sqrt)") + 
scale_color_manual(values=c("black","red")) + scale_alpha_manual(values=c(0.2,1)) + scale_size_manual(values=c(0.5,0.5)) +
coord_cartesian(ylim=c(0,30))
dev.off()
