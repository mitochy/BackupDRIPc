#!/usr/bin/perl

use strict; use warnings; use mitochy;

# Get wig value for each rep for each pos and neg
# Parse intron from both pos and neg from same rep, count value for each, and calculate ratio of intron (forward) / intron (reverse)
# Calculate average ratio from rep1 and rep2 and print out average in bed format
# Exon will be used as control. Parse exon pos and neg from same rep, count average value for each gene from all exon with same strand
# Ignore exon with diff strandedness for now, but this will be used as a control. Exon with diff strandedness shouldn't have bigger expression

# Input files.
my $intron_pos1 = "H1_pos_rep1_hg19Gencode_intron.txt";
my $intron_pos2 = "H1_pos_rep2_hg19Gencode_intron.txt";
my $exon_pos1   = "H1_pos_rep1_hg19Gencode_exon.txt";
my $exon_pos2   = "H1_pos_rep2_hg19Gencode_exon.txt";
my $intron_neg1 = "H1_neg_rep1_hg19Gencode_intron.txt";
my $intron_neg2 = "H1_neg_rep2_hg19Gencode_intron.txt";
my $exon_neg1   = "H1_neg_rep1_hg19Gencode_exon.txt";
my $exon_neg2   = "H1_neg_rep2_hg19Gencode_exon.txt";

# Parse input files and calculate ratio of forward / reverse
my %intron1 = %{parse($intron_pos1, $intron_neg1)};
my %intron2 = %{parse($intron_pos2, $intron_neg2)};
#my %exon1   = %{parse($exon_pos1, $exon_neg1)};
#my %exon2   = %{parse($exon_pos2, $exon_neg2)};

# Output the intron
my %intron;
foreach my $gene (keys %intron1) {
	$intron{$gene} += $intron1{$gene} / 2;
}
foreach my $gene (keys %intron2) {
	$intron{$gene} += $intron2{$gene} / 2;
}
open (my $out, ">", "H1_hg19Gencode_intronNORMTABLE.tsv") or die "Cannot write to H1_hg19Gencode_intronNORMTABLE.tsv: $!\n";
foreach my $gene (keys %intron) {
	print $out "$gene\t$intron{$gene}\n";
}
close $out;

print "Output: H1_hg19Gencode_intronNORMTABLE.tsv\n";
# Output the exon
#my %exon;
#foreach my $gene (keys %exon1) {
#	$exon{$gene} += %exon1{$gene} / 2;
#}
#foreach my $gene (keys %exon2) {
#	$exon{$gene} += %exon2{$gene} / 2;
#}
#open (my $out, ">", "H1_hg19Gencode_exon.tsv") or die "Cannot write to H1_hg19Gencode_exon.tsv: $!\n";
#foreach my $gene (keys %exon) {
#	print $out "$gene\t$exon{$gene}\n";
#}
#close $out;


sub parse {
	my ($input1, $input2, $type) = @_;
	my $strand1 = $input1 =~ /pos/ ? "+" : "-";
	my $strand2 = $input2 =~ /pos/ ? "+" : "-";

	my %bed;
	#my $linecount = 0;
	## Parse input 1
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $val, $gene, $zero, $strand) = split("\t", $line);
		#print "1. $input1: $gene with strand $strand: ";
		if ($strand ne $strand1) {
			$bed{$gene}{chr}     = $chr;
			$bed{$gene}{start}   = $start;
			$bed{$gene}{end}     = $end;
			$bed{$gene}{strand}  = $strand;
			$bed{$gene}{reverse} += $val;
			$bed{$gene}{revcount}   ++;
			#print "reverse\tvalue $bed{$gene}{reverse} plus $val\tcount $bed{$gene}{revcount}\n";
		}
		else {
			$bed{$gene}{chr}     = $chr;
			$bed{$gene}{start}   = $start;
			$bed{$gene}{end}     = $end;
			$bed{$gene}{strand}  = $strand;
			$bed{$gene}{forward} += $val;
			$bed{$gene}{forcount}   ++;
			#print "forward\tvalue $bed{$gene}{forward} plus $val\tcount $bed{$gene}{forcount}\n";
		}
		#$linecount++;
		#next if $linecount < 242776;
		#next if $linecount < 287251;
		#last if $linecount > 10;
	}
	close $in1;

	## Parse input 2
	#$linecount = 0;
	open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
	while (my $line = <$in2>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $val, $gene, $zero, $strand) = split("\t", $line);
		#print "2. $input2: $gene with strand $strand: ";
		if ($strand ne $strand2) {
			$bed{$gene}{chr}     = $chr;
			$bed{$gene}{start}   = $start;
			$bed{$gene}{end}     = $end;
			$bed{$gene}{strand}  = $strand;
			$bed{$gene}{reverse} += $val;
			$bed{$gene}{revcount}   ++;
			#print "reverse\tvalue $bed{$gene}{reverse} plus $val\tcount $bed{$gene}{revcount}\n";
		}
		else {
			$bed{$gene}{chr}     = $chr;
			$bed{$gene}{start}   = $start;
			$bed{$gene}{end}     = $end;
			$bed{$gene}{strand}  = $strand;
			$bed{$gene}{forward} += $val;
			$bed{$gene}{forcount}   ++;
			#print "forward\tvalue $bed{$gene}{forward} plus $val\tcount $bed{$gene}{forcount}\n";
		}
		#$linecount++;
		#last if $linecount > 10;
	}
	close $in2;

	#print "\n\nRESULT\n\n";
	my %res;
	foreach my $gene (keys %bed) {
		## If forward or reverse not defined then assign 0
		$bed{$gene}{forward}  = 0 if not defined($bed{$gene}{forward});
		$bed{$gene}{reverse}  = 0 if not defined($bed{$gene}{reverse});
		$bed{$gene}{forcount} = 1 if not defined($bed{$gene}{forcount});
		$bed{$gene}{revcount} = 1 if not defined($bed{$gene}{revcount});
	
		## If forward count or reverse count is 0 then assi

		## If both forward and reverse is zero, then value is zero
		if ($bed{$gene}{forward} == 0 and $bed{$gene}{reverse} == 0) {
			$bed{$gene}{val} = 0;
			#print "$gene: RATIO: $bed{$gene}{val}\n";
		}
	
		## Otherwise, value is the ratio fo forward divided by sum of both
		else {
			#print "$gene: FORWARD $bed{$gene}{forward} / $bed{$gene}{forcount} = ";
			$bed{$gene}{forward} = $bed{$gene}{forward} / $bed{$gene}{forcount};
			#print "$bed{$gene}{forward}, REVERSE $bed{$gene}{reverse} / $bed{$gene}{revcount} = ";
			$bed{$gene}{reverse} = $bed{$gene}{reverse} / $bed{$gene}{revcount};
			#print "$bed{$gene}{reverse}, RATIO ";
			#print "$gene at $input1, where $bed{$gene}{forward} + $bed{$gene}{reverse}\n" if $bed{$gene}{forward} + $bed{$gene}{reverse} == 0;
			$bed{$gene}{val} = $bed{$gene}{forward} / ($bed{$gene}{forward} + $bed{$gene}{reverse});
			$res{$gene} = int(1000000 * $bed{$gene}{val}) / 1000000;
			#print "$bed{$gene}{forward} / ($bed{$gene}{forward} + $bed{$gene}{reverse}) = $bed{$gene}{val}\n";
		}
	}
	return(\%res);
}

__END__
