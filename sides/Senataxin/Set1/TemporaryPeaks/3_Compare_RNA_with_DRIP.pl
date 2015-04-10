#!/usr/bin/perl
# This script find DRIP for significantly higher/lower expressed genes
# Process ALLPEAK (contain all DRIP peak) and prints out DRIP data. E.g. whether RNA is correlated with higher/lower DRIP
use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <ALLPEAK.BED>\n" unless @ARGV;

# Process RNA
my $RNA = "/data/mitochi/Work/Project/DRIPc/sides/Senataxin/Set1/RNA.rpkm";
my %rna;
open (my $in1, "<", $RNA) or die "Cannot read from $RNA: $!\n";
while (my $line = <$in1>) {
        chomp($line);
        next if $line =~ /#/;
        my ($name, $gene, $value, $qval) = split("\t", $line);
        $rna{$gene}{value} = $value;
        $rna{$gene}{name}  = $name;
        $rna{$gene}{qval} = $qval;
}
close $in1;

# Process ALLPEAK.BED
open (my $in2, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $names, $fold, $strand, $C, $D, $E, $F, $meanCD, $meanEF, $pval, $chr2, $start2, $end2) = split("\t", $line);
	my @names = split(";", $names);
	foreach my $gene (@names) {
		if (defined($rna{$gene})) {
			push(@{$rna{$gene}{fold}}, $fold);
			push(@{$rna{$gene}{pval}}, $pval);
			push(@{$rna{$gene}{meanwt}}, $meanCD);
			push(@{$rna{$gene}{meanKD}}, $meanEF);
			push(@{$rna{$gene}{chr}}, $chr2);
			push(@{$rna{$gene}{start}}, $start2);
			push(@{$rna{$gene}{end}}, $end2);
		}
	}
	my $number = keys %rna;
}
close $in2;

# Print out RNA data and its DRIP data
open (my $out1, ">", "RNA_High.RNA") or die "Cannot write to RNA_High.out: $!\n";
open (my $out2, ">", "RNA_Low.RNA") or die "Cannot write to RNA_Low.out: $!\n";
print $out1 "\#Gene\tmRNA\tmRNA Fold-Change (siSETX/siCtrl)\tmRNA q-value\tPeak Chr\tPeak Start\tPeak End\tDRIP Log2 Fold Change\tDRIP Mean Scramble\tDRIP Mean KD\tDRIP P-value (anova one way)\n";
print $out2 "\#Gene\tmRNA\tmRNA Fold-Change (siSETX/siCtrl)\tmRNA q-value\tPeak Chr\tPeak Start\tPeak End\tDRIP Log2 Fold Change\tDRIP Mean Scramble\tDRIP Mean KD\tDRIP P-value (anova one way)\n";
foreach my $gene (sort {$rna{$a}{value} <=> $rna{$b}{value}} keys %rna) {

	# Variables
	my $value = $rna{$gene}{value};
	my $name  = $rna{$gene}{name};
	my $qval  = $rna{$gene}{qval};
	
	# If this gene doesn't have DRIP, return NAs
	if (not defined($rna{$gene}{fold})) {
		print $out1 "$name\t$gene\t$value\t$qval\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" if $value > 1;
		print $out2 "$name\t$gene\t$value\t$qval\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" if $value <= 1;
	}

	# If this gene has DRIP, return all the 1kb DRIP peaks
	else {
		for (my $i = 0; $i < @{$rna{$gene}{fold}}; $i++) {
			my $fold   = $rna{$gene}{fold}[$i];
			my $pval   = $rna{$gene}{pval}[$i];
			my $meanwt = $rna{$gene}{meanwt}[$i];
			my $meanKD = $rna{$gene}{meanKD}[$i];
			my $chr    = $rna{$gene}{chr}[$i];
			my $start  = $rna{$gene}{start}[$i];
			my $end    = $rna{$gene}{end}[$i];
			print $out1 "$name\t$gene\t$value\t$qval\t$chr\t$start\t$end\t$fold\t$meanwt\t$meanKD\t$pval\n" if $value > 1;
			print $out2 "$name\t$gene\t$value\t$qval\t$chr\t$start\t$end\t$fold\t$meanwt\t$meanKD\t$pval\n" if $value <= 1;
		}
	}
}
close $out1;
close $out2;

__END__
