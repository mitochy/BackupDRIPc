#!/usr/bin/perl
# This script check if the names of final merged gencode result covers all gene names from main gencode gtf

use strict; use warnings; use mitochy;

my ($input, $input2) = @ARGV;
die "usage: $0 <Main Exon Bed File> <Final Merged File from MergeGene.pl>

Main GTF File: hg19_gencode_exon.bed
Final Merged File: hg19_gencode_exon_merged.out

./CheckNames.pl hg19_gencode_exon.bed hg19_gencode_exon_merged.out

" unless @ARGV;

my ($folder, $filename) = mitochy::getFilename($input, "folder");

# 1. Get all names from main exon bed file
print "Processing main exon bed file $input\n";
my %names;
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $gene, $zero, $strand) = split("\t", $line);
	$names{$gene} = 0;
}
close $in;

# 2. Get all names from final merged file
print "Processing final merged file $input2\n";
open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $gene, $zero, $strand) = split("\t", $line);
	my @names = split(",", $gene);

	# Check the status of names in @names
	for (my $i = 0; $i < @names; $i++) {
		my $name = $names[$i];

		# If defined name, then value of name is 1
		if (defined($names{$name}) and $names{$name} =~ /^\d+$/) {
			$names{$name} = 1;
		}
		
		# If not defined name then create one with value -1
		elsif (not defined($names{$name})) {
			$names{$name} = -1;
		}
	}
}
close $in2;

# Report how many names with 0 and how many with -1
my $number_of_zero  = 0;
my $number_of_minus = 0;
foreach my $name (sort keys %names) {
	$number_of_zero  ++ if $names{$name} == 0;
	$number_of_minus ++ if $names{$name} == -1;
}

print "Gene no grouped: $number_of_zero\nGene not found in main gtf database: $number_of_minus\n";
if ($number_of_zero > 0) {
	open (my $out, ">", "$filename.zero") or die "Cannot write to $filename.zero: $!\n";
	foreach my $name (sort keys %names) {
		print $out "$name\n" if $names{$name} == 0;
	}
	close $out;
}
if ($number_of_minus > 0) {
	open (my $out, ">", "$filename.minus") or die "Cannot write to $filename.minus: $!\n";
	foreach my $name (sort keys %names) {
		print $out "$name\n" if $names{$name} == -1;
	}
	close $out;
}
