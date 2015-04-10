#!/usr/bin/perl

use strict; use warnings; use mitochy; use FAlite;

my ($input, $RE) = @ARGV;
die "usage: $0 <input> <RE>\n" unless @ARGV == 2;

my ($folder, $name) = mitochy::getFilename($input, "folder");
$RE = uc($RE);
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
my $fasta = new FAlite($in);
my $lengthChunk = 1000;
my $lengthRE = length($RE);
open (my $out, ">", "$name\_$RE.bed") or die "Cannot write to $name\_$RE.bed: $!\n";
while (my $entry = $fasta->nextEntry()) {
	my ($def, $seq) = ($entry->def, $entry->seq);
	$def =~ s/>//;
	print "$input $RE processing chr $def\n";
	#my $sequence = $seq;
	my $start = 0;
	#while ($sequence =~ /$RE/ig) {

	# 100500 will be 101
	for (my $i = 0; $i <= length($seq); $i+= $lengthChunk) {
		my $lengthSubstr = $lengthChunk + $lengthRE > length($seq) ? length($seq) : $lengthChunk + $lengthRE;
		my $sequence = substr($seq, $i, $lengthSubstr);
		while ($sequence =~ /$RE/ig) {
			my $length = length($`)+1;
			my $currstart = $i + $length;
			my $currend   = $i + $length + length($RE) - 1;
			print $out "$def\t$currstart\t$currend\n";
		}
	}
	my $rev = uc(mitochy::revcomp($RE));
	next if ($RE eq $rev);
	for (my $i = 0; $i <= length($seq); $i+= $lengthChunk) {
		my $lengthSubstr = $lengthChunk + $lengthRE > length($seq) ? length($seq) : $lengthChunk + $lengthRE;
		my $sequence = substr($seq, $i, $lengthSubstr);
		while ($sequence =~ /$rev/ig) {
			my $length = length($`)+1;
			my $currstart = $i + $length;
			my $currend   = $i + $length + length($rev) - 1;
			print $out "$def\t$currstart\t$currend\n";
		}
	}
}

close $in;
close $out;


