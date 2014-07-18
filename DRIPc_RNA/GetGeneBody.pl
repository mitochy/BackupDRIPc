#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($DRIPfile, $input) = @ARGV;
die "usage: $0 <DRIP Peaks> <input>\n" unless @ARGV == 2;

my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in2, "<", $DRIPfile) or die;
my %DRIP;

while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /track/;
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	$DRIP{$chr}{$start} = $end;
}
close $in2;
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";

while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name, $chr, $strand, $start, $end) = split("\t", $line);
	
	# Filter short genes
	next if ($end - $start < 10000);
	
	# Get middle 30% of gene
	my $lengthThird = int(0.3*($end - $start));
	my ($Gstart, $Gend) = ($start + $lengthThird, $end - $lengthThird);
	
	# Find DRIP peak that is inside the gene
	foreach my $Dstart (sort keys %{$DRIP{$chr}}) {
		my $Dend = $DRIP{$chr}{$Dstart};
		
		if (within($Gstart, $Gend, $Dstart, $Dend) == 1) {
			print $out "$chr\t$Dstart\t$Dend\t$name\t0\t$strand\n";
			last;
		}
	}
}

close $in;
close $out;

sub within {
	my ($start1, $end1, $start2, $end2) = @_;
	if ($start2 >= $start1 and $end2 <= $end1) {
		#print "START1 $start1 <= START2 $start2,  END1 $end1 >= END2 $end2\n";
		return 1 
	}
	return 0;

}
