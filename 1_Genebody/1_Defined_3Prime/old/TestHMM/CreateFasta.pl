#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($bedFile, $dripcFile, $chromFile) = ("hg19_genecode_pos.bed", "dripc_pos.bed", "hg19_chrominfo.tsv");

print "Processing Chrominfo\n";
my %chr;
open (my $in, "<", $chromFile) or die "Cannot read from $chromFile: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $end) = split("\t", $line);
	$chr{$chr} = $end;
}
close $in;

print "Processing Bed\n";
my %bed = %{parse_bed($bedFile, "gene")};
print "Processing Dripc\n";
my %dripc = %{parse_bed($dripcFile, "dripc")};

print "Printing\n";
open (my $out, ">", "output.fa") or die "Cannot write to output.fa: $!\n";
foreach my $chr (keys %chr) {
	print $out ">$chr\n";
	my $end = $chr{$chr};
	for (my $i = 1; $i <= $end; $i++) {
		print "Done $i lines\n" if $i % 1000000 == 0;
		if (defined($dripc{$chr}{$i})) {
			if (defined($bed{$chr}{$i})) {
				print $out "B";
			}
			else {
				print $out "A";
			}
		}
		else {
			if (defined($bed{$chr}{$i})) {
				print $out "C";
			}
			else {
				print $out "N";
			}
		}
	}
	print $out "\n";
}


sub parse_bed {
	my ($input, $type) = @_;
	my %bed;
	open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
		if ($type eq "dripc") {
			for (my $i = $start; $i <= $end; $i++) {
				$bed{$chr}{$i} = 1;
			}
		}
		elsif ($type eq "gene") {
			my $pos = $strand eq "+" ? $end : $start;
			$bed{$chr}{$pos} = 1;
		}
	}
	
	close $in;
	return(\%bed);
}



__END__

N = not interesting
A = has dripc
B = has dripc and terminal
C = has terminal
