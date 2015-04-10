#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
#open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";

my ($chr, $span);
my $lastchr = "INIT";
my $total = 0;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	if ($line =~ /variable/) {
		($chr, $span) = $line =~ /chrom=(.+) span=(\d+)/;
		print "Processing chr $chr\n" if $chr ne $lastchr;
		$lastchr = $chr if $chr ne $lastchr;
	}
	else {
		#next if $chr ne "chr12";
		my ($pos, $val) = split("\t", $line);
		next if $chr =~ /chrM/;
		next if $val < 4 and $input =~ /E14/;
		next if $val < 8 and $input =~ /3T3/;
		next if $val < 4 and $input =~ /Fibro/;
		next if $val < 4 and $input =~ /NT2/;
		$total += $val*$span;
	}
	#last if $chr eq "chr13";
}

close $in;
print "$input: TOTAL = $total\n";
__END__
my $check = 0;
my $header;
open (my $in2, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	if ($line =~ /variable/) {
		($chr, $span) = $line =~ /chrom=(.+) span=(\d+)/;
		$header = $line;
		$check = 1;
	}
	else {
		my ($pos, $val) = split("\t", $line);
		next if $val < 5 and $input =~ /3T3/;
		next if $val < 4 and $input =~ /E14/;
		$val *= 1115037092/1388210226 if $input =~ /3T3/;
		if ($check == 1) {
			print $out "$header\n";
			$check = 0;
		}
		print $out "$pos\t$val\n";
	}
}

close $in2;

