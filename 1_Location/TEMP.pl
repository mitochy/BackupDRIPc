#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $name) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";

while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $names, $strand) = split("\t", $line);
	my @names = split(";", $names);
	my @newname;
	foreach my $name (@names) {
		push(@newname, $name) if not grep(/^$name$/, @newname);
	}
	$names = join(";", @newname);
	print $out "$chr\t$start\t$end\t$names\t0\t$strand\n";
}

close $in;
close $out;
