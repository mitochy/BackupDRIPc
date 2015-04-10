#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my @int = `cat /data/mitochi/Work/Project/DRIPc/gtf/gencode_v19_annotation_InternalProm.id`;
my %int;
foreach my $int (@int) {
	chomp($int);
	$int{$int} = 1;
}
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my $count = 0;
open (my $in, "<", $input1) or die;
while (my $line = <$in>) {
	chomp($line);
	my ($name) = $line =~ /^(\w+\.\d+)/;
	$count++ if defined($int{$name});
}
close $in;

print "$count\n";
