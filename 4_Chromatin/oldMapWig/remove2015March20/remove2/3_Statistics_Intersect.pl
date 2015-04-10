#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1, $input2) = @ARGV;
die "usage: $0 <input.orig> <histone.bed>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
my ($folder2, $fileName2) = mitochy::getFilename($input2, "folder");
my ($shuffled) = "MappedShuffled/$fileName1.shuffled";
my $total1 = `wc -l $input1`;
my $intersect1 = `bedtools intersect -u -a $input1 -b $input2 | wc -l`;
my $intersect2 = `bedtools intersect -u -a $shuffled -b $input2`;
($total1) = $total1 =~ /^(\d+)/;
my ($result1) = $intersect1 =~ /^(\d+)/;
my %result2;
my @intersect2 = split("\n", $intersect2);
foreach my $line (@intersect2) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name) = split("\t", $line);
	push(@{$result2{$name}}, 1);
	
}

foreach my $name (keys %result2) {
	for (my $i = @{$result2{$name}}; $i < 125; $i++) {
		push(@{$result2{$name}}, 0);
	}
}
foreach my $name (keys %result2) {
	@{$result2{$name}} = shuffle(@{$result2{$name}});
}

my @count;
foreach my $name (keys %result2) {
	for (my $i = 0; $i < @{$result2{$name}}; $i++) {
		$count[$i] += $result2{$name}[$i];
	}
}
my $totalShuf = keys %result2;
for (my $i = 0; $i < @count; $i++) {
	print "$i\t$count[$i] / $totalShuf\n";
}
print "Orig: $result1 / $total1\n";
#foreach my $name (keys %result2) {
#	my $count = @{$result2{$name}};
#	print "$name\t$count\n";
#}


sub shuffle {
	my (@arr) = @_;
	for (my $i = 0; $i < 10000; $i++) {
		my $rand1 = int(rand(@arr));
		my $rand2 = int(rand(@arr));
		my $val1 = $arr[$rand1];
		my $val2 = $arr[$rand2];
		$arr[$rand1] = $val2;
		$arr[$rand2] = $val1;
	}
	return(@arr);
}
