#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my @int = `cat /data/mitochi/Work/Project/DRIPc/bed/gencode_v19_annotation_InternalProm.id`;
my %int;
foreach my $line (@int) {
	chomp($line);
	my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
	$int{$name} = 1;
}
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

my %data;
#open (my $out1, ">", "$fileName1.NOINTERNAL") or die "Cannot write to $fileName1.NOINTERNAL: $!\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
	my ($orig) = $info =~ /ORIG=(E\w+\.\d+),/;
	($name) = $name =~ /^(ENS\w+\.\d+)_/ if $name =~ /_/;
	die "Died at line $line\n" if not defined($name);
	if (defined($int{$name})) {
		$data{intshuf} ++;
	}
	else {
		$data{extshuf} ++;
	}
	if (defined($int{$orig})) {
		$data{intorig}{$orig} = 1 
	}
	else {
		$data{extorig}{$orig} = 1 
	}
	#print $out1 "$line\n" if not defined($int{$name}) and not defined($int{$orig});
	#if (defined($int{$orig})) {
	#	$origInternal ++;
	#	$shufInternal ++;
	#	#$data{orig}{$orig} = -999;
#
#		#$data{shuf}{$line} = -999;
#		#push(@{$data{orig}{$orig}, $name);
	#}
	#if (defined($int{$name})) {# and not defined($int{$orig})) {
	#	$data{orig}{$orig} ++;
	#	$data{shuf}{$line} = -1;
	#}
	#elsif (not defined($data{orig}{$orig})) {
	#	$data{orig}{$orig} = 0;
	#	$data{orig}{$line} = 0;
	#}
}
close $in1;

my $countintshuf = $data{intshuf};
my $countintorig = (keys %{$data{intorig}});
my $countextshuf = $data{extshuf};
my $countextorig = (keys %{$data{extorig}});

printf "
Internal Promoter at orig: $countintorig (%.2f %%)
Internal Promoter at shuf: $countintshuf (%.2f %%)
External Promoter at orig: $countextorig
External Promoter at shuf: $countextshuf
", $countintorig *100 / ($countintorig + $countextorig), $countintshuf *100 / ($countextshuf + $countintshuf);

__END__
#close $out1;

my %count;
my $shufdie;
for (my $i = 1; $i <= 5; $i++) {
	$count{$i} = 0;
}
foreach my $orig (keys %{$data{orig}}) {
	my $count = $data{orig}{$orig};
	$count{$count} ++;
	$shufdie += $count;
}

$origexternal /= 5;
printf "
Internaled orig = $origInternal
Internaled shuf = $shufInternal
5 = $count{5}
4 = $count{4}
3 = $count{3}
2 = $count{2}
1 = $count{1}
0 = $count{0}
total Shuf die = $shufdie
";

