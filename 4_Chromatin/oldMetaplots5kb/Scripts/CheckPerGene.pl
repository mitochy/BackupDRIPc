#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($orig, $shuf) = @ARGV;
die "usage: $0 <orig.txt> <shuf.txt>\n" unless @ARGV == 2;

my %data;
my $length = 10000;
my @names;
print "Processing orig file $orig\n";
my ($folder1, $fileName1) = mitochy::getFilename($orig, "folder");
open (my $in1, "<", $orig) or die "Cannot read from $orig: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $val, $name, $val2, $strand, $info) = split("\t", $line);
	push(@names, $name);
	if ($val eq 0) {
		next;
	}
	my (@val) = split("_", $val);
	foreach my $val (@val) {
		my ($pos, $value) = split(",", $val);
		push(@{$data{$name}}, $value);
	}
	@{$data{$name}} = sort {$b <=> $a} @{$data{$name}};
}
close $in1;

print "Processing shuf file $shuf\n";
my %result;
my $longest = 0;
my ($folder2, $fileName2) = mitochy::getFilename($shuf, "folder");
open (my $in2, "<", $shuf) or die "Cannot read from $shuf: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $val, $name, $val2, $strand, $info) = split("\t", $line);
	my ($origName) = $info =~ /ORIG=(ENST\w+\.\d+),chr/;
	die "Died at line $line\n" if not defined($origName);
	if ($val eq "0") {
		if (not defined($data{$origName})) {
			$result{$origName}{type} = 0;
		}
		else {
			my $iter = 0;
			foreach my $value (sort {$b <=> $a} @{$data{$origName}}) {
				push(@{$result{$origName}{val}[$iter]}, $value);
				$longest = @{$data{$origName}} if @{$data{$origName}} > $longest;
				$result{$origName}{type} = 1 if not defined($result{$origName}{type}) or $result{$origName}{type} == 1;
				$iter ++;
			}
		}
		next;
	}

	my (@val) = split("_", $val);
	my @temp;
	foreach my $val (@val) {
		my ($pos, $value) = split(",", $val);
		push(@temp, $value);
	}
	my $iter = 0;
	$longest = @temp if $longest < @temp;
	foreach my $temp (sort {$b <=> $a} @temp) {
		if (not defined($data{$origName})) {
			push(@{$result{$origName}{val}[$iter]}, 0 - $temp);
			$result{$origName}{type} = 2 if not defined($result{$origName}{type}) or $result{$origName}{type} <= 2;
		#	print "ITER $iter NAME $origName TYPE $result{$origName}{type} TEMP $temp\n";
		}
		else {
			if (defined($data{$origName}[$iter])) {
				push(@{$result{$origName}{val}[$iter]}, $data{$origName}[$iter] - $temp);
				$result{$origName}{type} = 3 if not defined($result{$origName}{type}) or $result{$origName}{type} <= 3;

			}
			else {
				push(@{$result{$origName}{val}[$iter]}, 0 - $temp);
				$result{$origName}{type} = 3 if not defined($result{$origName}{type}) or $result{$origName}{type} <= 3;

			}
		}
		$iter ++;
	}
}
close $in2;

print "Calculating and print\n";
open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
foreach my $name (sort keys %result) {
	if ($result{$name}{type} == 0) {
		print $out1 "$name\_UNDEF";
		for (my $i = 0; $i < $longest; $i++) {
			print $out1 "\t0";
		}
		print $out1 "\n";
	}
	#elsif ($result{$name}{type} == 1) {
	#	print $out1 "$name\_ORIG";
#
#		for (my $i = 0; $i < $longest; $i++) {
#			print $out1 "\t$result{$name}{val}[$i]" if defined($result{$name}{val}[$i]);
#			print $out1 "\t0" if not defined($result{$name}{val}[$i]);
#		}
#		print $out1 "\n";
#	}
	#else ($result{$name}{type} == 3) {
	#	for (my $i = 0; $i < $longest; $i++) {
	#		print $out1 "\t$result{$name}{val}[$i]" if defined($result{$name}{val}[$i];
	#		print $out1 "\t0" if not defined($result{$name}{val}[$i]);
	#	}
	#	print $out1 "\n";
	#}
	else {
		print $out1 "$name\_ORIG" if $result{$name}{type} == 1;
		print $out1 "$name\_SHUF" if $result{$name}{type} == 2;
		print $out1 "$name\_BOTH" if $result{$name}{type} == 3;
		#print $out1 "$name\_BOTHORIGMIS" if $result{$name}{type} == 3;

		for (my $i = 0; $i < $longest; $i++) {
			print $out1 "\t0" and next if not defined($result{$name}{val}[$i]);
			my $median = median(@{$result{$name}{val}[$i]});
			print $out1 "\t$median";
		}
		print $out1 "\n";
	}
}
close $out1;


sub median {
	my @arr = @_;
	@arr = sort(@arr);
	return($arr[int(@arr/2)]);
}
