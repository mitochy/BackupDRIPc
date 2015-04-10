#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my @features = qw(promoter terminal genebody);
for (my $I = 0; $I < @features; $I++) {
	my $feature = $features[$I];
	print "Processing $feature\n";
	my %data;
	my $orig = "dripc_$feature\_orig.peak";
	my $shuf = "dripc_$feature\_shuf.peak";
	%data = %{process_bed($orig, "orig")};
	%data = %{process_bed($shuf, "shuf", \%data)};
	my (@files) = <./Intersect/dripc_$feature\_*.peak\_*.bed>;
	my @chip;
	foreach my $file (@files) {
		next if $file =~ /_rev_/;
		next if $file =~ /shuf/;
		my $chiporig = $file;
		my $chipshuf = $chiporig; $chipshuf =~ s/orig/shuf/;
		my ($chip) = $file =~ /.peak\_(.+).bed/; die "Died at $file not defined chip\n" unless defined($chip);
		push(@chip, $chip);
		print "\t$chip\n";
		%data = %{process_bed($chiporig, "$chip\_orig", \%data)};
		%data = %{process_bed($chipshuf, "$chip\_shuf", \%data)};
	}
	open (my $out1, ">", "Intersect/ResultPeak_$feature\_orig.txt") or die;
	open (my $out2, ">", "Intersect/ResultPeak_$feature\_shuf.txt") or die;
	foreach my $name (keys %data) {
		print $out1 "ID";
		print $out2 "ID";
		foreach my $chip (@chip) {
			next if $chip eq "init";
			print $out1 "\t$chip";
		}
		print $out1 "\n";
		foreach my $chip (@chip) {
			next if $chip eq "init";
			print $out2 "\t$chip";
		}
		print $out2 "\n";
		last;
	}
	foreach my $name (keys %data) {
		print $out1 "$name";
		print $out2 "$name";
		my $rand = rand() < 0.5 ? -1 : 1;
		foreach my $chip (sort @chip) {
			next if $chip eq "init";
			my $intersect = defined($data{$name}{orig}{$chip}) ? $data{$name}{orig}{$chip} : 0 + rand()*0.0001 *$rand;
			print $out1 "\t$intersect";
		}
		print $out1 "\n";
		foreach my $chip (@chip) {
			next if $chip eq "init";
			my $intersect = defined($data{$name}{shuf}{$chip}) ? $data{$name}{shuf}{$chip} : 0 + rand()*0.0001 * $rand;
			print $out2 "\t$intersect";
		}
		print $out2 "\n";
	}
	close $out1;
	close $out2;
}


sub process_bed {
	my ($input1, $type, $data) = @_;
	my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
	my %data;
	%data = %{$data} if defined($data);
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
		($name) = $info =~ /DRIP=(.+);ORIG/; die "Died at $input1 undefined name\n" unless defined($name);
		my $rand = rand() < 0.5 ? -1 : 1;
		if ($type eq "orig") {
			$data{$name}{orig}{init} = 0 + rand()*0.000001 * $rand;
		}
		elsif ($type eq "shuf") {
			$data{$name}{shuf}{init} = 0 + rand()*0.000001 * $rand;
		}
		elsif ($type =~ /_orig/) {
			my ($chip) = $type =~ /^(.+)\_orig$/;
			$data{$name}{orig}{$chip} = 1 + rand()*0.00001 * $rand;
		}
		elsif ($type =~ /_shuf/) {
			my ($chip) = $type =~ /^(.+)\_shuf$/;
			$data{$name}{shuf}{$chip} += 1/125;
		}
		else {
			die;
		}
		
	}
	close $in1;
	return(\%data);
}

#open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
#close $out1;
