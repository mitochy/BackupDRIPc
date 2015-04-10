#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my @features = qw(promoter genebody terminal);
foreach my $feature (@features) {
	print "Processing $feature\n";
	my %data;
	my $orig = "dripc_$feature\_orig.peak";
	my $shuf = "dripc_$feature\_shuf.peak";
	%data = %{process_bed($orig, "orig")};
	%data = %{process_bed($shuf, "shuf", \%data)};
	my (@files) = <./Intersect/dripc_$feature\_*.peak\_*.bed>;
	foreach my $file (@files) {
		next if $file =~ /_rev_/;
		my ($chip) = $file =~ /.peak\_(.+).bed/; die "Died at $file not defined chip\n" unless defined($chip);
		print "\t$chip\n";
		%data = %{process_bed($orig, "$chip\_orig", \%data)};
		%data = %{process_bed($shuf, "$chip\_shuf", \%data)};
	}
	open (my $out, ">", "Intersect/ResultPeak_$feature.txt") or die;
	foreach my $name (keys %data) {
		print $out "ID";
		foreach my $chip (sort keys %{$data{$name}}) {
			next if $chip eq "init";
			print $out "\t$chip";
		}
		print $out "\n";
		last;
	}
	foreach my $name (keys %data) {
		print $out "$name";
		foreach my $chip (sort keys %{$data{$name}}) {
			next if $chip eq "init";
			print $out "\t$data{$name}{$chip}";
		}
		print $out "\n";
	}
	close $out;
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
		if ($type eq "orig") {
			$data{$name}{orig}{init} = 0;
		}
		elsif ($type eq "shuf") {
			$data{$name}{shuf}{init} = 0;
		}
		elsif ($type =~ /_orig/) {
			my ($chip) = $type =~ /^(.+)\_orig$/;
			$data{$name}{orig}{$chip} = 1;
		}
		elsif ($type =~ /_shuf/) {
			my ($chip) = $type =~ /^(.+)\_shuf$/;
			$data{$name}{shuf}{$chip} += 0.2;
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
