#!/usr/bin/perl

use strict; use warnings; use mitochy; use Cache::FileCache;

my ($folder, $stringent, $boolean) = @ARGV;
die "usage: $0 <Folder containing all enhancer> <stringency (0.5, 0.75)> <boolean to redo cache>\n" unless @ARGV >= 2;
my @input = <$folder/*.merged>;

my $cache = new Cache::FileCache;
$cache->set_cache_root("/data/mitochi/Work/Cache/");

my $respond;
if (defined($boolean)) {
	print "Are you sure you want to erase all enhancer cache? (\"y\" to continue)\n";
	$respond = <STDIN>;
}
my %data;
my @chr = qw(chr1 chr2 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY);

if (defined($respond) and $respond eq "y") {
	foreach my $chr (@chr) {
		print "Are you sure to delete ALL cache?\n";
		$respond = <STDIN>;
		$data{$chr}{0} = 1;
		$cache->set("ENHANCER\.$chr", \%data);
	}
	my $fileNumber = @input;
	for (my $i = 0; $i < @input; $i++) {
		my $input = $input[$i];
		my $currentNumber = $i + 1;
		print "Processing $currentNumber / $fileNumber $input\n";
		open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
	
		my $lastchr = "INIT";
		my $check = 0;
		while (my $line = <$in>) {
			chomp($line);
			next if $line =~ /#/;
			my ($chr, $start, $end) = split("\t", $line);
			if ($chr ne $lastchr) {
				print "\t$chr\n";
				if ($check == 1) {
					$cache->set("ENHANCER\.$lastchr", \%data);
				}
				%data = ();
				my $data = $cache->get("ENHANCER\.$chr");
				$check = defined($data) ? 1 : 0;
				$lastchr = $chr;
				%data = %{$data} if defined($data);
			}
			$lastchr = $chr;
			next if $check == 0;
			for (my $i = $start; $i < $end; $i++) {
				$data{$chr}{$i} ++;
			}
		}
		$cache->set("ENHANCER\.$lastchr", \%data) if $check == 1;
		%data = ();
		close $in;
	}
}
open (my $out, ">", "ALL_$stringent.BED") or die;
foreach my $chr (sort @chr) {
	print "Processing chr $chr\n";
	my $data = $cache->get("ENHANCER\.$chr");
	print "Undefined cache at chromosome $chr (did you do cache yet?)\n" and next if not defined($data);
	%data = %{$data};
	foreach my $pos (sort {$a <=> $b} keys %{$data{$chr}}) {
		my $end = $pos + 1;
		print $out "$chr\t$pos\t$end\n" if $data{$chr}{$pos} >= $stringent * @input;
	}
}

close $out;
