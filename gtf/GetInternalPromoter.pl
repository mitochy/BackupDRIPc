#!/usr/bin/perl

use strict; use warnings; use mitochy; use Cache::FileCache; use Getopt::Std;
use vars qw($opt_c);
getopts("c");

my ($input) = @ARGV;
die "
This script gets gene with earliest promoter (not necessary longest though..)

Usage: $0 <gencode_v19_annotation_All.gtf>\n

$0 gencode_v19_annotation_All.gtf

" unless @ARGV;

my ($folder, $filename) = mitochy::getFilename($input, "folder");
$filename =~ s/_All//i;
my $cache = new Cache::FileCache;
$cache -> set_cache_root("/data/mitochi/Work/Cache/");

parse_gtf($input);

#################
## Subroutines ##
#################

sub parse_gtf {
	my ($input) = @_;
	my ($folder, $filename) = mitochy::getFilename($input, "folder");
	$filename =~ s/_All//i;
	my %data;

	# 1. Check if cache $filename exist

	# 2a. If not exist or forced to re-cache
	open (my $outINTERNAL, ">", "../bed/$filename\_InternalProm.id") or die "Cannot write to ../bed/$filename\_InternalProm.id: $!\n";
	open (my $outEXTERNAL, ">", "../bed/$filename\_ExternalProm.id") or die "Cannot write to ../bed/$filename\_ExternalProm.id: $!\n";
	#if (not defined($data) or $opt_c) {
	open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
	my $lineCount = 0;
	my $geneCount = 0;
	while (my $line = <$in>) {
		chomp($line);
		next if $line =~ /#/;
		print "Done $lineCount\n" if $lineCount % 100000 == 0;
		$lineCount++;
		
		# Parse each line
		my ($chr, $source, $feat, $start, $end, $dot, $strand, $dot2, $name) = split("\t", $line);

		# Only take gene with KNOWN and NOVEL gene status
		next if $name !~ /transcript_status \"KNOWN\"\;/ and $name !~ /transcript_status \"NOVEL\"\;/;
		next if $name !~ /transcript_type \"protein_coding\"\;/ and $name !~ /transcript_type \"lincRNA\"\;/ and $name !~ /transcript_type \"antisense\"\;/;

		# For gtf, take all features but only protein_coding, lincRNA, or antisense
		my ($transcriptType) = $name =~ /transcript_type \"(\w+)\"\; transcript_status \"(\w+)\"; transcript/;

		# For bed take only transcripts
		next if $feat ne "transcript";
		my ($gid, $tid, $type) = $name =~ /^gene_id "(.+)"; transcript_id "(.+)"; gene_type.+ transcript_type "(.+)"; transcript_status/;
		die "\nDied at getting gene_id from name at line:\n$line\n\n" unless defined($gid) and defined($tid) and defined($type);

		$geneCount ++;
		$data{$gid}{$tid}{chr}    = $chr;
		$data{$gid}{$tid}{start}  = $start;
		$data{$gid}{$tid}{end}    = $end;
		$data{$gid}{$tid}{strand} = $strand;

	}
	close $in;
	print "Total Gene Count = $geneCount\n";

	foreach my $gid (keys %data) {
		my $bestname = "NA";
		my ($bestpos, $bestlength) = (0,0);
		foreach my $tid (keys %{$data{$gid}}) {
			my $strand = $data{$gid}{$tid}{strand};
			if (not defined($bestpos)) {
				$bestpos = $data{$gid}{$tid}{start} if $strand eq "+";
				$bestpos = $data{$gid}{$tid}{end} if $strand eq "-";
			}
			elsif ($strand eq "+" and $bestpos > $data{$gid}{$tid}{start} and $bestlength > $data{$gid}{$tid}{end} - $data{$gid}{$tid}{start}) {
				$bestname   = $tid;
				$bestpos    = $data{$gid}{$tid}{start};
				$bestlength = $data{$gid}{$tid}{end} - $data{$gid}{$tid}{start}
			}
			elsif ($strand eq "-" and $bestpos < $data{$gid}{$tid}{end} and $bestlength > $data{$gid}{$tid}{end} - $data{$gid}{$tid}{start}) {
				$bestname   = $tid;
				$bestpos    = $data{$gid}{$tid}{end};
				$bestlength = $data{$gid}{$tid}{end} - $data{$gid}{$tid}{start}
			}
		}
		foreach my $tid (keys %{$data{$gid}}) {
			my $chr    = $data{$gid}{$tid}{chr};
			my $start  = $data{$gid}{$tid}{start};
			my $end    = $data{$gid}{$tid}{end};
			my $strand = $data{$gid}{$tid}{strand};
			my $print  = "$chr\t$start\t$end\t$tid\t0\t$strand";
			print $outEXTERNAL "$print\n" and next if $tid eq $bestname;
			print $outEXTERNAL "$print\n" if $data{$gid}{$tid}{start} <= $bestpos + 2000 and $strand eq "+";
			print $outEXTERNAL "$print\n" if $data{$gid}{$tid}{end} >= $bestpos - 2000 and $strand eq "-";
			print $outINTERNAL "$print\n" if $data{$gid}{$tid}{start} > $bestpos + 2000 and $strand eq "+";
			print $outINTERNAL "$print\n" if $data{$gid}{$tid}{end} < $bestpos - 2000 and $strand eq "-";
		}
	}
	close $outINTERNAL;
	close $outEXTERNAL;
	system("cat ../bed/$filename\_InternalProm.id | sort -k1,1 -k2,2n > ../bed/$filename\_InternalProm.id_TEMP && mv ../bed/$filename\_InternalProm.id_TEMP ../bed/$filename\_InternalProm.id");
	system("cat ../bed/$filename\_ExternalProm.id | sort -k1,1 -k2,2n > ../bed/$filename\_ExternalProm.id_TEMP && mv ../bed/$filename\_ExternalProm.id_TEMP ../bed/$filename\_ExternalProm.id");
}

sub group_gene {
	my ($data) = @_;
	my %data = %{$data};
	foreach my $tid (keys %data) {
		my @exon = @{$data{$tid}};
	}
}





__END__
open (my $out, ">", "$name.out") or die "Cannot write to $name.out: $!\n";
close $out;
