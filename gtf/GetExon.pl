#!/usr/bin/perl

use strict; use warnings; use mitochy; use Cache::FileCache; use Getopt::Std;
use vars qw($opt_c);
getopts("c");

my ($input) = @ARGV;
die "usage: $0 <gencode_v19_annotation_ALL.gtf>\n" unless @ARGV;

my ($folder, $filename) = mitochy::getFilename($input, "folder");
$filename =~ s/_All//i;
my $cache = new Cache::FileCache;
$cache -> set_cache_root("/data/mitochi/Work/Cache/");

# 1. Parsing
#my %data = %{parse_gtf($input)};
parse_gtf($input);
# 2. Print out genes in a BED file in /data/mitochi/Work/Project/DRIPc/bed/
#my $output = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
#open (my $outBed, ">", $output) or die "Cannot write to $output: $!\n";
#foreach my $type (sort {(keys %{$data{by_type}{$b}}) <=> (keys %{$data{by_type}{$a}})} keys %{$data{by_type}}) {

	# Only take protein coding, lincRNA, and antisense genes
#	next if $type ne "protein_coding" and $type ne "lincRNA" and $type ne "antisense";

	# Print Out
#	foreach my $gid (keys %{$data{by_type}{$type}}) {
#		foreach my $tid (keys %{$data{by_type}{$type}{$gid}}) {
#			my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $data{by_type}{$type}{$gid}{$tid}[0]);
#			print $outBed "$chr\t$start\t$end\t$name\t$zero\t$strand\n";
#		}
#	}
#}
#close $outBed;

#system("cat $output | sort -k1,1 -k2,2n > $output\_TEMP && mv $output\_TEMP $output");
#print "Use this as cache: $filename\n";

#################
## Subroutines ##
#################

sub parse_gtf {
	my ($input) = @_;
	my ($folder, $filename) = mitochy::getFilename($input, "folder");
	$filename =~ s/_All//i;
	my %data;

	# 1. Check if cache $filename exist
	my $data = $cache->get("$filename");

	# 2a. If not exist or forced to re-cache
	open (my $outGTF, ">", "$filename\_exon.gtf") or die "Cannot write to $filename\_exon.gtf: $!\n";
	open (my $outEXO, ">", "$filename\_exon.bed") or die "Cannot write to $filename\_exon.bed: $!\n";
	#if (not defined($data) or $opt_c) {
		open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
		my $lineCount = 0;
		while (my $line = <$in>) {
			chomp($line);
		
			next if $line =~ /#/;
			print "Done $lineCount\n" if $lineCount % 100000 == 0;
			$lineCount++;
		
			# Parse each line
			my ($chr, $source, $feat, $start, $end, $dot, $strand, $dot2, $name) = split("\t", $line);

			# Only take gene with KNOWN and NOVEL gene status
			next if $name !~ /transcript_status \"KNOWN\"\;/ and $name !~ /transcript_status \"NOVEL\"\;/;

			# For gtf, take all features but only protein_coding, lincRNA, or antisense
			my ($transcriptType) = $name =~ /transcript_type \"(\w+)\"\; transcript_status \"(\w+)\"; transcript/;
			print $outGTF "$line\n" if $transcriptType eq "protein_coding" or $transcriptType eq "lincRNA" or $transcriptType eq "antisense";

			# For bed take only transcripts
			next if $feat ne "exon";
			my ($gid, $tid, $type, $exon) = $name =~ /^gene_id "(.+)"; transcript_id "(.+)"; gene_type.+ transcript_type "(.+)"; transcript_status.+exon_number (\d+);/;
			die "\nDied at getting gene_id from name at line:\n$line\n\n" unless defined($gid) and defined($tid) and defined($type) and defined($exon);
	
			print $outEXO "$chr\t$start\t$end\t$tid\_$exon\t0\t$strand\n";
			# Group transcripts to gene and type
			#push(@{$data{by_gene}{$gid}{tid}{$tid}}, "$chr\t$start\t$end\t$tid\t0\t$strand");
			#$data{by_gene}{$gid}{type}        = $type;	
			#$data{by_tid}{$tid} = $gid;
			#push(@{$data{by_type}{$type}{$gid}{$tid}}, "$chr\t$start\t$end\t$tid\t0\t$strand");
		}
	
		# Group exons
		#foreach my $gid (keys %{$data{by_gene}}) {
		#	my %tid = %{$data{by_gene}{$gid}{tid}};
		#	$data{by_gene}{$gid}{coor} = group_gene(\%tid);
		#}
		#$cache->set($filename, \%data);
		close $in;
		close $outGTF;
		system("cat $filename\_exon.bed | sort -k1,1 -k2,2n -T /data/mitochi/Work/SortTMP/ > $filename\_exon\_TEMP && mv $filename\_exon\_TEMP $filename\_exon.bed");
		
	#}

	# 2b. Else just get data
	#else {
	#	%data = %{$data};
	#}
	
	#return(\%data);
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
