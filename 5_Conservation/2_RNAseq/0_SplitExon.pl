#!/usr/bin/perl

use strict; use warnings; use mitochy qw(:all); use Cache::FileCache;

my $input = "/data/mitochi/Work/Project/DRIPc/gtf/gencode_vM2_annotation.gtf";

my ($folder, $filename) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
my $lineCount = 0;
my ($totalLine) = `wc -l $input` =~ /^(\d+)/;
my %data;
my %type;
my %gID;
while (my $line = <$in>) {
	chomp($line);
	$lineCount ++;
	printf "Processed $lineCount / $totalLine (%.2f %%)\n", $lineCount * 100 / $totalLine if $lineCount % 100000 == 0;
	next if $line =~ /#/;
	my ($chr, $source, $feature, $start, $end, $dot, $strand, $dot2, $name) = split("\t", $line);
	my ($gID, $tID, $transcriptType, $transcriptStatus) = $name =~ /gene_id \"(.+)\"; transcript_id \"(.+)\"; gene_type.+transcript_type \"(.+)\"; transcript_status \"(.+)\"; transcript_name/;
	die "Died undef ID or transcript at line $lineCount: $line\n" if not defined($gID) or not defined($tID) or not defined($transcriptType) or not defined($transcriptStatus);

	# Get Exon
	next if $feature ne "exon";
	# Get PC/Linc/Anti
	next if $transcriptType ne "protein_coding" and $transcriptType ne "lincRNA" and $transcriptType ne "antisense";
	# Get NOVEL and KNOWN only
	next if $transcriptStatus ne "NOVEL" and $transcriptStatus ne "KNOWN";
	
	push(@{$data{$gID}{$tID}}, "$chr\t$start\t$end\t$tID\t0\t$strand\n");
	$type{$tID} = $transcriptType;
	$gID{$tID} = $gID;
}

close $in;

mkdir "temp" if not -d "temp";

foreach my $gID (keys %data) {
	open (my $out, ">", "temp/$gID.temp") or die "Cannot write to temp/$gID.temp: $!\n";
	foreach my $tID (keys %{$data{$gID}}) {
		for (my $i = 0; $i < @{$data{$gID}{$tID}}; $i++) {
			print $out "$data{$gID}{$tID}[$i]\n";
		}
	}
	close $out;
}

# Merge each file
my @files = <temp/*.temp>;
open (my $out2, ">", "mm10_gencode_vM2_exon.gtf");
$lineCount = 0;
foreach my $file (@files) {
	$lineCount ++;
	printf "Processed $lineCount / @files (%.2f %%\n", $lineCount * 100 / @files if $lineCount % 5000 == 0;
	system("bedtools merge -s -nms -i $file > $file.merged");
	open (my $in2, "<", "$file.merged") or die "Cannot read from $file.merged: $!\n";
	while (my $line = <$in2>) {
		chomp($line);
		my ($chr, $start, $end, $names, $strand) = split("\t", $line);
		my @names = split(";", $names);
		my @newnames;
		my @gID;
		my @transcriptType;
		foreach my $name (@names) {
			push(@newnames, $name) if not grep(/^$names$/, @newnames);
			push(@gID, $gID{$name}) if not grep(/^$gID{$name}$/, @gID);
			push(@transcriptType, $type{$name}) if not grep(/^$type{$name}$/, @transcriptType);
			
		}
		my $newnames = join(";", @newnames);
		my $gIDnames = join(";", @gID);
		my $transcriptType = "protein_coding" if @transcriptType == 1 and $transcriptType[0] eq "protein_coding";
		$transcriptType = "mixed" if @transcriptType > 1 and grep(/protein_coding/, @transcriptType) and grep(/lincRNA/, @transcriptType);
		$transcriptType = "mixed" if @transcriptType > 1 and grep(/protein_coding/, @transcriptType) and grep(/antisense/, @transcriptType);
		$transcriptType = "noncoding" if @transcriptType > 1 and grep(/lincRNA/, @transcriptType) and grep(/antisense/, @transcriptType) and not grep(/protein_coding/, @transcriptType);
		$transcriptType = "lincRNA" if @transcriptType == 1 and grep(/lincRNA/, @transcriptType);
		$transcriptType = "antisense" if @transcriptType == 1 and grep(/antisense/, @transcriptType);

		print $out2 "$chr\tSOURCE\texon\t$start\t$end\t.\t$strand\t.\tgene_id \"$gIDnames\"; transcript_id \"$newnames\"; transcript_type \"$transcriptType\"\n";
	}
	close $in2;
}


__END__
my $cache = new Cache::FileCache;
$cache->set_cache_root("/data/mitochi/Work/Cache/");

my $data = $cache->get("gencode_v19_annotation.gtf");
my %data;
if (defined($data)) {
	$data = %{$data};
}
else {
	die "Undefined cache, please run gtf/Group.pl first\n";
}

