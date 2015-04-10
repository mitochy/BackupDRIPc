#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <gtf from Ensembl>\n" unless @ARGV == 1;
my ($outName) = mitochy::getFilename($input);
my $output = "$outName.bed";
my %bed;
my %gene;
open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", $output) or die "Cannot write to $output: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /^\#/;
	my ($chr, $junk0, $junk1, $start, $end, $dot, $strand, $dot2, $name) = split("\t", $line);
	my @name = split(";", $name);
	my $gene_id;
	my $tran_id;
	my $type;
	foreach my $names (@name) {
		$names =~ s/\s{1,10}//g;
		($gene_id) = $names =~ /^gene_id"(.+)"$/ if $names =~ /^gene_id/;
		($tran_id) = $names =~ /^transcript_id"(.+)"$/ if $names =~ /^transcript_id/;
		($type)    = $names =~ /^gene_biotype"(.+)"$/ if $names =~ /^gene_biotype/;
	}
	$type = $junk0 if not defined($type);# and print "Undef type: using $junk0 at $line\n" if not defined($type);
	die "Undef gene_id at $line\n" if not defined($gene_id);
	die "Undef gene_id at $line\n" if $gene_id =~ /^$/;
	#if ($type eq "protein_coding") {
		$gene{$gene_id} = 1;
		if (defined($bed{$gene_id}{val})) {
			print "Previous type: $bed{$gene_id}{val}. Current: $type\n" if $bed{$gene_id}{val} ne $type;
			die "$line\n" if $bed{$gene_id}{val} ne $type;
		}
	#}
	my $chr_type;
	if ($chr =~ /^\d+$/) {
		$chr_type = "numeric";
	}
	else {
		$chr_type = "alphabet";
	}
	die if not defined($chr_type);
	$bed{$chr_type}{$chr}{$gene_id}{start}  = $start if not defined($bed{$chr_type}{$chr}{$gene_id}{start}) or $start < $bed{$chr_type}{$chr}{$gene_id}{start};
	$bed{$chr_type}{$chr}{$gene_id}{end}    = $end   if not defined($bed{$chr_type}{$chr}{$gene_id}{end})   or $end > $bed{$chr_type}{$chr}{$gene_id}{end}    ;
	$bed{$chr_type}{$chr}{$gene_id}{val}    = $type;
	$bed{$chr_type}{$chr}{$gene_id}{strand} = $strand;
	$bed{$chr_type}{$chr}{$gene_id}{tID}    = $tran_id;
}
close $in;
my $count = 0;
my $total = 0;
foreach my $chr (sort {$bed{numeric}{$a} <=> $bed{numeric}{$b}} keys %{$bed{numeric}}) {
	foreach my $gene_id (sort {$bed{numeric}{$chr}{$a}{start} <=> $bed{numeric}{$chr}{$b}{start}} keys %{$bed{numeric}{$chr}}) {
		my $start   = $bed{numeric}{$chr}{$gene_id}{start};
		my $end     = $bed{numeric}{$chr}{$gene_id}{end};
		my $val     = $bed{numeric}{$chr}{$gene_id}{val};
		my $strand  = $bed{numeric}{$chr}{$gene_id}{strand};
		my $tID     = $bed{numeric}{$chr}{$gene_id}{tID};
		print $out "$chr\t$start\t$end\t$gene_id\_$tID\t$val\t$strand\n";
		$total++;
		$count++ if $val eq "protein_coding";
	}
}
foreach my $chr (sort {$bed{alphabet}{$a} cmp $bed{alphabet}{$b}} keys %{$bed{alphabet}}) {
	foreach my $gene_id (sort {$bed{alphabet}{$chr}{$a}{start} <=> $bed{alphabet}{$chr}{$b}{start}} keys %{$bed{alphabet}{$chr}}) {
		my $start   = $bed{alphabet}{$chr}{$gene_id}{start};
		my $end     = $bed{alphabet}{$chr}{$gene_id}{end};
		my $val     = $bed{alphabet}{$chr}{$gene_id}{val};
		my $strand  = $bed{alphabet}{$chr}{$gene_id}{strand};
		my $tID     = $bed{alphabet}{$chr}{$gene_id}{tID};
		print $out "$chr\t$start\t$end\t$gene_id\_$tID\t$val\t$strand\n";
		$total++;
		$count++ if $val eq "protein_coding";
	}
}
print "$input Total = $total, Protein_coding = $count\n";
