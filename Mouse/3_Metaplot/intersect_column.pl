#!/usr/bin/perl
use strict; use warnings; use Getopt::Std;
use vars qw($opt_o $opt_e $opt_n $opt_b $opt_t);
getopts("o:enbt:");

my ($query, $query_col, $subject, $subject_col) = @ARGV;
die "usage: $0 [options] <queryfile> <que data col> <subject file> <subj data col>
options:
-o: output
-e: exists only
-n: print not exists only
-b: print both
-t: delimiter (space or tab) [default: tab]
" unless @ARGV == 4;
my $delimiter = defined($opt_t) ? $opt_t : "tab";
my $output = "$query.intersected_exists.bed" if defined($opt_e) or defined($opt_b);
my $output2 = "$query.intersected_notexists.bed" if defined($opt_n) or defined($opt_b);

die "delimiter must be space or tab\n" unless $delimiter eq "space" or $delimiter eq "tab";
open (my $que, "<", $query) or die "Cannot read from $query: $!\n";
open (my $sub, "<", $subject) or die "Cannot read from $subject: $!\n";
open (my $out, ">", $output) or die "Cannot write to $output: $!\n" if (defined($output));
open (my $out3, ">", "$output.exist2") or die "Cannot write to $output: $!\n" if (defined($output));
open (my $out2, ">", $output2) or die "Cannot write to $output2: $!\n" if (defined($output2));
my $print = 0;
my %query;
while (my $line = <$que>) {
	chomp($line);
	next if $line =~ /^"V1/i;
	next if $line =~ /Ensembl/;
	$line =~ s/"//ig;
	my @arr = split("\t", $line) if $delimiter eq "tab";
	@arr = split(" ", $line) if $delimiter eq "space";
	print "$line undefined column at query\n" and next if not defined($arr[$query_col]);
	$query{data}{$arr[$query_col]}++;
	$query{line}{$arr[$query_col]} = $line;
}
close $que;
my %sub;
while (my $line = <$sub>) {
	chomp($line);
	next if $line =~ /^"V1/i;
	next if $line =~ /Ensembl/;
	$line =~ s/"//ig;
	my @arr = split("\t", $line) if $delimiter eq "tab";
	@arr = split(" ", $line) if $delimiter eq "space";
	print "$line undefined column at subject\n" and next if not defined($arr[$subject_col]);
	$sub{data}{$arr[$subject_col]}++;
 	print "Example subject: $arr[$subject_col]\n" if $print == 0;
        $print = 1;
}

$print = 0;
my $count = 0;
my %data;
foreach my $que (sort keys %{$query{data}}) {
	print "Example query: $que\n" if $print == 0;
	$print = 1;
	$count++ if exists($sub{data}{$que});
	#print "$que = $sub{$que}\n";
	$data{notexist}{data}{$que}++ if not exists($sub{data}{$que});
	$data{notexist}{line}{$que} = $query{line}{$que} if not exists($sub{data}{$que});
	$data{exist}{data}{$que}++ if exists($sub{data}{$que});

	$data{exist}{line}{$que} = $query{line}{$que} if exists($sub{data}{$que});
}

my $querycount = (keys %{$query{data}});
my $subjectcount = (keys %{$sub{data}});
print "query: $query subj: $subject\n";
print "consensus = $count\n";
print "query = $querycount\n";
print "subject = $subjectcount\n";

$count = 0;
print "Exists: \n";
foreach my $que (keys %{$data{exist}{data}}) {
	$count++;
	print "$que\n" if $count <= 10;
}
if (defined($output)) {
	open (my $que2, "<", $query) or die "Cannot read from $query: $!\n";
	while (my $line = <$que2>) {
		chomp($line);
		next if $line =~ /^"V1/i;
		next if $line =~ /Ensembl/;
		$line =~ s/"//ig;
		my @arr = split("\t", $line) if $delimiter eq "tab";
		@arr = split(" ", $line) if $delimiter eq "space";
		next if not defined($arr[$query_col]);
		print "LINE $line\n" if $line =~ /ENSG00000160766/;
		if (defined($data{exist}{line}{$arr[$query_col]}) and $data{exist}{line}{$arr[$query_col]} !~ /^$/) {
			print $out "$data{exist}{line}{$arr[$query_col]}\n" 
		}
	}
	close $que2;
	open (my $sub2, "<", $subject) or die "Cannot read from $subject: $!\n";
	while (my $line = <$sub2>) {
		chomp($line);
		next if $line =~ /^"V1/i;
		next if $line =~ /Ensembl/;
		$line =~ s/"//ig;
		my @arr = split("\t", $line) if $delimiter eq "tab";
		@arr = split(" ", $line) if $delimiter eq "space";
		next if not defined($arr[$query_col]);
		print "LINE $line\n" if $line =~ /ENSG00000160766/;
		if (defined($data{exist}{line}{$arr[$query_col]}) and $data{exist}{line}{$arr[$query_col]} !~ /^$/) {
			print $out3 "$line\n";
		}
	}
	close $sub2;
}
print "Not exists: \n";
$count = 0;
foreach my $que (keys %{$data{notexist}{data}}) {
	$count++;
	print "$que\n" if $count == 10;
	print $out2 "$data{notexist}{line}{$que}\n" if defined($output2);
}
if (defined($output2)) {
	open (my $que2, "<", $query) or die "Cannot read from $query: $!\n";
	while (my $line = <$que2>) {
		chomp($line);
		next if $line =~ /^"V1/i;
		next if $line =~ /Ensembl/;
		$line =~ s/"//ig;
		my @arr = split("\t", $line) if $delimiter eq "tab";
		@arr = split(" ", $line) if $delimiter eq "space";
		next if not defined($arr[$query_col]);
		if (defined($data{notexist}{line}{$arr[$query_col]}) and $data{notexist}{line}{$arr[$query_col]} !~ /^$/) {
			print $out2 "$data{notexist}{line}{$arr[$query_col]}\n" 
		}
	}
	close $que2;
}
