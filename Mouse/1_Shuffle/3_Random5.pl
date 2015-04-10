#!/usr/bin/perl
# This script get random 5 for metaplot
use strict; use warnings; use mitochy;

my ($input1, $input2) = @ARGV;
die "usage: $0 <dripc_promoter_shuf.metaplot> <dripcTSS_promoter_shuf.metaplot>\n" unless @ARGV;
die "usage: $0 <dripc_promoter_shuf.metaplot> <dripcTSS_promoter_shuf.metaplot>\n" unless -e $input1;
die "usage: $0 <dripc_promoter_shuf.metaplot> <dripcTSS_promoter_shuf.metaplot>\n" if $input1 =~ /orig/;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
my ($folder2, $fileName2) = mitochy::getFilename($input2, "folder");

my %data;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $value, $strand, $info) = split("\t", $line);
	my ($orig) = $info =~ /DRIP=(.+);ORIG=/;die if not defined($orig);
	push(@{$data{$orig}}, $line);
}
close $in1;

my %used;
open (my $out1, ">", "$fileName1\.metaplotout") or die "Cannot write to $fileName1\.metaplotout: $!\n";
foreach my $orig (keys %data) {
	my @random = shuffle(@{$data{$orig}});
	for (my $i = 0; $i < 5; $i++) {
		print $out1 "$random[$i]\n";
		my ($chr, $start, $end, $name, $value, $strand, $info) = split("\t", $random[$i]);
		$used{$info} = 1 if $input1 =~ /promoter\_shuf/i or $input1 =~ /terminal\_shuf/ or $input1 =~ /antisense_shuf/i;
	}
}
close $out1;

die "Done for non-prom/term/anti\n" if ($input1 !~ /promoter_shuf/i and $input1 !~ /terminal_shuf/i and $input1 !~ /antisense_shuf/i);
open (my $out2, ">", "$fileName2\.metaplotout") or die "Cannot write to $fileName2\.metaplotout: $!\n";
open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $value, $strand, $info) = split("\t", $line);
	my ($orig) = $info =~ /DRIP=(.+);ORIG=/;die if not defined($orig);
	if (defined($used{$info})) {
		if ($used{$info} == 1) {
			print $out2 "$line\n";
		}
	}
	
}
close $in2;
close $out2;

sub shuffle {
        my (@value) = @_;
        #print "Before: @value\n";
        for (my $i = 0; $i < 10000; $i++) {
                my $rand1 = int(rand(@value));
                my $rand2 = int(rand(@value));
                my $val1 = $value[$rand1];
                my $val2 = $value[$rand2];
                $value[$rand1] = $val2;
                $value[$rand2] = $val1;
        }
        #print "After: @value\n";
        return(@value);
}
