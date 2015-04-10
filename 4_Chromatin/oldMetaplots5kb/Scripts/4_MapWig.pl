#!/usr/bin/perl

use strict; use warnings; use mitochy;

my (@input) = <./*.input>;
foreach my $input (@input) {
	#next if $input =~ /genebody/;
	#next if $input =~ /intergenic/;
	next if $input !~ /antisense_other/ and $input !~ /intergenic/;
	#next if $input =~ /both/;
	runbash("run_script_in_paralel2.pl -v \"map_wig_to_bed_BIG.pl -w FILENAME -p $input\" /data/mitochi/Work/Project/DRIPc/4_Chromatin/Chromatin/signal/promoter_Set/ wig 15");
}
#foreach my $input (@input) {
#	next if $input !~ /genebody/;
#	runbash("run_script_in_paralel2.pl -v \"map_wig_to_bed_BIG.pl -w FILENAME -p $input\" /data/mitochi/Work/Project/DRIPc/4_Chromatin/Chromatin/signal/promoter_Set/ wig 15");
#}

sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}
