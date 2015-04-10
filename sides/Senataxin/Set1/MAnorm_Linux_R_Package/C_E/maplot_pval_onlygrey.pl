#!/usr/bin/env perl
use warnings;
use strict;

my $usage = "
usage: <MAnorm outfile .xls>

Warnings:
REQUIRES pre-processing, see comments in while loop
";

my ($ma_out) = @ARGV;
die $usage unless @ARGV;

open (IN, "<", $ma_out) or die "Could not open input $ma_out: $!\n";
open (OUT, ">", "manorm.out") or die "Could not open output manorm.out: $!\n";

while (<IN>) {
  my $line = $_;
  chomp $line;
	next if $line =~ /^chr\t/;

  #9 normal fields of MAnorm
  #cut and paste these new fields:
  # bedtools -c -a DRIPc -b NEAT1 followed by a -a DRIP #10 (check sample order) (requires DRIP/DRIPc to be split out and resorted)
  # "" with MALAT1 #11
  my @fields = split(/\s+/, $line);
  my $newline = $line;
  if ($fields[8] > 2) {
    #p-value cut offs

    if ($fields[8] > 2 && $fields[8] <= 4) {
      $newline.="\t1";
    }
    if ($fields[8] > 4 && $fields[8] <= 6) {
      $newline.="\t2";
    }
    if ($fields[8] > 6 && $fields[8] <= 8) {
      $newline.="\t3";
    }
    if ($fields[8] > 8 && $fields[8] <= 10) {
      $newline.="\t4";
    }
    if ($fields[8] > 10) {
      $newline.="\t5";
    }
  }
  else {
    $newline.="\t0";
  }
  print OUT "$newline\n";
} 

close IN;
close OUT;

system("run_Rscript.pl plot_manorm.R");
