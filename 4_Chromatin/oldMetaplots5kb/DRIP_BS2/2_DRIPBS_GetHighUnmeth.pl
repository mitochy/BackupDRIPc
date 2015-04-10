#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1, $input2) = @ARGV;
die "usage: $0 <input_orig.txt from map_wig -p> <input_shuf.txt>\n" unless @ARGV;

my %rna = %{parse_rna("/data/mitochi/Work/Project/DRIPc/data/NT2.rpkm")};
my ($folder, $fileName) = mitochy::getFilename($input1, "folder");


my %data;
my %change;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
	my ($chr, $start, $end, $vals, $name) = split("\t", $line);
	next if $vals eq "NA" or $vals eq 0;
	die "Died at $line\n" if not defined($vals);
	my (@vals) = split("_", $vals);
	my $count = 0;
	foreach my $val (@vals) {
		my ($pos, $val) = split(",", $val);
		$data{$name}{orig}{count}{$pos}{10}++ if $val < 10;
		$data{$name}{orig}{count}{$pos}{20}++ if $val >= 10 and $val < 20;
		$data{$name}{orig}{count}{$pos}{30}++ if $val >= 20 and $val < 30;
		$data{$name}{orig}{count}{$pos}{40}++ if $val >= 30 and $val < 40;
		$data{$name}{orig}{count}{$pos}{50}++ if $val >= 40 and $val < 50;
		$data{$name}{orig}{count}{$pos}{90}++ if $val >= 90;
	}
	
	# Get RNAseq
	my $rna = $rna{$name};
	die "Died at $name not defined RNA\n" if not defined($rna);
	my $label;
	$label = "0_zero"  if $rna < 10;
	$label = "1_low"   if $rna >= 10 and $rna < 50;
	$label = "2_med"   if $rna >= 50 and $rna < 100;
	$label = "3_high"  if $rna >= 100 and $rna < 200;
	$label = "4_super" if $rna > 200;
	$change{total} ++ if $rna < 10;

	$data{$name}{orig}{coor} = "$chr\t$start\t$end";
	$data{$name}{rna} = $label;
}
close $in1;
open (my $in2, "<", $input2) or die "Cannot read from $input2: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
	my ($chr, $start, $end, $vals, $name) = split("\t", $line);
	next if $vals eq "0" or $vals eq "NA";
	my (@vals) = split("_", $vals);
	my $count = 0;
	foreach my $val (@vals) {
		my ($pos, $val) = split(",", $val);
		next if not defined($data{$name}{orig}{count}{$pos});
		$data{$name}{shuf}{count}{$pos}{10} ++ if $val < 10;
		$data{$name}{shuf}{count}{$pos}{50} ++ if $val <= 60;
		$data{$name}{shuf}{count}{$pos}{90}++ if $val >= 90;
		$data{$name}{total} ++;
	}
	$data{$name}{shuf}{coor} = "$chr\t$start\t$end";
}
close $in2;

open (my $out1, ">", "$fileName\_hyper10.out") or die "Cannot write to $fileName.out: $!\n";
open (my $out7, ">", "$fileName\_bothhigh10.out") or die "Cannot write to $fileName.out: $!\n";
open (my $out8, ">", "$fileName\_bothlow10.out") or die "Cannot write to $fileName.out: $!\n";
open (my $out2, ">", "$fileName\_hyper20.out") or die "Cannot write to $fileName.out: $!\n";
open (my $out3, ">", "$fileName\_hyper30.out") or die "Cannot write to $fileName.out: $!\n";
open (my $out4, ">", "$fileName\_hyper40.out") or die "Cannot write to $fileName.out: $!\n";
open (my $out5, ">", "$fileName\_hyper50.out") or die "Cannot write to $fileName.out: $!\n";
open (my $out6, ">", "$fileName\_hypo.out") or die "Cannot write to $fileName.out: $!\n";
#open (my $out3, ">", "$fileName\_none.out") or die "Cannot write to $fileName.out: $!\n";

foreach my $name (sort {$data{$a}{rna} cmp $data{$b}{rna}} keys %data) {
	my %count;
	my $count_total = 0;
	next if not defined($data{$name}{orig}{count}) or keys %{$data{$name}{orig}{count}} == 0;
	foreach my $pos (keys %{$data{$name}{orig}{count}}) {
		my $count1_10 = defined($data{$name}{orig}{count}{$pos}{10}) ? $data{$name}{orig}{count}{$pos}{10} : 0;
		my $count1_20 = defined($data{$name}{orig}{count}{$pos}{20}) ? $data{$name}{orig}{count}{$pos}{20} : 0;
		my $count1_30 = defined($data{$name}{orig}{count}{$pos}{30}) ? $data{$name}{orig}{count}{$pos}{30} : 0;
		my $count1_40 = defined($data{$name}{orig}{count}{$pos}{40}) ? $data{$name}{orig}{count}{$pos}{40} : 0;
		my $count1_50 = defined($data{$name}{orig}{count}{$pos}{50}) ? $data{$name}{orig}{count}{$pos}{50} : 0;
		my $count1_90 = defined($data{$name}{orig}{count}{$pos}{90}) ? $data{$name}{orig}{count}{$pos}{90} : 0;
		my $count2_10 = defined($data{$name}{shuf}{count}{$pos}{10}) ? $data{$name}{shuf}{count}{$pos}{10} : 0;
		my $count2_50 = defined($data{$name}{shuf}{count}{$pos}{50}) ? $data{$name}{shuf}{count}{$pos}{50} : 0;
		my $count2_90 = defined($data{$name}{shuf}{count}{$pos}{90}) ? $data{$name}{shuf}{count}{$pos}{90} : 0;
		$count{orig}{10} += $count1_10;
		$count{orig}{20} += $count1_20;
		$count{orig}{30} += $count1_30;
		$count{orig}{40} += $count1_40;
		$count{orig}{50} += $count1_50;
		$count{orig}{90} += $count1_90;
		$count{shuf}{10} += $count2_10;
		$count{shuf}{50} += $count2_50;
		$count{shuf}{90} += $count2_90;
		$count_total ++;
	}
	next if $count_total / $data{$name}{total} < 0.3; # Next if DRIPBS has 70% lower C assessed than Bulk
	next if $data{$name}{total} < 3;
	my $rna = $data{$name}{rna};
	my $rnaval = $rna{$name};
	my $ratio_10 = $count{shuf}{50} == 0 ? $count{orig}{10} : $count{orig}{10} / $count{shuf}{50};
	my $ratio_20 = $count{shuf}{50} == 0 ? $count{orig}{20} : $count{orig}{20} / $count{shuf}{50};
	my $ratio_30 = $count{shuf}{50} == 0 ? $count{orig}{30} : $count{orig}{30} / $count{shuf}{50};
	my $ratio_40 = $count{shuf}{50} == 0 ? $count{orig}{40} : $count{orig}{40} / $count{shuf}{50};
	my $ratio_50 = $count{shuf}{50} == 0 ? $count{orig}{50} : $count{orig}{50} / $count{shuf}{50};

	if ($count{orig}{10} >= 5 and $count{shuf}{50} <= 1) {
		#print "\tout1: $name $count{orig}{10} $count{shuf}{50}\n";
		print $out1 "$data{$name}{orig}{coor}\t$name\t$rna\t$ratio_10\t$data{$name}{total}\n";
		$change{hyper} ++ if $rnaval < 10;
	}
	elsif ($count{orig}{10} >= 5 and $count{shuf}{50} > 1 and $ratio_10 > 2) {
		#print "\tout1: $name $count{orig}{10} $count{shuf}{50}\n";
		print $out1 "$data{$name}{orig}{coor}\t$name\t$rna\t$ratio_10\t$data{$name}{total}\n";
		$change{hyper} ++ if $rnaval < 10;
	}
	elsif ($count{orig}{10} >= 5 and $count{orig}{10} / $data{$name}{total} > 0.5 and $count{shuf}{50} > 1 and $ratio_10 > 1.2) {
		#print "\tout1: $name $count{orig}{10} $count{shuf}{50}\n";
		print $out1 "$data{$name}{orig}{coor}\t$name\t$rna\t$ratio_10\t$data{$name}{total}\n";
		$change{hyper} ++ if $rnaval < 10;
	}
	elsif ($count{orig}{10} >= 15 and $count{shuf}{50} > 1 and $ratio_10 > 1.2) {
		#print "\tout1: $name $count{orig}{10} $count{shuf}{50}\n";
		print $out1 "$data{$name}{orig}{coor}\t$name\t$rna\t$ratio_10\t$data{$name}{total}\n";
		$change{hyper} ++ if $rnaval < 10;
	}
	#elsif ($count{orig}{10} >= 15 and $count{shuf}{50} > 1 and $count{orig}{10} / $count{shuf}{50} > 1.2) {
	#	#print "\tout1: $name $count{orig}{10} $count{shuf}{50}\n";
	#	print $out1 "$data{$name}{orig}{coor}\t$name\t$rna\t$ratio_10\t$data{$name}{total}\n";
	#	$change{hyper} ++ if $rnaval < 10;
	#}


	elsif ($count{orig}{10} >= 5 and $count{shuf}{10} >= 5) {
		$change{bothhigh} ++ if $rnaval < 10;
		print $out8 "$data{$name}{orig}{coor}\t$rna\t$ratio_10\t$data{$name}{total}\n";
	}
	elsif ($count{orig}{90} >= 5 and $count{shuf}{90} >= 5 and $count{orig}{90} / $data{$name}{total} > 0.9) {
		$change{bothlow} ++ if $rnaval < 10;
		print $out7 "$data{$name}{orig}{coor}\t$rna\t$ratio_10\t$data{$name}{total}\n";
	}
	elsif ($count{shuf}{10} >= 5 and $count{orig}{50} <= 1) {
		$change{hypo} ++ if $rnaval < 10;
	}
	elsif ($count{shuf}{10} >= 5 and $count{orig}{50} > 1 and $count{shuf}{10} / $count{shuf}{50} > 2) {
		$change{hypo} ++ if $rnaval < 10;
	}
	elsif ($count{orig}{90} >= 5 and $count{shuf}{50} >= 5) {
		$change{hypo} ++ if $rnaval < 10;
	}

	elsif ($count{orig}{20} >= 5 and $count{shuf}{50} <= 1) {
		#print "\tout2: $name $count{orig}{20} $count{shuf}{50}\n";
		print $out2 "$data{$name}{orig}{coor} $rna $ratio_20\n";
	}
	elsif ($count{orig}{20} >= 5 and $count{shuf}{50} > 1 and $count{orig}{20} / $count{shuf}{50} > 2) {
		#print "\tout2: $name $count{orig}{20} $count{shuf}{50}\n";
		print $out2 "$data{$name}{orig}{coor} $rna $ratio_20\n";
	}
	elsif ($count{orig}{30} >= 5 and $count{shuf}{50} <= 1) {
		#print "\tout3: $name $count{orig}{30} $count{shuf}{50}\n";
		print $out3 "$data{$name}{orig}{coor} $rna $ratio_30\n";
	}
	elsif ($count{orig}{30} >= 5 and $count{shuf}{50} > 1 and $count{orig}{30} / $count{shuf}{50} > 2) {
		#print "\tout3: $name $count{orig}{30} $count{shuf}{50}\n";
		print $out3 "$data{$name}{orig}{coor} $rna $ratio_30\n";
	}
	elsif ($count{orig}{40} >= 5 and $count{shuf}{50} <= 1) {
		#print "\tout4: $name $count{orig}{40} $count{shuf}{50}\n";
		print $out4 "$data{$name}{orig}{coor} $rna $ratio_40\n";
	}
	elsif ($count{orig}{40} >= 5 and $count{shuf}{50} > 1 and $count{orig}{40} / $count{shuf}{50} > 2) {
		#print "\tout4: $name $count{orig}{40} $count{shuf}{50}\n";
		print $out4 "$data{$name}{orig}{coor} $rna $ratio_40\n";
	}
	elsif ($count{orig}{50} >= 5 and $count{shuf}{50} <= 1) {
		#print "\tout5: $name $count{orig}{50} $count{shuf}{50}\n";
		print $out5 "$data{$name}{orig}{coor} $rna $ratio_50\n";
	}
	elsif ($count{orig}{50} >= 5 and $count{shuf}{50} > 1 and $count{orig}{50} / $count{shuf}{50} > 2) {
		#print "\tout5: $name $count{orig}{50} $count{shuf}{50}\n";
		print $out5 "$data{$name}{orig}{coor} $rna $ratio_50\n";
	}
}
print "Output:
1. Hypermethylated in $input1: $fileName\_hyper.out
2. Hypomethylated in $input1: $fileName\_hypo.out
";

foreach my $type (keys %change) {
	print "$type\t$change{$type}\n";
}
sub parse_rna {
        my ($input) = @_;
        my %data;
        open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
        while (my $line = <$in>) {
                chomp($line);
                next if $line =~ /#/;
                my ($gene, $val) = split("\t", $line);
                $data{$gene} = $val;
        }
        close $in;
        return(\%data);
}

