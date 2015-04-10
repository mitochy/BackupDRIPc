#!/usr/bin/perl
# 0. Parse data/NT2.rpkm to get RNA-seq of each genes
# 1. Input is bed/hg19_gencode.bed, run bedtools_bed_change.pl
# 2. Divide into two: Those with DRIPc and those without using bedtools intersect
# 3. Divide genes into expression quintiles: 17, 56, 142 log(0-3, 3-4, 4-5, 5-max)

use strict; use warnings; use mitochy; use Getopt::Std;
use vars qw($opt_r $opt_g $opt_d $opt_s);
getopts("r:g:d:s");

die "usage: $0 -r <rnaseq> -g <hg19_gencode.bed> -d <dripc.bed>\nOptional:\n-s:Strandedness is off\n" unless defined($opt_g) and defined($opt_d);

my $geneFile = $opt_g;
my $dripFile = $opt_d;
my $rnaFile  = $opt_r;
my ($folder, $filename) = mitochy::getFilename($geneFile, "folder");

#my $geneFile = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";

# 0. Parse NT2.rpkm
my %rna = %{parse_rna($rnaFile)};

# 1. Divide into promoter and terminal
runbash("bedtools_bed_change.pl -a -x -2000 -y 2000 -i $geneFile -o $filename\_promoter.temp");
runbash("bedtools_bed_change.pl -b -x -2000 -y 2000 -i $geneFile -o $filename\_terminal.temp");

# 2. Divide promoter and terminal into those with DRIPc and not
runbash("bedtools intersect -u -s -a $filename\_promoter.temp -b $dripFile > $filename\_promoter_drip.temp");
runbash("bedtools intersect -v -s -a $filename\_promoter.temp -b $dripFile > $filename\_promoter_nodrip.temp");
runbash("bedtools intersect -u -s -a $filename\_terminal.temp -b $dripFile > $filename\_terminal_drip.temp");
runbash("bedtools intersect -v -s -a $filename\_terminal.temp -b $dripFile > $filename\_terminal_nodrip.temp");

# 3. Divide into log(RNA-seq+1): 1 is zero, 1-3 is low, 3-4 is med, 4-5 is high, 5-max is super high
divide_by_exp("$filename\_promoter_drip.temp", "$filename\_promoter_nodrip.temp", "$filename\_terminal_drip.temp", "$filename\_terminal_nodrip.temp");

sub divide_by_exp {
	my (@input) = @_;
	foreach my $input (@input) {
		my ($folder, $fileName) = mitochy::getFilename($input, "folder");
		open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
		open (my $outzero,  ">", "$fileName\_zero.temp")  or die "Cannot write to $fileName\_zero.temp: $!\n";
		open (my $outlow,   ">", "$fileName\_low.temp")   or die "Cannot write to $fileName\_low.temp: $!\n";
		open (my $outmed,   ">", "$fileName\_med.temp")   or die "Cannot write to $fileName\_med.temp: $!\n";
		open (my $outhigh,  ">", "$fileName\_high.temp")  or die "Cannot write to $fileName\_high.temp: $!\n";
		open (my $outsuper, ">", "$fileName\_super.temp") or die "Cannot write to $fileName\_super.temp: $!\n";

		# If the exact same position, then use it once. %data to keep track
		# If same position, use the one with highest RNA-seq
		my %data;
		while (my $line = <$in>) {
			chomp($line);
			next if $line =~ /#/;
			my ($chr, $start, $end, $name, $val, $strand) = split("\t", $line);
			my $rna = defined($rna{$name}) ? $rna{$name} : 0;
			
			if (defined($data{$chr}{$start}{$end}{rna})) {
				my $prev = $data{$chr}{$start}{$end}{line};
				my ($chr0, $start0, $end0, $name0) = split("\t", $prev);
				my $newname = "$name0;$name";
				$data{$chr}{$start}{$end}{line} = "$chr\t$start\t$end\t$newname\t$val\t$strand";
				$data{$chr}{$start}{$end}{rna} = log($rna+1) if $data{$chr}{$start}{$end}{rna} < log($rna+1);

			}
			else {
				$data{$chr}{$start}{$end}{line} = "$chr\t$start\t$end\t$name\t$val\t$strand";
				$data{$chr}{$start}{$end}{rna} = log($rna+1);
			}
		}

		# Print out
		foreach my $chr (sort keys %data) {
			foreach my $start (sort {$a <=> $b} keys %{$data{$chr}}) {
				foreach my $end (sort keys %{$data{$chr}{$start}}) {
					my $data = $data{$chr}{$start}{$end}{line};
					my $rna = $data{$chr}{$start}{$end}{rna};
					print $outzero  "$data\n" if $rna <= 1;
					print $outlow   "$data\n" if $rna > 1 and $rna <= 3;
					print $outmed   "$data\n" if $rna > 3 and $rna <= 4;
					print $outhigh  "$data\n" if $rna > 4 and $rna <= 5;
					print $outsuper "$data\n" if $rna > 5;
				}
			}
		}
		close $in;
	}
}

sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
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


__END__
sub {
	my ($input) = @_;
	my %data;
	open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		next if $line =~ /#/;
		my @arr = split("\t", $line);
	}
	close $in;
}
