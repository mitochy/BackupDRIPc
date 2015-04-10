#!/usr/bin/perl

use strict; use warnings; use mitochy; use Cache::FileCache;

my ($inputs) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my @wig = qw(
/data/mitochi/Work/Project/DRIPc/5_Conservation/E14_hg19.wig
/data/mitochi/Work/Project/DRIPc/5_Conservation/3T3_hg19.wig
/data/mitochi/Work/Project/DRIPc/5_Conservation/NT2_hg192.wig
/data/mitochi/Work/Project/DRIPc/5_Conservation/Fibro_hg192.wig
);

my $cache = new Cache::FileCache;
$cache->set_cache_root("/data/mitochi/Work/Cache/");

my @filename = qw(E14_hg19 3T3_hg19 NT2_hg192 Fibro_hg192); 
my @chr = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY);

# TO ERASE CACHE UNHASH THESE
#foreach my $chr (@chr) {
#	my %data;
#	$cache->set("ConsWig.$chr", \%data);
#}

my $data = $cache->get("ConsWig.chr22");
my %data;
if (not defined($data)) {
	# Parse each wig file into %data
	foreach my $input (@wig) {
		my ($folder, $fileName) = mitochy::getFilename($input, "folder");
		open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
		my $chr = "INIT";
		my $linecount = 0;
		while (my $line = <$in>) {
			chomp($line);
			$linecount++;
			next if $line =~ /track/;
			if ($line =~ /variable/) {
				if ($chr ne "INIT") {
					$cache->set("ConsWig.$chr", \%data);
					%data = ();
				}
				($chr) = $line =~ /chrom=(.+) span/;
				my $data = $cache->get("ConsWig.$chr");
				if (defined($data)) {
					%data = %{$data};
				}
				print "Processing $input chromosome $chr\n";
				die "Died at line $line not defined chromosome\n" unless defined($chr);
			}
			else {
				my ($pos, $val) = split("\t", $line);
				my $pos_index = int($pos / 500);
				$data{$chr}{$pos_index}{$fileName}{value} += $val;
				$data{$chr}{$pos_index}{$fileName}{count} ++;
			}
		}
		$cache->set("ConsWig.$chr", \%data);
		%data = ();
		
		close $in;
	}
}

#=comments

open (my $out1, ">", "E14_3T3.compare");
open (my $out2, ">", "E14_NT2.compare");
open (my $out3, ">", "E14_Fibro.compare");
open (my $out4, ">", "3T3_NT2.compare");
open (my $out5, ">", "3T3_Fibro.compare");
open (my $out6, ">", "NT2_Fibro.compare");
foreach my $chr (sort @chr) {
	print "Processing chr $chr\n";
	$data = $cache->get("ConsWig.$chr");
	%data = %{$data};
	foreach my $pos (keys %{$data{$chr}}) {
		# Discard if all data is les than 3
		my $check = 0;
		my @used;
		for (my $i = 0; $i < @filename; $i++) {
			my $filename1 = $filename[$i];
			for (my $j = 0; $j < @filename; $j++) {
				my $filename2 = $filename[$j];
				next if grep(/^$filename1\_$filename2$/, @used);
				next if grep(/^$filename2\_$filename1$/, @used);
				#print "$filename1 $chr $pos VALUE $data{$chr}{$pos}{$filename1}{value} COUNT $data{$chr}{$pos}{$filename1}{count}\n" if $pos == 1104601;
				my $val1 = defined($data{$chr}{$pos}{$filename1}{value}) ? $data{$chr}{$pos}{$filename1}{value} / $data{$chr}{$pos}{$filename1}{count} : 0;
				my $val2 = defined($data{$chr}{$pos}{$filename2}{value}) ? $data{$chr}{$pos}{$filename2}{value} / $data{$chr}{$pos}{$filename2}{count} : 0;
				next if $val1 <= 3 or $val2 <= 3;
				#next if $val1 <= 3 and $val2 <= 3;
				my $start = $pos * 500;
				my $end = $pos * 500 + 499;
				print $out1 "$chr\t$start\t$end\t$val1\t$val2\t+\n" if $filename1 =~ /E14/ and $filename2 =~ /3T3/;
				print $out2 "$chr\t$start\t$end\t$val1\t$val2\t+\n" if $filename1 =~ /E14/ and $filename2 =~ /NT2/;
				print $out3 "$chr\t$start\t$end\t$val1\t$val2\t+\n" if $filename1 =~ /E14/ and $filename2 =~ /Fibro/;
				print $out4 "$chr\t$start\t$end\t$val1\t$val2\t+\n" if $filename1 =~ /3T3/ and $filename2 =~ /NT2/;
				print $out5 "$chr\t$start\t$end\t$val1\t$val2\t+\n" if $filename1 =~ /3T3/ and $filename2 =~ /Fibro/;
				print $out6 "$chr\t$start\t$end\t$val1\t$val2\t+\n" if $filename1 =~ /NT2/ and $filename2 =~ /Fibro/;
				push(@used, "$filename1\_$filename2");
			}
		}
	}
}

__END__
open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";
