#!/usr/bin/perl
# Steps:
# 1. grep -P "\texon\t" main.gtf | cut -f1,4,5,7,9 | perl -pi -e 's/gene_id \"(.+)\"\; transcript_id.+$/$1/' | perl -pi -e 's/^(.+)\t(.+)\t(.+)\t(.+)\t(.+)$/$1\t$2\t$3\t$5\t0\t$4/' > exon.bed
# 2. bedtools merge -s -nms -i exon.bed > exon_merged.bed
# 3. ./1_MergeGene.pl exon_merged.bed
# 4. bed2gtf exon_merged.out

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";

my %group;
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $names, $strand) = split("\t", $line);
	my @names = split(";", $names);
	
	my $designated_group_number = "-1"; #Does not have any group yet

	# Check in %group if any of the name exist; if yes then put all genes in the same group
	foreach my $name (@names) {
		foreach my $group_number (keys %{$group{$chr}}) {
			if (defined($group{$chr}{$group_number}{names}{$name})) {
				#print "$name exist in group number $group_number\n";
				$designated_group_number = $group_number;
			}
			last if $designated_group_number != -1;
		}
		last if $designated_group_number != -1;
	}	

	if ($designated_group_number == -1) {
		$designated_group_number = (keys %{$group{$chr}})+1;
		#print "Group name @names is $designated_group_number\n";
		# Create new group
	}
	foreach my $name (@names) {
		$group{$chr}{$designated_group_number}{names}{$name} = 1;
		#print "$designated_group_number: ";
		foreach my $name (keys %{$group{$chr}{$designated_group_number}{names}}) {
		#	print "$name,";
		}
		#print "\n\n";
		#die if $designated_group_number > 10;
		$group{$chr}{$designated_group_number}{exon}{$start} = $end;
		$group{$chr}{$designated_group_number}{strand} = $strand;
	}
}

close $in;

print "Doing re-checking of groups to make sure all associated genes are grouped together\n";
my $check = 1;
my %chr;
my $last_group_number = 0;
my $last_chr;
while ($check != 0) {
	$check = 0;
	foreach my $chr (sort keys %group) {
		$last_chr = $chr if not defined($last_chr);
		print "Processing chr $chr\n" and $chr{$chr} = 2 if not defined($chr{$chr});
		print "\tChr $chr is done\n" and next if $chr{$chr} == 1;
		foreach my $group_number (sort {$a <=> $b} keys %{$group{$chr}}) {
			$last_group_number = $group_number if $group_number > $last_group_number;
			$last_group_number = $group_number if $chr ne $last_chr;
			$last_chr = $chr if $last_chr ne $chr;
			next if $group_number < $last_group_number - 200;
			next if $group_number > $last_group_number + 200;
			foreach my $name (keys %{$group{$chr}{$group_number}{names}}) {
				foreach my $group_number2 (sort {$a <=> $b} keys %{$group{$chr}}) {
					next if $group_number < $last_group_number - 200;
					last if $group_number > $last_group_number + 200;
					next if $group_number == $group_number2;
					if (defined($group{$chr}{$group_number2}{names}{$name})) {
						$check = 2;
						$chr{$chr} = 2;
						foreach my $start (keys %{$group{$chr}{$group_number2}{exon}}) {
							my $end = $group{$chr}{$group_number2}{exon}{$start};
							$group{$chr}{$group_number}{exon}{$start} = $end;
						}
						foreach my $name2 (keys %{$group{$chr}{$group_number2}{names}}) {
							$group{$chr}{$group_number}{names}{$name2} = 1;
						}
						print "Group chr $chr number $group_number name $name has another name defined in $group_number2 (deleted)\n";
						delete($group{$chr}{$group_number2});
					}
					last if $check == 2;
				}
				last if $check == 2;
			}
			last if $check == 2;
		}
		last if $check == 2;
		$chr{$chr} = 1;
	}
}

foreach my $chr (sort keys %group) {
	foreach my $group_number (sort {$a <=> $b} keys %{$group{$chr}}) {
		my $strand = $group{$chr}{$group_number}{strand};
		my $group_name;
		foreach my $name (sort keys %{$group{$chr}{$group_number}{names}}) {
			$group_name .= "$name,";
		}
		$group_name =~ s/,$//;
		foreach my $start (sort {$a <=> $b} keys %{$group{$chr}{$group_number}{exon}}) {
			my $end = $group{$chr}{$group_number}{exon}{$start};
			print $out "$chr\t$start\t$end\t$group_name\t0\t$strand\n";
		}
	}
}
close $out;


