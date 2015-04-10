#!/usr/bin/perl
# To make sure gene list used after Pipeline.pl is representative of the whole gene list
# Check Chr
# Check Length

use strict; use warnings; use Cache::FileCache; use Statistics::Basic qw(:all); use R_toolbox;

my ($input) = @ARGV;
die "Usage: $0 <TEMP_promoter_nogene.bed>\n" unless @ARGV;

my $cache = new Cache::FileCache();
$cache->set_cache_root("/data/mitochi/Work/Cache/");
check_gid($input);

sub check_gid {
        my ($input) = @_;
        # Check how many genes are used out of protein coding genes
        my %data = %{$cache->get("gencode_v19_annotation")};
        my %gid;
	my %stats;
        my $total_type = 0;
	my ($total_chr1, $total_chr2) = (0,0);
        open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
        while (my $line = <$in>) {
                next if $line =~ /^\#/;
                my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);

                # Get Gene ID of transcript (or unknown)
                my ($gid, $type) = ("$name\_UNK", "UNKNOWN");
                if (defined($data{by_tid}{$name})) {
                        $gid  = $data{by_tid}{$name};
                        $type = $data{by_gene}{$gid}{type};
                }
                $gid{name}{$name} = $gid;
                $gid{type}{$type}++;
                $total_type++;

                # Check protein coding gene statistics
                if ($type eq "protein_coding") {
                        my @coor = @{$data{by_gene}{$gid}{tid}{$name}};
                        my ($chr2, $start2, $end2) = split("_", $coor[0]);
                        my $length = $end2 - $start2;
			$stats{2}{chr}{$chr2} ++;
			$total_chr2 ++;
			push(@{$stats{2}{length}}, $length);
                }
        }

	foreach my $gid (keys %{$data{by_type}{protein_coding}}) {
		foreach my $tid (keys %{$data{by_type}{protein_coding}{$gid}}) {
			my @coor = @{$data{by_type}{protein_coding}{$gid}{$tid}};
			my ($chr, $start, $end) = split("_", $coor[0]);
			my $length = $end - $start;
			$stats{1}{chr}{$chr} ++;
			$total_chr1 ++;
			push(@{$stats{1}{length}}, $length);
		}
	}
	# Check Chromosome: Genes set is not enriched at certain chr.
	foreach my $chr (sort keys %{$stats{1}{chr}}) {
		my $chr1 = int($stats{1}{chr}{$chr}/$total_chr1 * 10000) / 100;
		my $chr2 = defined($stats{2}{chr}{$chr}) ? int($stats{2}{chr}{$chr}/$total_chr2 * 10000) / 100 : 0;
		print "$chr\t$chr1\t$chr2\n";
	}
	
	# Check Length: Length of genes are not too high/low
	my @length1 = @{$stats{1}{length}};
        my $mean = mean(@length1);
        my $median = median(@length1);
        my $mode = mode(@length1);
        my $stddev = stddev(@length1);

	my @length2 = @{$stats{2}{length}};
        my $mean2 = mean(@length2);
        my $median2 = median(@length2);
        my $mode2 = mode(@length2);
        my $stddev2 = stddev(@length2);

	my $length1 = R_toolbox::newRArray(\@length1, "length1");
	my $length2 = R_toolbox::newRArray(\@length2, "length2");
	my $Rscript = "
	$length1
	$length2

	pdf(\"CheckLength.pdf\")
	boxplot(log(length1),log(length2),outline=FALSE)
	dev.off()
	";
	R_toolbox::execute_Rscript($Rscript);
        print "mean\t$mean\t$mean2\n";
        print "median\t$median\t$median2\n";
        print "mode\t$mode\t$mode2\n";
        print "stddev\t$stddev\t$stddev2\n";
        close $in;

        #foreach my $type (sort keys %{$gid{type}}) {
        #        my $perc = int($gid{type}{$type} / $total_type * 10000)/100;
        #        print "$type\t$gid{type}{$type} ($perc)\n";
        #}
}


