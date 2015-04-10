#!/usr/bin/perl
# Intersect DRIP/DRIPc peak with all the files here and count enrichment
use strict; use warnings; use mitochy; use Getopt::Std;
use vars qw($opt_a $opt_r);
getopts("ar:");

my ($folder) = @ARGV;
die "usage: $0 [-a to redo everything, -r for result] <folder>\n" unless @ARGV;

#my $DRIP1 = "/data/mitochi/Work/Project/DRIPc/bed/DRIP1_Peak.bed";
my $DRIP1 = "/data/mitochi/Work/Project/DRIPc/bed/DRIP1_Peak_NoPromoter.bed";
my $DRIP2 = "/data/mitochi/Work/Project/DRIPc/bed/DRIP2_Peak_NoPromoter.bed";
my $DRIPc = "/data/mitochi/Work/Project/DRIPc/bed/dripc.bed";
my $DRIP1_shuf = "/data/mitochi/Work/Project/DRIPc/sides/TFBS/Shuffle/DRIP1_Peak_Shuffle.bed";
my $DRIP2_shuf = "/data/mitochi/Work/Project/DRIPc/sides/TFBS/Shuffle/DRIP2_Peak_Shuffle.bed";
my $DRIPc_shuf = "/data/mitochi/Work/Project/DRIPc/sides/TFBS/Shuffle/DRIPc_Peak_Shuffle.bed";
my @drip = ($DRIP1, $DRIP2, $DRIPc, $DRIP1_shuf, $DRIP2_shuf, $DRIPc_shuf);
my (@peak) = <$folder/*.bed>;
my @cell = qw(A549 Ecc1 Gm12878 Gm12891 Gm12892 H1hesc H1neurons Hct116 Helas3 Hepg2 Hl60 Huvec K562 Mcf7 Panc1 Pfsk1 Sknmc Sknsh T47d U87);
my @protocol = qw(V0422111 V0416101 V0416102 Pcr1x Pcr2x);
my %stat;

my $resultfile = defined($opt_r) ? $opt_r : "RESULT.txt";
if (-e "$resultfile") {
	open (my $in, "<", "$resultfile") or die "Cannot read from $resultfile: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		my ($file, $md5, $count1_total, $count2_total, $length_total) = split("\t", $line);
		$stat{$file}{md5} = $md5;
		$stat{$file}{count1_total} = $count1_total;
		$stat{$file}{count2_total} = $count2_total;
		$stat{$file}{length_total} = $length_total;
		#print "$line\n";
	}
}

=comment
foreach my $drip (@drip) {
	print "Getting stats on $drip\n";
	my $check = getStats($drip);
	next if $check eq "NA";
	%{$stat{$drip}} = %{$check};
}
foreach my $peak (@peak) {
	print "Getting stats on $peak\n";
	# Check how many treatment exist, use treatment: 
	my $check = getStats($peak);
	next if $check eq "NA";
	%{$stat{$peak}} = %{$check};
}

my $total_cell = @peak;
my $total_drip = @drip;
my $total_file = $total_cell * $total_drip;
my $filecount = 0;
foreach my $peak (@peak) {
	my ($cell, $tf);
	foreach my $cell_candidate (@cell) {
		$cell = $cell_candidate if $peak =~ /$cell_candidate/i;
	}
	die "Fatal: No Cell Type Specified at file $peak\n" unless defined($cell);
	($tf) = $peak =~ /Pcr\dx/ ? $peak =~ /wgEncodeHaibTfbs$cell(\w+)Pcr\dx/ : $peak =~ /wgEncodeHaibTfbs$cell(\w+)V0/;
	die "Fatal: No TF Type Specified at file $peak\n" unless defined($tf);

	print "Intersecting $cell $tf $peak\n";
	# Intersect with DRIP1
	# Intersect with DRIP2
	# Intersect with DRIPc
	# And shuffles
	foreach my $drip (@drip) {
		$filecount ++;
		my ($fileNamedrip) = mitochy::getFilename($drip, "fullname");
		printf "\t(%.2f %%): Intersecting $cell $tf with $fileNamedrip\n", int(1000*$filecount / $total_file)/10;
		my ($count_peak1, $count_peak2, $length_peak, $md5) = intersect($drip, $peak);
		my $name = "$drip\_$peak";
		next if $count_peak1 eq "NA";
		$stat{$name}{count1_total} = $count_peak1;
		$stat{$name}{count2_total} = $count_peak2;
		$stat{$name}{length_total} = $length_peak;
		$stat{$name}{md5} = $md5;
	}
}
my $time = mitochy::getTime();
open (my $out, ">", "RESULT\_$time.txt") or die;
foreach my $file (sort keys %stat) {
	#my ($file, $md5, $count1_total, $count2_total, $length_total) = split("\t", $line);
	my $md5 = $stat{$file}{md5};
	my $count1_total = $stat{$file}{count1_total};
	my $count2_total = $stat{$file}{count2_total};
	my $length_total = $stat{$file}{length_total};
	print $out "$file\t$md5\t$count1_total\t$count2_total\t$length_total\n";
}
close $out;
system("rm RESULT.txt && ln -s RESULT\_$time.txt RESULT.txt");
=cut
foreach my $peak (@peak) {
	foreach my $drip (@drip) {
		next if $drip =~ /Shuf/;
		my $drip_shuf = $drip =~ /DRIP1/ ? $DRIP1_shuf : $drip =~ /DRIP2/ ? $DRIP2_shuf : $DRIPc_shuf;
		my ($filename) = $drip =~ /\/(\w+)_Peak/; $filename = "DRIPc" if not defined($filename);
		my ($cell, $tf) = getType($peak);
		my $lengthcell = length($cell);
		for (my $i = 0; $i < (13 - $lengthcell); $i++) {
			$cell .= " ";
		}
		my $lengthtf = length($tf);
		for (my $i = 0; $i < (13 - $lengthtf); $i++) {
			$tf .= " ";
		}
		die "DIED AT $tf\n" if length($tf) < 12;
		my $drip_count_total 	  = $stat{$drip}{count1_total};
		my $drip_shuf_count_total = $stat{$drip_shuf}{count1_total};
		my $drip_count_drip 	  = $stat{"$drip\_$peak"}{count1_total};
		die "Died at peak $peak drip $drip\n" if not defined($drip_count_drip);
		my $drip_shuf_count_drip  = $stat{"$drip_shuf\_$peak"}{count1_total};
		my $drip_count_peak 	  = $stat{"$drip\_$peak"}{count2_total};
		my $drip_shuf_count_peak  = $stat{"$drip_shuf\_$peak"}{count2_total};
		my $peak_count_total      = $stat{$peak}{count1_total};

		my $drip_length_total 	   = $stat{$drip}{length1_total};
		my $drip_shuf_length_total = $stat{$drip_shuf}{length1_total};
		my $drip_length 	   = $stat{"$drip\_$peak"}{length_total};
		my $drip_shuf_length       = $stat{"$drip_shuf\_$peak"}{length_total};
		my $peak_length_total      = $stat{$peak}{length1_total};
		
		printf "$cell\t$tf\t$filename\tDRIP\t%.1f\tSHUFDRIP\t%.1f\tPEAK\t%.1f\tSHUFPEAK\t%.1f\n",
100*$drip_count_drip / $drip_count_total,
100*$drip_shuf_count_drip / $drip_shuf_count_total,
100*$drip_count_peak / $peak_count_total,
100*$drip_shuf_count_peak / $peak_count_total



		;
	}
}

sub getType {
	my ($file) = @_;

	my ($cell, $tf);
	foreach my $cell_candidate (@cell) {
		$cell = $cell_candidate if $file =~ /$cell_candidate/i;
	}
	die "Fatal: No Cell Type Specified at file $file\n" unless defined($cell);
	($tf) = $file =~ /Pcr\dx/ ? $file =~ /wgEncodeHaibTfbs$cell(\w+)Pcr\dx/ : $file =~ /wgEncodeHaibTfbs$cell(\w+)V0/;
	die "Fatal: No TF Type Specified at file $file\n" unless defined($tf);
	return($cell, $tf);
}
sub getStats {
	my ($file) = @_;
	my %data;
	my ($md5) = `md5sum $file` =~ /^(.+) $file/;
	return ("NA") if defined($stat{$file}) and $stat{$file}{md5} eq $md5 and not $opt_a;
	my ($count_total)   = `wc -l $file` =~ /^(\d+)/;
	my ($length_total)  = `bedtools_bed_stats.pl $file` =~ /\t(\d+) bp \(/;
	$data{count1_total}  = $count_total;
	$data{count2_total}  = $count_total;
	$data{length_total}  = $length_total;
	$data{md5}	     = $md5;
	return(\%data);
}
sub intersect {
	my ($file1, $file2) = @_;
	my ($folder1, $name1) = mitochy::getFilename($file1, "folder");
	my ($folder2, $name2) = mitochy::getFilename($file2, "folder");
	my ($md51) = `md5sum $file1` =~ /^(.+) $file1/;
	my ($md52) = `md5sum $file2` =~ /^(.+) $file2/;
	my $name = "$file1\_$file2";
	my $md5  = "$md51\_$md52";
	return("NA") if defined($stat{$name}) and $stat{$name}{md5} eq $md5 and not $opt_a;
	my ($count_peak1) = `bedtools intersect -u -a $file1 -b $file2 | wc -l` =~ /^(\d+)/;
	my ($count_peak2) = `bedtools intersect -u -a $file2 -b $file1 | wc -l` =~ /^(\d+)/;
	my ($length_peak) = `bedtools intersect -a $file1 -b $file2 > $name1\_$name2\_TEMP && bedtools_bed_stats.pl $name1\_$name2\_TEMP` =~ /\t(\d+) bp \(/;
	`mv $name1\_$name2\_TEMP /data/mitochi/Recycle_bin/`;
	print "return($count_peak1, $count_peak2, $length_peak, $md5)\n" if $count_peak1 == 0;
	return($count_peak1, $count_peak2, $length_peak, $md5);
}
__END__
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");

open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}
close $in1;

open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;

	my $DRIP1_count_total = $stat{$DRIP1}{count1_total};
	my $DRIP2_count_total = $stat{$DRIP2}{count1_total};
	my $DRIPc_count_total = $stat{$DRIPc}{count1_total};
	my $DRIP1_shuf_count_total = $stat{$DRIP1_shuf}{count1_total};
	my $DRIP2_shuf_count_total = $stat{$DRIP2_shuf}{count1_total};
	my $DRIPc_shuf_count_total = $stat{$DRIPc_shuf}{count1_total};
	my $count_peak_total = $stat{$peak}{count1_total};
	my $DRIP1_count_drip = $stat{"$DRIP1\_$peak"}{count1_total};
	my $DRIP1_count_peak = $stat{"$DRIP1\_$peak"}{count2_total};
	my $DRIP2_count_drip = $stat{"$DRIP2\_$peak"}{count1_total};
	my $DRIP2_count_peak = $stat{"$DRIP2\_$peak"}{count2_total};
	my $DRIPc_count_drip = $stat{"$DRIPc\_$peak"}{count1_total};
	my $DRIPc_count_peak = $stat{"$DRIPc\_$peak"}{count2_total};
	my $DRIP1_shuf_count_drip = $stat{"$DRIP1_shuf\_$peak"}{count1_total};
	my $DRIP1_shuf_count_peak = $stat{"$DRIP1_shuf\_$peak"}{count2_total};
	my $DRIP2_shuf_count_drip = $stat{"$DRIP2_shuf\_$peak"}{count1_total};
	my $DRIP2_shuf_count_peak = $stat{"$DRIP2_shuf\_$peak"}{count2_total};
	my $DRIPc_shuf_count_drip = $stat{"$DRIPc_shuf\_$peak"}{count1_total};
	my $DRIPc_shuf_count_peak = $stat{"$DRIPc_shuf\_$peak"}{count2_total};

