#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($folders) = @ARGV;
die "Usage: $0 <folder containing .shuffled>\n" unless @ARGV;
# Get shuffled regions and turns them into +/- 5kb
my %bad;
open (my $ins, "<", "Bidirectional.bed") or die;
while (my $line = <$ins>) {
	chomp($line);
	$bad{$line} =1;
}
close $ins;
#my @files = </data/mitochi/Work/Project/DRIPc/4_Chromatin/1_Shuffle/Result_125/dripc*.shuffled>;
my @files = <$folders/*.shuffled>;
foreach my $file (@files) {
	print "Processing file $file\n";
	# Promoter/term/antisense need to get fixed (+/-2kb of TSS/TTS) and +/- 5kb
	# Genebody just take 5 shuffles and print out
	# Intergenic just print out 
	#next if $file !~ /antisense\./;
	#next if $file !~ /promoter\./ and $file !~ /antisense\./ and $file !~ /terminal\./;
	#next if $file !~ /NODRIP_terminal/ and $file !~ /NODRIP_antisense/ and $file !~ /NODRIP_genebody/;
	#next if $file !~ /NODRIP_promoter/;
	getCoorFixed($file);
}

sub getCoorFixed {
	my ($input) = @_;
	my ($folder, $fileName) = mitochy::getFilename($input, "folder");
	my ($out1, $out2);
	my $gencode;# = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";
	# CENTERED ON TSS/TTS 
	#$gencode = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_promoter.bed"  if $input =~ /promoter/;
	#$gencode = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_terminal.bed"  if $input =~ /terminal/;
	#$gencode = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_antisense.bed" if $input =~ /antisense/;
	# CENTERED ON DRIPC
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_promoter.bed" if $input =~ /promoter/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_terminal.bed" if $input =~ /terminal/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_antisense.bed" if $input =~ /antisense/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_genebody.bed" if $input =~ /genebody/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_antisense_ext.bed" if $input =~ /antisense_ext/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_antisense_other.bed" if $input =~ /antisense_other/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_promoter_ext.bed" if $input =~ /promoter_ext/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_terminal_ext.bed" if $input =~ /terminal_ext/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_both.bed" if $input =~ /both/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_both_ext.bed" if $input =~ /both_ext/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/dripc_intergenic.bed" if $input =~ /intergenic/;

	if ($gencode =~ /gencode/) {
		my $outname = $fileName;
		($outname) =~ s/dripc/dripcFIX/;
		open ($out1, ">", "$outname\_orig.tmp") or die "Cannot write to $outname\_orig.tmp: $!\n";
		open ($out2, ">", "$outname\_shuf.tmp") or die "Cannot write to $outname\_shuf.tmp: $!\n";

	}
	else {
		open ($out1, ">", "$fileName\_orig.tmp") or die "Cannot write to $fileName\_orig.tmp: $!\n";
		open ($out2, ">", "$fileName\_shuf.tmp") or die "Cannot write to $fileName\_shuf.tmp: $!\n";

	}
	my %data;
	open (my $in1, "<", $gencode) or die "Cannot read from $gencode: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $name, $val, $strand) = split("\t", $line);
		$data{$name} = $line;
		#print "$name $line\n";
	}
	close $in1;

	# Get all lines so we can get random 15 shuffles	
	my %lines;
	open (my $in2, "<", $input) or die "Cannot read from $input: $!\n";
	while (my $line = <$in2>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
		my ($orig, $shuf) = $info =~ /ORIG=(\w+\.?\d*),chr.+TWIN=(\w+\.?\d*),chr/;
		next if defined($bad{$orig}) and $bad{$orig} == 1;
		next if defined($bad{$shuf}) and $bad{$shuf} == 1;
		next if $chr eq "chrM";
		
		push(@{$lines{$name}}, $line);
	}
	close $in2;
	
	my %usedOrig;
	foreach my $name (keys %lines) {
		my @val = @{$lines{$name}};
		my @line = @val;
		for (my $i = 0; $i < @line; $i++) {
			my $line = $line[$i];
			my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
			die "Died at line $line\n" if not defined($info);
			my ($orig, $shuf) = $info =~ /ORIG=(\w+\.?\d*),chr.+TWIN=(\w+\.?\d*),chr/;
			die "Died at getting original/shuffled gene name at file $input: Failed at getCoorFixed\n" if not defined($orig) or not defined($shuf);
	
			if ($gencode =~ /gencode/ and ($input =~ /dripc_promoter\./ or $input =~ /dripc_terminal\./ or $input =~ /dripc_antisense\./)) {
				my $origCoor = $data{$orig}; my ($chrOrig, $startOrig, $endOrig, $nameOrig, $valOrig, $strandOrig) = split("\t", $origCoor);
				my $shufCoor = $data{$shuf}; my ($chrShuf, $startShuf, $endShuf, $nameShuf, $valShuf, $strandShuf) = split("\t", $shufCoor);
				die "1 Died at getting original/shuffled gene coordinates at file $input: Failed at getCoorFixed\n" if not defined($origCoor) or not defined($shufCoor);
				print $out1 "$chrOrig\t$startOrig\t$endOrig\t$orig\t$valOrig\t$strandOrig\t$info\n" if not defined($usedOrig{$orig});
				print $out2 "$chrShuf\t$startShuf\t$endShuf\t$shuf\t$valShuf\t$strand\t$info\n";
				$usedOrig{$orig} = 1;
			}
			else {
				#print "NAME $name\n";
				my $origCoor = $data{$name}; $origCoor = $data{$orig} if not defined($data{$name});
				my ($chrOrig, $startOrig, $endOrig, $nameOrig, $valOrig, $strandOrig) = split("\t", $origCoor);
				die "2 Died at getting original/shuffled gene coordinates at file $input\n$line\n: Failed at getCoorFixed\n" if not defined($origCoor);
				print $out1 "$chrOrig\t$startOrig\t$endOrig\t$orig\t$valOrig\t$strandOrig\t$info\n" if not defined($usedOrig{$name});
				print $out2 "$chr\t$start\t$end\t$shuf\t$val\t$strand\t$info\n";
				$usedOrig{$name} = 1;
			}
		}
	}
	
	runbash("bedtools_bed_change.pl -c -x -5000 -y 5000 -i $fileName\_shuf.tmp -o $fileName\_shuf.input");
	runbash("bedtools_bed_change.pl -c -x -5000 -y 5000 -i $fileName\_orig.tmp -o $fileName\_orig.input");
}

sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}














__END__
my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}

close $in;
close $out;
